#pragma once

#include "LATfield2.hpp"
#include <array>
#include <cmath>

namespace gevolution
{

template<class F_type>
double show_msq(const F_type& F,int i=-1, int j=-1)
{
    double rms =0 ;
    LATfield2::Site x(F.lattice());
    
    if(i<0)
    for(x.first();x.test();x.next())
    {
        double v = F(x);
        rms += v*v;
    }else if(j<0)
    for(x.first();x.test();x.next())
    {
        double v = F(x,i);
        rms += v*v;
    }else
    for(x.first();x.test();x.next())
    {
        double v = F(x,i,j);
        rms += v*v;
    }
    LATfield2::parallel.sum(rms);
    return rms;
}
template<class F_type>
double show_mean(const F_type& F,int i=-1, int j=-1)
{
    double mean =0 ;
    LATfield2::Site x(F.lattice());
    
    if(i<0)
    for(x.first();x.test();x.next())
    {
        double f = F(x);
        mean += f;
    }else if(j<0)
    for(x.first();x.test();x.next())
    {
        double f = F(x,i);
        mean += f;
    }else
    for(x.first();x.test();x.next())
    {
        double f = F(x,i,j);
        mean += f;
    }
    LATfield2::parallel.sum(mean);
    const long N = F.lattice().size(0);
    return mean/N/N/N;
}

template<typename complex_type, typename particle_container>
class particle_mesh
{
    public:
    using real_type = typename complex_type::value_type;
    using particle_type = typename particle_container::value_type;
    
    using real_field_type = LATfield2::Field<real_type>;
    using complex_field_type = LATfield2::Field<complex_type>;
    using site_type = LATfield2::Site;
    using fft_plan_type = LATfield2::PlanFFT<complex_type>;
    
    
    std::size_t my_size;
    LATfield2::Lattice lat,latFT;
    
    std::size_t size() const { return my_size;  }
    
    particle_mesh(int N):
        my_size{N},
        lat(/* dims        = */ 3,
            /* size        = */ N,
            /* ghost cells = */ 2),
        latFT(lat,0,LATfield2::Lattice::FFT::RealToComplex)
    {
    }
    
    const LATfield2::Lattice& lattice() const 
    {
        return lat;
    }
    LATfield2::Lattice& lattice()
    {
        return lat;
    }
    
    void scalar_to_zero(real_field_type& F)
    {
        site_type x(lat);
        for(x.first();x.test();x.next())
            F(x) = 0.0;
        F.updateHalo();
    }
    
    
    std::array<double,3> test_velocities(const particle_container& pcls) const
    {
        double mass{},massvel{},masspos{};
        int count{};
        pcls.for_each(
            [&mass,&massvel,&masspos,&count]
            (const particle_type& part, const site_type& /*xpart*/)
            {
               double v2 = 0,p2=0;
               for(int i=0;i<3;++i)
               {
                   v2 += part.vel[i]*part.vel[i];
                   p2 += part.pos[i]*part.pos[i];
               }
               masspos += p2*part.mass;
               massvel += v2*part.mass;
               mass += part.mass;
               count++;
            }
            );
        LATfield2::parallel.sum(count);
        LATfield2::parallel.sum(masspos);
        LATfield2::parallel.sum(massvel);
        LATfield2::parallel.sum(mass);
        massvel /= mass;
        masspos /= mass;
        mass /= count;
        return {mass,std::sqrt(masspos),std::sqrt(massvel)};
    }
    
    
    std::array<real_type,3> gradient(
        const real_field_type& F, 
        const site_type& x,
        const std::array<real_type,3>& pos)const
    // First order CIC gradient
    {
        const int N = size();
        const real_type dx = 1.0 / N;
        
        std::array<real_type,3> ref_dist{0,0,0};
        for(int i=0;i<3;++i)
            ref_dist[i] = pos[i]/dx - x.coord(i);
            
        std::array<real_type,3> grad{0,0,0};
        for(int i=0;i<3;++i)
        {
            const int j=(i+1)%3,k=(j+1)%3;
            grad[i] += (1. - ref_dist[j]) * (1. - ref_dist[k])
                         * (F (x + i) - F (x));
            grad[i] += ref_dist[j] * (1. - ref_dist[k])
                          * (F (x + i + j) - F (x + j));
            grad[i] += (1. - ref_dist[j]) * ref_dist[k]
                          * (F (x + i + k) - F (x + k));
            grad[i] += ref_dist[j] * ref_dist[k]
                          * (F (x + i + j + k) - F (x + j + k));
        }
        return grad;
    }
    
    virtual std::string report() const = 0;
    virtual double density() const = 0;
    virtual void clear_sources() = 0 ;
    virtual void sample(const particle_container& pcls, double a) = 0;
    virtual void compute_potential(double fourpiG, double a, double Hc,double dt, double Omega) = 0;
    virtual void compute_forces(particle_container& pcls, double fourpiG, double a) const = 0;
    virtual ~particle_mesh(){}
    
    virtual std::array<real_type,3> momentum_to_velocity(
                          const std::array<real_type,3>& momentum,
                          const std::array<real_type,3>& position,
                          const LATfield2::Site& xpart,
                          const real_type a) const = 0;
};

}

