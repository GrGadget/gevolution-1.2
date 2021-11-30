#pragma once

#include "LATfield2.hpp"
#include <array>
#include <cmath>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>

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
template<class T, class field_type, class bin_op_type,class un_op_type>
T reduce_field(
    const ::boost::mpi::communicator& com,
    const field_type& F,
    T initial,
    bin_op_type op1,
    un_op_type op2)
{
    LATfield2::Site x(F.lattice());
    for(x.first();x.test();x.next())
    {
        T f = op2(F(x));
        initial = op1(initial,f);
    }
    ::boost::mpi::all_reduce(com,initial,op1);
    return initial;
}
template<class T, class field_type, class bin_op_type, class un_op_type>
T reduce_field_halo(
    const ::boost::mpi::communicator& com,
    const field_type& F,
    T initial,
    bin_op_type op1,
    un_op_type op2)
{
    LATfield2::Site x(F.lattice());
    for(x.haloFirst();x.haloTest();x.haloNext())
    {
        T f = op2(F(x));
        initial = op1(initial,f);
    }
    ::boost::mpi::all_reduce(com,initial,op1);
    return initial;
}

// TODO
// we need these because the current boost::mpi::all_reduce template does not
// compile with lambdas
template<typename T>
struct my_max_func
{
    T operator()(const T& x, const T& y)
    {
        using std::max;
        return max(x,y);
    }
};
template<typename T>
struct my_sum_func
{
    T operator()(const T& x, const T& y)
    {
        return x+y;
    }
};
template<typename T, typename T2 = T>
struct my_abs_func
{
    T operator()(const T2& x)
    {
        using std::abs;
        return abs(x);
    }
};
template<typename T>
struct my_sqr_func
{
    T operator()(const T& x)
    {
        return x*x;
    }
};
    
// TODO: simplify this template
template<class functor_type, class field_type1, class field_type2, class fft_plan_type >
void apply_filter_kspace(
    field_type1 &phi, 
    field_type2 &phi_FT,
    fft_plan_type &plan,
    functor_type f)
{
    plan.execute(::LATfield2::FFT_FORWARD);
    ::LATfield2::rKSite k(phi_FT.lattice());
    const double N = phi.lattice().size(0);
    const double inv_N3 = 1.0/N/N/N;
    for (k.first(); k.test(); k.next())
    {
        phi_FT(k) *= f({k.coord(0),k.coord(1),k.coord(2)}) * inv_N3;
    }
    phi_FT.updateHalo();
    plan.execute(LATfield2::FFT_BACKWARD);
    phi.updateHalo();
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
    
    
    ::boost::mpi::communicator com;
    std::size_t my_size;
    LATfield2::Lattice lat,latFT;
    
    std::size_t size() const { return my_size;  }
    
    particle_mesh(unsigned int N,const MPI_Comm& that_com):
        com(that_com,::boost::mpi::comm_create_kind::comm_duplicate),
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
    
    
    std::array<double,4> test_velocities(const particle_container& pcls) const
    {
        double mass{},massvel{},masspos{},massacc{};
        int count{};
        pcls.for_each(
            [&mass,&massvel,&masspos,&massacc,&count]
            (const particle_type& part, const site_type& /*xpart*/)
            {
               double v2 = 0,p2=0,a2=0;
               for(int i=0;i<3;++i)
               {
                   v2 += part.momentum[i]*part.momentum[i];
                   p2 += part.pos[i]*part.pos[i];
                   a2 += part.force[i]*part.force[i];
               }
               masspos += p2*part.mass;
               massvel += v2*part.mass;
               massacc += a2*part.mass;
               mass += part.mass;
               count++;
            }
            );
        LATfield2::parallel.sum(count);
        LATfield2::parallel.sum(masspos);
        LATfield2::parallel.sum(massvel);
        LATfield2::parallel.sum(massacc);
        LATfield2::parallel.sum(mass);
        massvel /= mass;
        masspos /= mass;
        mass /= count;
        massacc /= mass;
        using std::sqrt;
        return {mass,sqrt(masspos),sqrt(massvel),sqrt(massacc)};
    }
    
    
    std::array<real_type,3> gradient(
        const real_field_type& F, 
        const site_type& x,
        const std::array<real_type,3>& pos)const
    // First order CIC gradient
    // precondition: F has valid ghost cells
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
    
    virtual std::string report(const particle_container& pcls, double a) const = 0;
    virtual double density() const = 0;
    virtual double sum_phi() const = 0;
    virtual void clear_sources() = 0 ;
    virtual void sample(const particle_container& pcls, double a) = 0;
    virtual void compute_potential(double fourpiG, double a, double Hc,double Omega) = 0;
    
    enum class force_reduction {
        assign, plus, minus
    };
    
    virtual void compute_forces(
        particle_container& pcls, double fourpiG, double a, 
        force_reduction = force_reduction::assign) const = 0;
    virtual ~particle_mesh(){}
    
    virtual void save_to_file( std::string  ) const = 0;
    
    virtual std::array<real_type,3> momentum_to_velocity(
                          const std::array<real_type,3>& momentum,
                          const std::array<real_type,3>& position,
                          const LATfield2::Site& xpart,
                          const real_type a) const = 0;
    virtual std::array<real_type,3> velocity_to_momentum(
                          const std::array<real_type,3>& velocity,
                          const std::array<real_type,3>& position,
                          const LATfield2::Site& xpart,
                          const real_type a) const = 0;
    
    void compute_velocities(particle_container& pcls, double a=1)
    {
        pcls.for_each(
            [&](particle_type &part, const site_type& xpart)
            {
                std::array<real_type,3> velocity =
                    momentum_to_velocity(
                        part.momentum,
                        {part.pos[0],part.pos[1],part.pos[2]},
                        xpart,
                        a);
                for(int i=0;i<3;++i)
                {
                    part.vel[i] = velocity[i];
                }
            }
        );
    }
};

}

