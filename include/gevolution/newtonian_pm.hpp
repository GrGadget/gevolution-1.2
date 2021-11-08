#pragma once

#include "gevolution/config.h"
#include "gevolution/particle_mesh.hpp"
#include "LATfield2.hpp"
#include "gevolution/gevolution.hpp"

namespace gevolution
{

template<typename complex_type, typename particle_container>
class newtonian_pm : public particle_mesh<complex_type,particle_container>
{
    public:
    using base_type = particle_mesh<complex_type,particle_container>;
    using typename base_type::real_type;
    using typename base_type::real_field_type;
    using typename base_type::complex_field_type;
    using typename base_type::site_type;
    using base_type::size;
    using base_type::scalar_to_zero;
    using base_type::gradient;
    using typename base_type::fft_plan_type;
    
    
    real_field_type source,phi;
    complex_field_type phi_FT;
    fft_plan_type plan_source, plan_phi;
    
    public:
    newtonian_pm(int N):
        base_type(N),
        
        source(base_type::lat,1),
        phi(base_type::lat,1),
        
        phi_FT(base_type::latFT,1),
        
        plan_source(&source,&phi_FT),
        plan_phi (&phi, &phi_FT)
    {
        scalar_to_zero(source);
        scalar_to_zero(phi);
    }
    
    void clear_sources() override
    {
        scalar_to_zero(source);
    }
    
    /*
        sample particle masses into the source field
    */
    void sample(const particle_container& pcls, double /* a */=0) override
    {
        scalarProjectionCIC_project (&pcls, &source); // samples
        scalarProjectionCIC_comm (&source); // communicates the ghost cells
    }
    
    /*
        computes the potential
    */
    void update_kspace()
    {
        plan_source.execute (LATfield2::FFT_FORWARD); // Newton: directly go to k-space
        phi_FT.updateHalo (); // update ghost cells
    }
    void update_rspace()
    {
        plan_phi.execute (LATfield2::FFT_BACKWARD); // go back to position space
        phi.updateHalo (); // update ghost cells
    }
    void solve_poisson_eq(double factor=1)
    {
        solveModifiedPoissonFT (phi_FT, phi_FT,factor); // Newton: in k-space
        // (4 pi G)/a = 1
    }
    void compute_potential(
        double fourpiG    =1, 
        double a          =1, 
        double /* Hc */   =0, 
        double /* dt */   =0, 
        double /* Omega*/ =0) override
    {
        update_kspace();
        solve_poisson_eq(fourpiG/a);
        update_rspace();
    }
    
    /*
        compute forces
    */
    void compute_forces(
        particle_container& pcls, 
        double fourpiG =1, 
        double a = 1) 
        const override
    {
    #ifdef GEVOLUTION_OLD_VERSION
        const double dx = 1.0/size();
        fourpiG /= dx;
        
        site_type xpart(pcls.lattice());
        
        for(xpart.first();xpart.test();xpart.next())
        {
            for(auto& part : pcls.field()(xpart).parts )
            {
                std::array<real_type,3> pos{part.pos[0],part.pos[1],part.pos[2]};
                std::array<real_type,3> gradphi=gradient(phi,xpart,pos);
                for (int i=0;i<3;i++)
                {
                    part.acc[i] = -gradphi[i] * fourpiG * a;
                }
            }
        }
    #else
        /*
        Let's do like in Gadget4:
        1. compute Fx field from phi at 4th order FD
        2. interpolate Fx at particle's position using CIC
        */
        const double dx = 1.0/pcls.lattice().size()[0];
        fourpiG /= dx;
        
        real_field_type Fx(base_type::lat);
        
        site_type x(base_type::lat);
        site_type xpart(pcls.lattice());
        
        // phi.updateHalo();
        for(int i=0;i<3;++i)
        {
            for(x.first();x.test();x.next())
            {
                Fx(x)
                = (-1)*a*fourpiG*( 
                        2.0/3 * (phi(x+i) - phi(x-i)) 
                        - 1.0/12 * (phi(x+i+i) - phi(x-i-i))  );
            }
            Fx.updateHalo();
            for(xpart.first();xpart.test();xpart.next())
            {
                for(auto& part : pcls.field()(xpart).parts )
                {
                    std::array<double,3> ref_dist;
                    for(int l=0;l<3;++l)
                        ref_dist[l] = part.pos[l]/dx - xpart.coord(l);
                    
                    part.acc[i] = 0.0;
                    
                    part.acc[i] +=
                    (1-ref_dist[0])*(1-ref_dist[1])*(1-ref_dist[2])*Fx(xpart);
                    
                    part.acc[i] +=
                    (ref_dist[0])*(1-ref_dist[1])*(1-ref_dist[2])*Fx(xpart+0);
                    
                    part.acc[i] +=
                    (1-ref_dist[0])*(ref_dist[1])*(1-ref_dist[2])*Fx(xpart+1);
                    
                    part.acc[i] +=
                    (ref_dist[0])*(ref_dist[1])*(1-ref_dist[2])*Fx(xpart+1+0);
                    
                    part.acc[i] +=
                    (1-ref_dist[0])*(1-ref_dist[1])*(ref_dist[2])*Fx(xpart+2);
                    
                    part.acc[i] +=
                    (ref_dist[0])*(1-ref_dist[1])*(ref_dist[2])*Fx(xpart+2+0);
                    
                    part.acc[i] +=
                    (1-ref_dist[0])*(ref_dist[1])*(ref_dist[2])*Fx(xpart+2+1);
                    
                    part.acc[i] +=
                    (ref_dist[0])*(ref_dist[1])*(ref_dist[2])*Fx(xpart+2+1+0);
                }
            }
        }
    #endif
    } 
    ~newtonian_pm() override {}
    
    std::string report() const override
    {
        std::stringstream ss;
        // sources
        double std_t00 = show_msq(source);
        // potentials
        double std_phi = show_msq(phi);
        ss << "RMS(T00) = " << std_t00 << '\n';
        ss << "RMS(Phi) = " << std_phi << '\n';
        return ss.str();
    }
    double density() const override
    {
        return show_mean(source);
    }
    double sum_phi() const override
    {
        const double N = size();
        return show_mean(phi) * N * N *N;
    }
    std::array<real_type,3> momentum_to_velocity(
                          const std::array<real_type,3>& momentum,
                          const std::array<real_type,3>& /* position */,
                          const site_type& /* xpart */,
                          const real_type a) const override
    {
        std::array<real_type,3> velocity{0,0,0};
        
        const real_type inv_a =  1/a;
        
        for(int i=0;i<3;++i)
        {
            velocity[i] = momentum[i] * inv_a ;
        }
        
        return velocity;
    }
    std::array<real_type,3> velocity_to_momentum(
                          const std::array<real_type,3>& velocity,
                          const std::array<real_type,3>& /* position */,
                          const site_type& /* xpart */,
                          const real_type a) const override
    {
        std::array<real_type,3> momentum{0,0,0};
        for(int i=0;i<3;++i)
        {
            momentum[i] = velocity[i] * a ;
        }
        return momentum;
    }
};

template<class functor_type, typename complex_type, typename particle_container>
void apply_filter_kspace(
    newtonian_pm<complex_type,particle_container> &pm,
    functor_type f)
{
    apply_filter_kspace(pm.phi,pm.phi_FT,pm.plan_phi,f);
}

}
