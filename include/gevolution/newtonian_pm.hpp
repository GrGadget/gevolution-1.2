#pragma once

#include "gevolution/config.h"
#include "LATfield2.hpp"
#include "gevolution/real_type.hpp"
#include "gevolution/gevolution.hpp"
#include "gevolution/Particles_gevolution.hpp"

namespace gevolution
{

class newtonian_pm
{
    public:
    using real_field_type = LATfield2::Field<Real>;
    using site_type = LATfield2::Site;
    std::size_t my_size;
    LATfield2::Lattice lat,latFT;
    LATfield2::Field<Real> source,phi;
    LATfield2::Field<Cplx> phi_FT;
    LATfield2::PlanFFT<Cplx> plan_source;
    LATfield2::PlanFFT<Cplx> plan_phi;
    
    public:
    std::size_t size() const { return my_size;  }
    newtonian_pm(int N):
        my_size{N},
        lat(/* dims        = */ 3,
            /* size        = */ N,
            /* ghost cells = */ 2),
        latFT(lat,0,LATfield2::Lattice::FFT::RealToComplex),
        
        source(lat,1),
        phi(lat,1),
        
        phi_FT(latFT,1),
        
        plan_source(&source,&phi_FT),
        plan_phi (&phi, &phi_FT)
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
    
    void clear_sources()
    {
        scalar_to_zero(source);
    }
    
    /*
        sample particle masses into the source field
    */
    void sample(const Particles_gevolution& pcls)
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
    void compute_potential(double factor=1)
    {
        update_kspace();
        solve_poisson_eq(factor);
        update_rspace();
    }
    
    template<class Functor>
    void apply_filter_kspace(Functor f)
    {
        rKSite k(phi_FT.lattice());
        for (k.first(); k.test(); k.next())
        {
            phi_FT(k) *= f({k.coord(0),k.coord(1),k.coord(2)});
        }
        phi_FT.updateHalo();
    }
    template<class Functor>
    void apply_filter_rspace(Functor f)
    {
        Site x(phi.lattice());
        for (x.first(); x.test(); x.next())
        {
            phi(x) *= f({x.coord(0),x.coord(1),x.coord(2)});
        }
        phi.updateHalo();
    }
    
    std::array<Real,3> gradient(
        const real_field_type& F, 
        const LATfield2::Site& x,
        const std::array<Real,3>& pos)const
    // First order CIC gradient
    {
        const int N = size();
        const Real dx = 1.0 / N;
        
        std::array<Real,3> ref_dist{0,0,0};
        for(int i=0;i<3;++i)
            ref_dist[i] = pos[i]/dx - x.coord(i);
            
        std::array<Real,3> grad{0,0,0};
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
    
    /*
        compute forces
        factor = 4 pi G
    */
    void compute_forces(Particles_gevolution& pcls, double factor = 1.0) const
    {
    #ifdef GEVOLUTION_OLD_VERSION
        const double dx = 1.0/size();
        factor /= dx;
        
        LATfield2::Site x(lat);
        LATfield2::Site xpart(pcls.lattice());
        
        for(xpart.first();xpart.test();xpart.next())
        {
            for(auto& part : pcls.field()(xpart).parts )
            {
                std::array<Real,3> pos{part.pos[0],part.pos[1],part.pos[2]};
                std::array<Real,3> gradphi=gradient(phi,xpart,pos);
                for (int i=0;i<3;i++)
                {
                    part.acc[i] = -gradphi[i] * factor;
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
        factor /= dx;
        
        LATfield2::Field<Real> Fx(lat);
        
        LATfield2::Site x(lat);
        LATfield2::Site xpart(pcls.lattice());
        
        // phi.updateHalo();
        for(int i=0;i<3;++i)
        {
            for(x.first();x.test();x.next())
            {
                Fx(x)
                = (-1)*factor*( 
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
    virtual ~newtonian_pm(){}
};

}
