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
    LATfield2::Lattice lat,latFT;
    LATfield2::Field<Real> source,phi;
    LATfield2::Field<Cplx> phi_FT;
    LATfield2::PlanFFT<Cplx> plan_source;
    LATfield2::PlanFFT<Cplx> plan_phi;
    
    public:
    newtonian_pm(int N):
        lat(/* dims        = */ 3,
            /* size        = */ N,
            /* ghost cells = */ 1),
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
    
    /*
        sample particle masses into the source field
    */
    void sample(const Particles_gevolution& pcls)
    {
        projection_init (&source); // sets to zero the field
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
    void solve_poisson_eq()
    {
        solveModifiedPoissonFT (phi_FT, phi_FT,1.0); // Newton: in k-space
        // (4 pi G)/a = 1
    }
    void compute_potential()
    {
        update_kspace();
        solve_poisson_eq();
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
    }
    template<class Functor>
    void apply_filter_rspace(Functor f)
    {
        Site x(phi.lattice());
        for (x.first(); x.test(); x.next())
        {
            phi(x) *= f({x.coord(0),x.coord(1),x.coord(2)});
        }
    }
    
    /*
        compute forces
        factor = 4 pi G
    */
    void compute_forces(Particles_gevolution& pcls, double factor = 1.0)const
    {
        /*
        Let's do like in Gadget4:
        1. compute Fx field from phi at 4th order FD
        2. interpolate Fx at particle's position using CIC
        3. repeat for y,z
        */
        LATfield2::Site xpart(pcls.lattice());
        const double dx = 1.0/pcls.lattice().size()[0];
        factor /= dx;
        for(xpart.first();xpart.test();xpart.next())
        {
            for(auto& part : pcls.field()(xpart).parts )
            {
                std::array<double,3> ref_dist;
                for(int l=0;l<3;++l)
                    ref_dist[l] = part.pos[l]/dx - xpart.coord(l);
                
                std::array<double,3> gradphi{ 0, 0, 0 };
                gradphi[0] = (1. - ref_dist[1]) * (1. - ref_dist[2])
                             * (phi (xpart + 0) - phi (xpart));
                gradphi[1] = (1. - ref_dist[0]) * (1. - ref_dist[2])
                             * (phi (xpart + 1) - phi (xpart));
                gradphi[2] = (1. - ref_dist[0]) * (1. - ref_dist[1])
                             * (phi (xpart + 2) - phi (xpart));
                gradphi[0] += ref_dist[1] * (1. - ref_dist[2])
                              * (phi (xpart + 1 + 0) - phi (xpart + 1));
                gradphi[1] += ref_dist[0] * (1. - ref_dist[2])
                              * (phi (xpart + 1 + 0) - phi (xpart + 0));
                gradphi[2] += ref_dist[0] * (1. - ref_dist[1])
                              * (phi (xpart + 2 + 0) - phi (xpart + 0));
                gradphi[0] += (1. - ref_dist[1]) * ref_dist[2]
                              * (phi (xpart + 2 + 0) - phi (xpart + 2));
                gradphi[1] += (1. - ref_dist[0]) * ref_dist[2]
                              * (phi (xpart + 2 + 1) - phi (xpart + 2));
                gradphi[2] += (1. - ref_dist[0]) * ref_dist[1]
                              * (phi (xpart + 2 + 1) - phi (xpart + 1));
                gradphi[0] += ref_dist[1] * ref_dist[2]
                              * (phi (xpart + 2 + 1 + 0) - phi (xpart + 2 + 1));
                gradphi[1] += ref_dist[0] * ref_dist[2]
                              * (phi (xpart + 2 + 1 + 0) - phi (xpart + 2 + 0));
                gradphi[2] += ref_dist[0] * ref_dist[1]
                              * (phi (xpart + 2 + 1 + 0) - phi (xpart + 1 + 0));
                for(int i=0;i<3;++i)
                    part.acc[i] =  (-1) * gradphi[i] * factor;
                    
               // if(part.ID==1)
               //     std::cout << "from PM Part ID 1: " 
               //     << " acc[0] "<< part.acc[0] 
               //     << "\n";
            }
        }
    }
    
    virtual ~newtonian_pm(){}
};

}
