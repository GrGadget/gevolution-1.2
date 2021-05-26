#pragma once

#include "LATfield2.hpp"
#include "gevolution.hpp"
#include "Particles_gevolution.hpp"

namespace gevolution
{

class newtonian_pm
{
    LATfield2::Lattice lat,latFT;
    LATfield2::Field<Real> source,phi;
    LATfield2::Field<Cplx> scalarFT;
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
        
        scalarFT(latFT,1),
        
        plan_source(&source,&scalarFT),
        plan_phi (&phi, &scalarFT)
    {
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
    void compute_potential()
    {
        plan_source.execute (LATfield2::FFT_FORWARD); // Newton: directly go to k-space
        solveModifiedPoissonFT (scalarFT, scalarFT,1.0); // Newton: in k-space
        // (4 pi G)/a = 1
        plan_phi.execute (LATfield2::FFT_BACKWARD); // go back to position space
        phi.updateHalo (); // update ghost cells
    }
    
    /*
        compute forces
    */
    void compute_forces(Particles_gevolution& pcls)
    {
        LATfield2::Site xpart(pcls.lattice());
        const double dx = 1.0/pcls.lattice().size()[0];
        for(xpart.first();xpart.test();xpart.next())
        {
            for(auto part : pcls.field()(xpart).parts )
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
                std::array<double,3> acc;
                for (int i = 0; i < 3; i++)
                    acc[i] =  (-1) * gradphi[i] / dx;
            }
        }
    }
};

}
