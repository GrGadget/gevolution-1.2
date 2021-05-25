#pragma once

#include "LATfield2.hpp"
#include "Particles_gevolution.hpp"

namespace gevolution
{

class newtonian_pm
{
    using Real = LATfield2::Real;
    using Cplx = LATfield2::Imag;
    
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
        latFT(lat,0),
        
        source(lat,1),
        phi(lat,1),
        
        scalarFT(latFT,1),
        
        plan_source(&source,&scalarFT),
        plan_phi (&phi, &scalarFT)
    {
        // latFT.initializeRealFFT (lat, 0);
        // source.initialize (lat, 1);
        // phi.initialize (lat, 1);
        // scalarFT.initialize (latFT, 1);
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
        solveModifiedPoissonFT (scalarFT, scalarFT,1.0);
                            //    cosmo.fourpiG
                            //        / a); // Newton: phi update (k-space)
        plan_phi.execute (LATfield2::FFT_BACKWARD); // go back to position space
    }
};

}
