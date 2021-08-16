#pragma once

#include "gevolution/config.h"
#include "LATfield2.hpp"
#include "gevolution/real_type.hpp"
#include "gevolution/gevolution.hpp"
#include "gevolution/Particles_gevolution.hpp"

/*
    TODO:
    in a first approximation we will only use T00, hence the scalar fields to
    compute the particle dynamics. Later on we will add T0i and Tij.
*/

namespace gevolution
{
class relativistic_pm
{
    using real_field_type = LATfield2::Field<Real>;
    using complex_field_type = LATfield2::Field<Cplx>;
    using fft_plan_type = LATfield2::PlanFFT<Cplx>;
    
    // lattice for the real and transform spaces
    LATfield2::Lattice lat,latFT;
    
    // metric perturbations
    real_field_type phi,chi,Bi;
    complex_field_type phi_FT, chi_FT, Bi_FT;
    
    // source fields
    real_field_type T00,T0i,Tij;
    complex_field_type T00_FT, T0i_FT, Tij_FT;
    
    // FT plans
    fft_plan_type plan_phi, plan_chi, plan_Bi;
    fft_plan_type plan_T00, plan_T0i, plan_Tij;
    
    
    public:
    newtonian_pm(int N):
        lat(/* dims        = */ 3,
            /* size        = */ N,
            /* ghost cells = */ 2),
        latFT(lat,0,LATfield2::Lattice::FFT::RealToComplex),
        
        // initialize fields, metric
        phi(lat,1),
        chi(lat,1),
        Bi (lat,3),
        
        // initialize fields, sources
        T00(lat,1),
        T0i(lat,3),
        Tij(lat,3,3,LATfield2::matrix_symmetry::symmetric),
        
        // initialize k-fields, metric
        phi_FT(latFT,1),
        chi_FT(latFT,1),
        Bi_FT (latFT,1),
        
        // initialize k-fields, sources
        T00_FT(latFT,1),
        T0i_FT(latFT,3),
        Tij_FT(latFT,3,3,LATfield2::matrix_symmetry::symmetric),
        
        // initialize plans
        plan_phi(&phi, &phi_FT),
        plan_chi(&chi, &chi_FT),
        plan_Bi (&Bi, &Bi_FT),
        
        plan_T00(&T00,&T00_FT),
        plan_T0i(&T0i,&T0i_FT),
        plan_Tij(&Tij,&Tij_FT)
    {
        // TODO set phi,chi,Bi to zero
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
    void sample(const Particles_gevolution& pcls, double a) 
    // TODO: we don't actually need the scale factor here if we use particles'
    // canonical momentum normalized q = p/mca.
    {
        // WARNING: has phi been initialized? 
        
        projection_init (&T00); // sets to zero the field
        projection_T00_project(&pcls, &T00, a, &phi); // samples
        projection_T00_comm (&T00); // communicates the ghost cells
        
        projection_init(&T0i);
        projection_T0i_project(&pcls,&T0i,&phi);
        projection_T0i_comm(&T0i);
        
        projection_init(&Tij);
        projection_Tij_project(&pcls,&Tij,a,&phi);
        projection_Tij_comm(&Tij);
    }
    
    /*
        computes the potential
    */
    void update_kspace()
    {
    }
    void update_rspace()
    {
    }
    void solve_poisson_eq()
    {
    }
    
    void compute_phi()
    {
        prepareFTsource<Real> (
            phi, 
            chi, 
            T00,
            cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm (a, cosmo),
            T00, 
            3. * Hconf (a, cosmo) * dx * dx / dtau_old,
            cosmo.fourpiG * dx * dx / a,
            3. * Hconf (a, cosmo) * Hconf (a, cosmo) * dx * dx); 
        
        plan_T00.execute (LATfield2::FFT_FORWARD);
        T00_FT.updateHalo ();
        
        solveModifiedPoissonFT (/* source = */ T00_FT, 
                                /* poten. = */ phi_FT, 
                                1. / (dx * dx),
                                3. * Hconf (a, cosmo)/ dtau_old);
                                
        plan_phi.execute (LATfield2::FFT_BACKWARD); // go back to position space
        phi.updateHalo (); // update ghost cells
    }
    void compute_chi()
    {
        prepareFTsource<Real> (
            phi, 
            Tij, 
            Tij,
            2. * cosmo.fourpiG * dx * dx / a);
        
        plan_Tij.execute (LATfield2::FFT_FORWARD);
        Tij_FT.updateHalo ();
        
        projectFTscalar (Tij_FT,chi_FT);
        
        plan_chi.execute(LATfield2::FFT_BACKWARD);
        chi.updateHalo();
    }
    void compute_Bi()
    {
        plan_Ti0.execute(LATfield2::FFT_FORWARD);
        Ti0_FT.updateHalo();
        
        projectFTvector (T0i_FT, Bi_FT, cosmo.fourpiG * dx * dx);
        
        plan_Bi.execute(LATfield2::FFT_BACKWARD);
        Bi.updateHalo();
    }
    
    void compute_potential()
    {
        compute_phi();
        compute_chi();
        compute_Bi();
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
    
    /*
        compute forces
        factor = 4 pi G
    */
    void compute_forces(Particles_gevolution& pcls, double factor = 1.0) const
    {
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
    }
    
    virtual ~newtonian_pm(){}
};
} // namespace gevolution
