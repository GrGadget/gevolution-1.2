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
    public:
    using real_field_type = LATfield2::Field<Real>;
    using complex_field_type = LATfield2::Field<Cplx>;
    using fft_plan_type = LATfield2::PlanFFT<Cplx>;
    using site_type = LATfield2::Site;
    
    // lattice for the real and transform spaces
    LATfield2::Lattice lat,latFT;
    
    // metric perturbations
    public:
    std::size_t my_size;
    real_field_type phi,chi,Bi;
    real_field_type T00,T0i,Tij;
    
    complex_field_type phi_FT, chi_FT, Bi_FT;
    
    // source fields
    complex_field_type T00_FT, T0i_FT, Tij_FT;
    
    // FT plans
    fft_plan_type plan_phi, plan_chi, plan_Bi;
    fft_plan_type plan_T00, plan_T0i, plan_Tij;
    
    void scalar_to_zero(real_field_type& F)
    {
        site_type x(lat);
        for(x.first();x.test();x.next())
            F(x) = 0.0;
        F.updateHalo();
    }
    void vector_to_zero(real_field_type& F)
    {
        site_type x(lat);
        for(x.first();x.test();x.next())
        {
            for(int i=0;i<3;++i)
                F(x,i) = 0.0;
        }
        F.updateHalo();
    }
    void tensor_to_zero(real_field_type& F)
    {
        site_type x(lat);
        for(x.first();x.test();x.next())
        {
            for(int i=0;i<3;++i)
                for(int j=0;j<3;++j)
                    F(x,i,j) = 0.0;
        }
        F.updateHalo();
    }
    
    public:
    relativistic_pm(int N):
        my_size{N},
        
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
        Bi_FT (latFT,3),
        
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
        scalar_to_zero(phi);
        scalar_to_zero(chi);
        vector_to_zero(Bi);
        
        scalar_to_zero(T00);
        vector_to_zero(T0i);
        tensor_to_zero(Tij);
    }
    
    const LATfield2::Lattice& lattice() const 
    {
        return lat;
    }
    LATfield2::Lattice& lattice()
    {
        return lat;
    }
    
    std::size_t size() const { return my_size;  }
    
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
    
    void compute_phi(
        double a, double Hc, double fourpiG, double dt, double Omega)
    {
        const double dx = 1.0/lat.size()[0];
        prepareFTsource<Real> (
            phi, 
            chi, 
            T00,
            Omega,
            T00, 
            3. * Hc * dx * dx / dt,
            fourpiG * dx * dx / a,
            3. * Hc * Hc * dx * dx); 
        plan_T00.execute (LATfield2::FFT_FORWARD);
        T00_FT.updateHalo ();
        solveModifiedPoissonFT (/* source = */ T00_FT, 
                                /* poten. = */ phi_FT, 
                                1. / (dx * dx),
                                3. * Hc/ dt);
        plan_phi.execute (LATfield2::FFT_BACKWARD); // go back to position space
        phi.updateHalo (); // update ghost cells
    }
    void compute_chi(double f = 1.0)
    {
        prepareFTsource<Real> (
            phi, 
            Tij, 
            Tij,
            2. * f);
        plan_Tij.execute (LATfield2::FFT_FORWARD);
        Tij_FT.updateHalo ();
        projectFTscalar (Tij_FT,chi_FT);
        plan_chi.execute(LATfield2::FFT_BACKWARD);
        chi.updateHalo();
    }
    void compute_Bi(double f = 1.0)
    {
        plan_T0i.execute(LATfield2::FFT_FORWARD);
        T0i_FT.updateHalo();
        
        projectFTvector (T0i_FT, Bi_FT, f);
        
        plan_Bi.execute(LATfield2::FFT_BACKWARD);
        Bi.updateHalo();
    }
    void compute_potential(
        double a, double Hc, double fourpiG, double dt, double Omega)
    // TODO: can we remove all of these dependencies?
    {
        const double dx = 1.0/lat.size()[0];
        compute_phi(a,Hc,fourpiG,dt,Omega);
        compute_chi(fourpiG*dx*dx/a);
        compute_Bi(fourpiG*dx*dx);
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
    */
    void compute_forces(Particles_gevolution& pcls, double a) const
    {
        const int N = size();
        const Real dx = 1.0 / N;
        
        LATfield2::Site xpart(pcls.lattice());
        for(xpart.first();xpart.test();xpart.next())
        {
            for(auto& part : pcls.field()(xpart).parts )
            {
                std::array<double,3> ref_dist;
                for(int l=0;l<3;++l)
                    ref_dist[l] = part.pos[l]/dx - xpart.coord(l);
                    
                Real gradphi[3] = { 0, 0, 0 }, 
                     pgradB[3] = { 0, 0, 0 },
                     gradchi[3] = {0,0,0};
                Real p2 =   part.vel[0] * part.vel[0] 
                          + part.vel[1] * part.vel[1] 
                          + part.vel[2] * part.vel[2];
                Real e2 = std::sqrt(p2 + a*a);

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

                gradchi[0] += (1. - ref_dist[1]) * (1. - ref_dist[2])
                              * (chi (xpart + 0) - chi (xpart));
                gradchi[1] += (1. - ref_dist[0]) * (1. - ref_dist[2])
                              * (chi (xpart + 1) - chi (xpart));
                gradchi[2] += (1. - ref_dist[0]) * (1. - ref_dist[1])
                              * (chi (xpart + 2) - chi (xpart));
                gradchi[0] += ref_dist[1] * (1. - ref_dist[2])
                              * (chi (xpart + 1 + 0) - chi (xpart + 1));
                gradchi[1] += ref_dist[0] * (1. - ref_dist[2])
                              * (chi (xpart + 1 + 0) - chi (xpart + 0));
                gradchi[2] += ref_dist[0] * (1. - ref_dist[1])
                              * (chi (xpart + 2 + 0) - chi (xpart + 0));
                gradchi[0] += (1. - ref_dist[1]) * ref_dist[2]
                              * (chi (xpart + 2 + 0) - chi (xpart + 2));
                gradchi[1] += (1. - ref_dist[0]) * ref_dist[2]
                              * (chi (xpart + 2 + 1) - chi (xpart + 2));
                gradchi[2] += (1. - ref_dist[0]) * ref_dist[1]
                              * (chi (xpart + 2 + 1) - chi (xpart + 1));
                gradchi[0] += ref_dist[1] * ref_dist[2]
                              * (chi (xpart + 2 + 1 + 0) - chi (xpart + 2 + 1));
                gradchi[1] += ref_dist[0] * ref_dist[2]
                              * (chi (xpart + 2 + 1 + 0) - chi (xpart + 2 + 0));
                gradchi[2] += ref_dist[0] * ref_dist[1]
                      * (chi (xpart + 2 + 1 + 0) - chi (xpart + 1 + 0));
        
                pgradB[0] = ((1. - ref_dist[2]) * (Bi (xpart + 0, 1) - Bi (xpart, 1))
                             + ref_dist[2] * (Bi (xpart + 2 + 0, 1) - Bi (xpart + 2, 1)))
                            * part.vel[1];
                pgradB[0] += ((1. - ref_dist[1]) * (Bi (xpart + 0, 2) - Bi (xpart, 2))
                              + ref_dist[1] * (Bi (xpart + 1 + 0, 2) - Bi (xpart + 1, 2)))
                             * part.vel[2];
                pgradB[0] += (1. - ref_dist[1]) * (1. - ref_dist[2])
                             * ((ref_dist[0] - 1.) * Bi (xpart - 0, 0)
                                + (1. - 2. * ref_dist[0]) * Bi (xpart, 0)
                                + ref_dist[0] * Bi (xpart + 0, 0))
                             * part.vel[0];
                pgradB[0] += ref_dist[1] * (1. - ref_dist[2])
                             * ((ref_dist[0] - 1.) * Bi (xpart + 1 - 0, 0)
                                + (1. - 2. * ref_dist[0]) * Bi (xpart + 1, 0)
                                + ref_dist[0] * Bi (xpart + 1 + 0, 0))
                             * part.vel[0];
                pgradB[0] += (1. - ref_dist[1]) * ref_dist[2]
                             * ((ref_dist[0] - 1.) * Bi (xpart + 2 - 0, 0)
                                + (1. - 2. * ref_dist[0]) * Bi (xpart + 2, 0)
                                + ref_dist[0] * Bi (xpart + 2 + 0, 0))
                             * part.vel[0];
                pgradB[0] += ref_dist[1] * ref_dist[2]
                             * ((ref_dist[0] - 1.) * Bi (xpart + 2 + 1 - 0, 0)
                                + (1. - 2. * ref_dist[0]) * Bi (xpart + 2 + 1, 0)
                                + ref_dist[0] * Bi (xpart + 2 + 1 + 0, 0))
                             * part.vel[0];

                pgradB[1] = ((1. - ref_dist[0]) * (Bi (xpart + 1, 2) - Bi (xpart, 2))
                             + ref_dist[0] * (Bi (xpart + 1 + 0, 2) - Bi (xpart + 0, 2)))
                            * part.vel[2];
                pgradB[1] += ((1. - ref_dist[2]) * (Bi (xpart + 1, 0) - Bi (xpart, 0))
                              + ref_dist[2] * (Bi (xpart + 1 + 2, 0) - Bi (xpart + 2, 0)))
                             * part.vel[0];
                pgradB[1] += (1. - ref_dist[0]) * (1. - ref_dist[2])
                             * ((ref_dist[1] - 1.) * Bi (xpart - 1, 1)
                                + (1. - 2. * ref_dist[1]) * Bi (xpart, 1)
                                + ref_dist[1] * Bi (xpart + 1, 1))
                             * part.vel[1];
                pgradB[1] += ref_dist[0] * (1. - ref_dist[2])
                             * ((ref_dist[1] - 1.) * Bi (xpart + 0 - 1, 1)
                                + (1. - 2. * ref_dist[1]) * Bi (xpart + 0, 1)
                                + ref_dist[1] * Bi (xpart + 0 + 1, 1))
                             * part.vel[1];
                pgradB[1] += (1. - ref_dist[0]) * ref_dist[2]
                             * ((ref_dist[1] - 1.) * Bi (xpart + 2 - 1, 1)
                                + (1. - 2. * ref_dist[1]) * Bi (xpart + 2, 1)
                                + ref_dist[1] * Bi (xpart + 2 + 1, 1))
                             * part.vel[1];
                pgradB[1] += ref_dist[0] * ref_dist[2]
                             * ((ref_dist[1] - 1.) * Bi (xpart + 2 + 0 - 1, 1)
                                + (1. - 2. * ref_dist[1]) * Bi (xpart + 2 + 0, 1)
                                + ref_dist[1] * Bi (xpart + 2 + 0 + 1, 1))
                             * part.vel[1];

                pgradB[2] = ((1. - ref_dist[1]) * (Bi (xpart + 2, 0) - Bi (xpart, 0))
                             + ref_dist[1] * (Bi (xpart + 2 + 1, 0) - Bi (xpart + 1, 0)))
                            * part.vel[0];
                pgradB[2] += ((1. - ref_dist[0]) * (Bi (xpart + 2, 1) - Bi (xpart, 1))
                              + ref_dist[0] * (Bi (xpart + 2 + 0, 1) - Bi (xpart + 0, 1)))
                             * part.vel[1];
                pgradB[2] += (1. - ref_dist[0]) * (1. - ref_dist[1])
                             * ((ref_dist[2] - 1.) * Bi (xpart - 2, 2)
                                + (1. - 2. * ref_dist[2]) * Bi (xpart, 2)
                                + ref_dist[2] * Bi (xpart + 2, 2))
                             * part.vel[2];
                pgradB[2] += ref_dist[0] * (1. - ref_dist[1])
                             * ((ref_dist[2] - 1.) * Bi (xpart + 0 - 2, 2)
                                + (1. - 2. * ref_dist[2]) * Bi (xpart + 0, 2)
                                + ref_dist[2] * Bi (xpart + 0 + 2, 2))
                             * part.vel[2];
                pgradB[2] += (1. - ref_dist[0]) * ref_dist[1]
                             * ((ref_dist[2] - 1.) * Bi (xpart + 1 - 2, 2)
                                + (1. - 2. * ref_dist[2]) * Bi (xpart + 1, 2)
                                + ref_dist[2] * Bi (xpart + 2 + 1, 2))
                             * part.vel[2];
                pgradB[2] += ref_dist[0] * ref_dist[1]
                             * ((ref_dist[2] - 1.) * Bi (xpart + 1 + 0 - 2, 2)
                                + (1. - 2. * ref_dist[2]) * Bi (xpart + 1 + 0, 2)
                                + ref_dist[2] * Bi (xpart + 1 + 0 + 2, 2))
                             * part.vel[2];
                
                for(int i=0;i<3;++i)
                {
                    part.acc[i] = (gradchi[i]*e2 - gradphi[i] * (e2 + p2/e2) 
                                    - pgradB[i]/a/a/N)/dx;
                }
            }
        }
    }
   
    std::array<Real,3> vector_at(
        const real_field_type& F, 
        const LATfield2::Site& x,
        const std::array<Real,3>& pos)const
    {
        const int N = size();
        const Real dx = 1.0 / N;
        
        std::array<Real,3> ref_dist{0,0,0};
        for(int i=0;i<3;++i)
            ref_dist[i] = pos[i]/dx - x.coord(i);
            
        std::array<Real,3> xF{0,0,0};
        
        // CIC interpolation
        for(int i=0;i<3;++i)
        {
            xF[i] = F(x,i) * (1.-ref_dist[0]) * (1.-ref_dist[1]) * (1.-ref_dist[2]);
            xF[i] += F(x+0,i) * ref_dist[0] * (1.-ref_dist[1]) * (1.-ref_dist[2]);
            xF[i] += F(x+1,i) * (1.-ref_dist[0]) * ref_dist[1] * (1.-ref_dist[2]);
            xF[i] += F(x+0+1,i) * ref_dist[0] * ref_dist[1] * (1.-ref_dist[2]);
            xF[i] += F(x+2,i) * (1.-ref_dist[0]) * (1.-ref_dist[1]) * ref_dist[2];
            xF[i] += F(x+0+2,i) * ref_dist[0] * (1.-ref_dist[1]) * ref_dist[2];
            xF[i] += F(x+1+2,i) * (1.-ref_dist[0]) * ref_dist[1] * ref_dist[2];
            xF[i] += F(x+0+1+2,i) * ref_dist[0] * ref_dist[1] * ref_dist[2];
        }
        return xF;
    }
    Real scalar_at(
        const real_field_type& F, 
        const LATfield2::Site& x,
        const std::array<Real,3>& pos)const
    {
        const int N = size();
        const Real dx = 1.0 / N;
        
        std::array<Real,3> ref_dist{0,0,0};
        for(int i=0;i<3;++i)
            ref_dist[i] = pos[i]/dx - x.coord(i);
            
        Real xF{0};
        
        // CIC interpolation
        xF = F(x) * (1.-ref_dist[0]) * (1.-ref_dist[1]) * (1.-ref_dist[2]);
        xF += F(x+0) * ref_dist[0] * (1.-ref_dist[1]) * (1.-ref_dist[2]);
        xF += F(x+1) * (1.-ref_dist[0]) * ref_dist[1] * (1.-ref_dist[2]);
        xF += F(x+0+1) * ref_dist[0] * ref_dist[1] * (1.-ref_dist[2]);
        xF += F(x+2) * (1.-ref_dist[0]) * (1.-ref_dist[1]) * ref_dist[2];
        xF += F(x+0+2) * ref_dist[0] * (1.-ref_dist[1]) * ref_dist[2];
        xF += F(x+1+2) * (1.-ref_dist[0]) * ref_dist[1] * ref_dist[2];
        xF += F(x+0+1+2) * ref_dist[0] * ref_dist[1] * ref_dist[2];
        
        return xF;
    }
   
    std::array<Real,3> momentum_to_velocity(
                          const std::array<Real,3>& momentum,
                          const std::array<Real,3>& position,
                          const LATfield2::Site& xpart,
                          const Real a) const
    {
        const int N = size();
        std::array<Real,3> velocity{0,0,0};
        Real xchi{0},xphi{0};
        std::array<Real,3> xBi{0,0,0};
        const Real momentum2 =    momentum[0]*momentum[0] 
                                + momentum[1]*momentum[1]
                                + momentum[2]*momentum[2];
        const Real e2 = std::sqrt(momentum2 + a*a);
                
        
        xphi = scalar_at(phi,xpart,{position[0],position[1],position[2]});
        xchi = scalar_at(chi,xpart,{position[0],position[1],position[2]});
        xBi  = vector_at(Bi,xpart,{position[0],position[1],position[2]});
        
        const Real velocity2 = 
            (1. + (3. - momentum2 / e2/e2) * xphi - xchi) / e2;
        
        for(int i=0;i<3;++i)
        {
            velocity[i] = momentum[i] * velocity2 + xBi[i]/a/a/N; // TODO: why this strange scaling of Bi?
        }
        
        return velocity;
    }
    
    virtual ~relativistic_pm(){}
};
} // namespace gevolution
