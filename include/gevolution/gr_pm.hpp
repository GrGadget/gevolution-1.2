#pragma once

#include "gevolution/config.h"
#include "gevolution/particle_mesh.hpp"
#include "LATfield2.hpp"
#include "gevolution/gevolution.hpp"

/*
    TODO:
    in a first approximation we will only use T00, hence the scalar fields to
    compute the particle dynamics. Later on we will add T0i and Tij.
*/

namespace gevolution
{
template<typename complex_type, typename particle_container>
class relativistic_pm : public particle_mesh<complex_type,particle_container>
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
    
    // metric perturbations
    public:
    
    // potentials
    real_field_type phi,chi,Bi;
    complex_field_type phi_FT, chi_FT, Bi_FT;
    
    // Energy-Momentum tensor
    real_field_type T00,T0i,Tij;
    // complex_field_type T00_FT, T0i_FT, Tij_FT;
    
    // Sources
    mutable real_field_type S00, S0i, Sij; // source fields
    mutable complex_field_type S00_FT, S0i_FT, Sij_FT; // source fields
    
    
    // FT plans
    fft_plan_type plan_phi, plan_chi, plan_Bi;
    // fft_plan_type plan_T00, plan_T0i, plan_Tij;
    fft_plan_type plan_S00, plan_S0i, plan_Sij;
    
    void scalar_to_zero(real_field_type& F)
    {
        site_type x(base_type::lat);
        for(x.first();x.test();x.next())
            F(x) = 0.0;
        F.updateHalo();
    }
    void vector_to_zero(real_field_type& F)
    {
        site_type x(base_type::lat);
        for(x.first();x.test();x.next())
        {
            for(int i=0;i<3;++i)
                F(x,i) = 0.0;
        }
        F.updateHalo();
    }
    void tensor_to_zero(real_field_type& F)
    {
        site_type x(base_type::lat);
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
        base_type(N),
        
        // initialize fields, metric
        phi(base_type::lat,1),
        chi(base_type::lat,1),
        Bi (base_type::lat,3),
        
        // initialize k-fields, metric
        phi_FT(base_type::latFT,1),
        chi_FT(base_type::latFT,1),
        Bi_FT (base_type::latFT,3),
        
        // initialize fields, Energy-momentum tensor
        T00(base_type::lat,1),
        T0i(base_type::lat,3),
        Tij(base_type::lat,3,3,LATfield2::matrix_symmetry::symmetric),
        
        // initialize k-fields, Energy-momentum tensor
        // T00_FT(base_type::latFT,1),
        // T0i_FT(base_type::latFT,3),
        // Tij_FT(base_type::latFT,3,3,LATfield2::matrix_symmetry::symmetric),
        
        // initialize fields, sources
        S00(base_type::lat,1),
        S0i(base_type::lat,3),
        Sij(base_type::lat,3,3,LATfield2::matrix_symmetry::symmetric),
        
        // initialize k-fields, sources
        S00_FT(base_type::latFT,1),
        S0i_FT(base_type::latFT,3),
        Sij_FT(base_type::latFT,3,3,LATfield2::matrix_symmetry::symmetric),
        
        
        // initialize plans
        plan_phi(&phi, &phi_FT),
        plan_chi(&chi, &chi_FT),
        plan_Bi (&Bi, &Bi_FT),
        
        // plan_T00(&T00,&T00_FT),
        // plan_T0i(&T0i,&T0i_FT),
        // plan_Tij(&Tij,&Tij_FT),
        
        plan_S00(&S00,&S00_FT),
        plan_S0i(&S0i,&S0i_FT),
        plan_Sij(&Sij,&Sij_FT)
    {
        scalar_to_zero(phi);
        scalar_to_zero(chi);
        vector_to_zero(Bi);
        
        scalar_to_zero(T00);
        vector_to_zero(T0i);
        tensor_to_zero(Tij);
        
        scalar_to_zero(S00);
        vector_to_zero(S0i);
        tensor_to_zero(Sij);
    }
    
    void clear_sources() override
    {
        scalar_to_zero(T00);
        vector_to_zero(T0i);
        tensor_to_zero(Tij);
    }
    
    /*
        sample particle masses into the source field
    */
    void sample(const particle_container& pcls, double a) override
    // TODO: we don't actually need the scale factor here if we use particles'
    // canonical momentum normalized q = p/mca.
    {
        // WARNING: has phi been initialized? 
        projection_T00_project(&pcls, &T00, a, &phi); // samples
        projection_T00_comm (&T00); // communicates the ghost cells
        
        projection_T0i_project(&pcls,&T0i,&phi);
        projection_T0i_comm(&T0i);
        
        projection_Tij_project(&pcls,&Tij,a,&phi);
        projection_Tij_comm(&Tij);
    }
    
    void compute_phi(
        double a, double Hc, double fourpiG, double Omega)
    {
        const double dx = 1.0/base_type::lat.size()[0];
        prepareFTsource<real_type> (
            phi, 
            chi, 
            T00,
            Omega,
            S00, 
            // 3. * Hc * dx * dx / dt,
            fourpiG * dx * dx / a,
            3. * Hc * Hc * dx * dx); 
        plan_S00.execute (LATfield2::FFT_FORWARD);
        S00_FT.updateHalo ();
        solveModifiedPoissonFT (/* source = */ S00_FT, 
                                /* poten. = */ phi_FT, 
                                1. / (dx * dx) );
                                // 3. * Hc/ dt);
        plan_phi.execute (LATfield2::FFT_BACKWARD); // go back to position space
        phi.updateHalo (); // update ghost cells
    }
    void compute_chi(double f = 1.0)
    {
        prepareFTsource<real_type> (
            phi, 
            Tij, 
            Sij,
            2. * f);
        plan_Sij.execute (LATfield2::FFT_FORWARD);
        Sij_FT.updateHalo ();
        projectFTscalar (Sij_FT,chi_FT);
        plan_chi.execute(LATfield2::FFT_BACKWARD);
        chi.updateHalo();
    }
    void compute_Bi(double f = 1.0)
    {
        // prepare source S0i from T0i
        Site x(base_type::lat);
        for(x.first();x.test();x.next())
        {
            for(int i=0;i<3;++i)
                S0i(x,i) = T0i(x,i);
        }
        
        plan_S0i.execute(LATfield2::FFT_FORWARD);
        S0i_FT.updateHalo();
        
        projectFTvector (S0i_FT, Bi_FT, f);
        
        plan_Bi.execute(LATfield2::FFT_BACKWARD);
        Bi.updateHalo();
    }
    void compute_potential(
        double fourpiG, 
        double a, 
        double Hc, 
        double Omega) override
    // TODO: can we remove all of these dependencies?
    {
        const double dx = 1.0/base_type::lat.size()[0];
        
        compute_phi(a,Hc,fourpiG,Omega);
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
        site_type x(phi.lattice());
        for (x.first(); x.test(); x.next())
        {
            phi(x) *= f({x.coord(0),x.coord(1),x.coord(2)});
        }
        phi.updateHalo();
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
    
    //std::array<real_type,3> gradient_vector(
    //    const real_field_type& F, 
    //    const std::array<real_type,3>& momentum,
    //    const LATfield2::Site& x,
    //    const std::array<real_type,3>& pos)const
    //// First order CIC gradient
    //{
    //    const int N = size();
    //    const real_type dx = 1.0 / N;
    //    
    //    std::array<real_type,3> ref_dist{0,0,0};
    //    for(int i=0;i<3;++i)
    //        ref_dist[i] = pos[i]/dx - x.coord(i);
    //        
    //    std::array<real_type,3> grad{0,0,0};
    //    for(int i=0;i<3;++i)
    //    {
    //        const int j=(i+1)%3,k=(j+1)%3;
    //        grad[i] += (1. - ref_dist[j]) * (1. - ref_dist[k])
    //                     * (F (x + i) - F (x));
    //        grad[i] += ref_dist[j] * (1. - ref_dist[k])
    //                      * (F (x + i + j) - F (x + j));
    //        grad[i] += (1. - ref_dist[j]) * ref_dist[k]
    //                      * (F (x + i + k) - F (x + k));
    //        grad[i] += ref_dist[j] * ref_dist[k]
    //                      * (F (x + i + j + k) - F (x + j + k));
    //    }
    //    return grad;
    //}
    
    /*
        compute forces
    */
    void compute_forces(particle_container& pcls, 
        double /* fourpiG */, 
        double a) const override
    {
        // const int N = size();
        // const real_type dx = 1.0 / N;
        
        site_type xpart(pcls.lattice());
        for(xpart.first();xpart.test();xpart.next())
        {
            for(auto& part : pcls.field()(xpart).parts )
            {
                std::array<real_type,3> pos{part.pos[0],part.pos[1],part.pos[2]};
                std::array<real_type,3> 
                    gradphi = gradient(phi,xpart,pos), 
                    gradchi = gradient(chi,xpart,pos), 
                    pgradB{0,0,0};
                    //pgradB  = gradient_vector( Bi,momentum,xpart,pos );
                    
                real_type p2 =   part.momentum[0] * part.momentum[0] 
                          + part.momentum[1] * part.momentum[1] 
                          + part.momentum[2] * part.momentum[2];
                real_type e2 = std::sqrt(p2 + a*a);
        
        const int N = size();
        const real_type dx = 1.0 / N;
        
        std::array<real_type,3> ref_dist{0,0,0};
        for(int i=0;i<3;++i)
            ref_dist[i] = pos[i]/dx - xpart.coord(i);


                pgradB[0] = ((1. - ref_dist[2]) * (Bi (xpart + 0, 1) - Bi (xpart, 1))
                             + ref_dist[2] * (Bi (xpart + 2 + 0, 1) - Bi (xpart + 2, 1)))
                            * part.momentum[1];
                pgradB[0] += ((1. - ref_dist[1]) * (Bi (xpart + 0, 2) - Bi (xpart, 2))
                              + ref_dist[1] * (Bi (xpart + 1 + 0, 2) - Bi (xpart + 1, 2)))
                             * part.momentum[2];
                pgradB[0] += (1. - ref_dist[1]) * (1. - ref_dist[2])
                             * ((ref_dist[0] - 1.) * Bi (xpart - 0, 0)
                                + (1. - 2. * ref_dist[0]) * Bi (xpart, 0)
                                + ref_dist[0] * Bi (xpart + 0, 0))
                             * part.momentum[0];
                pgradB[0] += ref_dist[1] * (1. - ref_dist[2])
                             * ((ref_dist[0] - 1.) * Bi (xpart + 1 - 0, 0)
                                + (1. - 2. * ref_dist[0]) * Bi (xpart + 1, 0)
                                + ref_dist[0] * Bi (xpart + 1 + 0, 0))
                             * part.momentum[0];
                pgradB[0] += (1. - ref_dist[1]) * ref_dist[2]
                             * ((ref_dist[0] - 1.) * Bi (xpart + 2 - 0, 0)
                                + (1. - 2. * ref_dist[0]) * Bi (xpart + 2, 0)
                                + ref_dist[0] * Bi (xpart + 2 + 0, 0))
                             * part.momentum[0];
                pgradB[0] += ref_dist[1] * ref_dist[2]
                             * ((ref_dist[0] - 1.) * Bi (xpart + 2 + 1 - 0, 0)
                                + (1. - 2. * ref_dist[0]) * Bi (xpart + 2 + 1, 0)
                                + ref_dist[0] * Bi (xpart + 2 + 1 + 0, 0))
                             * part.momentum[0];

                pgradB[1] = ((1. - ref_dist[0]) * (Bi (xpart + 1, 2) - Bi (xpart, 2))
                             + ref_dist[0] * (Bi (xpart + 1 + 0, 2) - Bi (xpart + 0, 2)))
                            * part.momentum[2];
                pgradB[1] += ((1. - ref_dist[2]) * (Bi (xpart + 1, 0) - Bi (xpart, 0))
                              + ref_dist[2] * (Bi (xpart + 1 + 2, 0) - Bi (xpart + 2, 0)))
                             * part.momentum[0];
                pgradB[1] += (1. - ref_dist[0]) * (1. - ref_dist[2])
                             * ((ref_dist[1] - 1.) * Bi (xpart - 1, 1)
                                + (1. - 2. * ref_dist[1]) * Bi (xpart, 1)
                                + ref_dist[1] * Bi (xpart + 1, 1))
                             * part.momentum[1];
                pgradB[1] += ref_dist[0] * (1. - ref_dist[2])
                             * ((ref_dist[1] - 1.) * Bi (xpart + 0 - 1, 1)
                                + (1. - 2. * ref_dist[1]) * Bi (xpart + 0, 1)
                                + ref_dist[1] * Bi (xpart + 0 + 1, 1))
                             * part.momentum[1];
                pgradB[1] += (1. - ref_dist[0]) * ref_dist[2]
                             * ((ref_dist[1] - 1.) * Bi (xpart + 2 - 1, 1)
                                + (1. - 2. * ref_dist[1]) * Bi (xpart + 2, 1)
                                + ref_dist[1] * Bi (xpart + 2 + 1, 1))
                             * part.momentum[1];
                pgradB[1] += ref_dist[0] * ref_dist[2]
                             * ((ref_dist[1] - 1.) * Bi (xpart + 2 + 0 - 1, 1)
                                + (1. - 2. * ref_dist[1]) * Bi (xpart + 2 + 0, 1)
                                + ref_dist[1] * Bi (xpart + 2 + 0 + 1, 1))
                             * part.momentum[1];

                pgradB[2] = ((1. - ref_dist[1]) * (Bi (xpart + 2, 0) - Bi (xpart, 0))
                             + ref_dist[1] * (Bi (xpart + 2 + 1, 0) - Bi (xpart + 1, 0)))
                            * part.momentum[0];
                pgradB[2] += ((1. - ref_dist[0]) * (Bi (xpart + 2, 1) - Bi (xpart, 1))
                              + ref_dist[0] * (Bi (xpart + 2 + 0, 1) - Bi (xpart + 0, 1)))
                             * part.momentum[1];
                pgradB[2] += (1. - ref_dist[0]) * (1. - ref_dist[1])
                             * ((ref_dist[2] - 1.) * Bi (xpart - 2, 2)
                                + (1. - 2. * ref_dist[2]) * Bi (xpart, 2)
                                + ref_dist[2] * Bi (xpart + 2, 2))
                             * part.momentum[2];
                pgradB[2] += ref_dist[0] * (1. - ref_dist[1])
                             * ((ref_dist[2] - 1.) * Bi (xpart + 0 - 2, 2)
                                + (1. - 2. * ref_dist[2]) * Bi (xpart + 0, 2)
                                + ref_dist[2] * Bi (xpart + 0 + 2, 2))
                             * part.momentum[2];
                pgradB[2] += (1. - ref_dist[0]) * ref_dist[1]
                             * ((ref_dist[2] - 1.) * Bi (xpart + 1 - 2, 2)
                                + (1. - 2. * ref_dist[2]) * Bi (xpart + 1, 2)
                                + ref_dist[2] * Bi (xpart + 2 + 1, 2))
                             * part.momentum[2];
                pgradB[2] += ref_dist[0] * ref_dist[1]
                             * ((ref_dist[2] - 1.) * Bi (xpart + 1 + 0 - 2, 2)
                                + (1. - 2. * ref_dist[2]) * Bi (xpart + 1 + 0, 2)
                                + ref_dist[2] * Bi (xpart + 1 + 0 + 2, 2))
                             * part.momentum[2];
                
                for(int i=0;i<3;++i)
                {
                    part.force[i] = (gradchi[i]*e2 - gradphi[i] * (e2 + p2/e2) 
                                    - pgradB[i]/a/a/N)/dx;
                }
            }
        }
    }
   
    std::array<real_type,3> vector_at(
        const real_field_type& F, 
        const site_type& x,
        const std::array<real_type,3>& pos)const
    {
        const int N = size();
        const real_type dx = 1.0 / N;
        
        std::array<real_type,3> ref_dist{0,0,0};
        for(int i=0;i<3;++i)
            ref_dist[i] = pos[i]/dx - x.coord(i);
            
        std::array<real_type,3> xF{0,0,0};
        
        // CIC interpolation
        for(int i=0;i<3;++i)
        {
            xF[i] += F(x,i) * (1.-ref_dist[0]) * (1.-ref_dist[1]) * (1.-ref_dist[2]);
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
    real_type scalar_at(
        const real_field_type& F, 
        const site_type& x,
        const std::array<real_type,3>& pos)const
    {
        const int N = size();
        const real_type dx = 1.0 / N;
        
        std::array<real_type,3> ref_dist{0,0,0};
        for(int i=0;i<3;++i)
            ref_dist[i] = pos[i]/dx - x.coord(i);
            
        real_type xF{0};
        
        // CIC interpolation
        xF += F(x) * (1.-ref_dist[0]) * (1.-ref_dist[1]) * (1.-ref_dist[2]);
        xF += F(x+0) * ref_dist[0] * (1.-ref_dist[1]) * (1.-ref_dist[2]);
        xF += F(x+1) * (1.-ref_dist[0]) * ref_dist[1] * (1.-ref_dist[2]);
        xF += F(x+0+1) * ref_dist[0] * ref_dist[1] * (1.-ref_dist[2]);
        xF += F(x+2) * (1.-ref_dist[0]) * (1.-ref_dist[1]) * ref_dist[2];
        xF += F(x+0+2) * ref_dist[0] * (1.-ref_dist[1]) * ref_dist[2];
        xF += F(x+1+2) * (1.-ref_dist[0]) * ref_dist[1] * ref_dist[2];
        xF += F(x+0+1+2) * ref_dist[0] * ref_dist[1] * ref_dist[2];
        
        return xF;
    }
   
    std::array<real_type,3> momentum_to_velocity(
                          const std::array<real_type,3>& momentum,
                          const std::array<real_type,3>& position,
                          const site_type& xpart,
                          const real_type a) const override
    {
        const int N = size();
        std::array<real_type,3> velocity{0,0,0};
        real_type xchi{0},xphi{0};
        std::array<real_type,3> xBi{0,0,0};
        const real_type momentum2 =    momentum[0]*momentum[0] 
                                + momentum[1]*momentum[1]
                                + momentum[2]*momentum[2];
        const real_type e2 = std::sqrt(momentum2 + a*a);
                
        
        xphi = scalar_at(phi,xpart,{position[0],position[1],position[2]});
        xchi = scalar_at(chi,xpart,{position[0],position[1],position[2]});
        xBi  = vector_at(Bi,xpart,{position[0],position[1],position[2]});
        
        const real_type velocity2 = 
            (1. + (3. - momentum2 / e2/e2) * xphi - xchi) / e2;
        
        for(int i=0;i<3;++i)
        {
            velocity[i] = momentum[i] * velocity2 + xBi[i]/a/a/N; // TODO: why this strange scaling of Bi?
        }
        
        return velocity;
    }
    // TODO: check if the inversion works
    std::array<real_type,3> velocity_to_momentum(
                          const std::array<real_type,3>& velocity,
                          const std::array<real_type,3>& position,
                          const site_type& xpart,
                          const real_type a) const override
    {
        const int N = size();
        std::array<real_type,3> momentum{0,0,0};
        real_type xchi{0},xphi{0};
        std::array<real_type,3> xBi{0,0,0};
        
        xphi = scalar_at(phi,xpart,{position[0],position[1],position[2]});
        xchi = scalar_at(chi,xpart,{position[0],position[1],position[2]});
        xBi  = vector_at(Bi,xpart,{position[0],position[1],position[2]});
        
        for(int i=0;i<3;++i) xBi[i] = xBi[i]/a/a/N; // TODO: why the strange scaling of Bi?
        
        real_type u2{0};
        for(int i=0;i<3;++i)
            u2 += velocity[i]*velocity[i];
        
        real_type Biui{0};
        
        for(int i=0;i<3;++i)
            Biui += velocity[i]*xBi[i];
            
        const real_type fact_B = 1.0/std::sqrt(1-u2);
        const real_type fact_A = 1 - fact_B*fact_B*(Biui - xchi + xphi*(3-u2));
        
        for(int i=0;i<3;++i)
        {
            momentum[i] = (velocity[i]*fact_A - xBi[i])*fact_B;
        }
        
        return momentum;
    }
    
    virtual ~relativistic_pm(){}
    std::string report() const override
    {
        std::stringstream ss;
        // sources
        double std_t00 = show_msq(T00);
        double std_t0i = show_msq(T0i,0);
        double std_tij = show_msq(Tij,0,0);
        // potentials
        double std_phi = show_msq(phi);
        double std_chi = show_msq(chi);
        double std_Bi = show_msq(Bi,0);
        ss << "RMS(T00) = " << std_t00 << '\n';
        ss << "RMS(T0i) = " << std_t0i << '\n';
        ss << "RMS(Tij) = " << std_tij << '\n';
        ss << "RMS(Phi) = " << std_phi << '\n';
        ss << "RMS(Chi) = " << std_chi << '\n';
        ss << "RMS(Bi)  = " << std_Bi << '\n';
        return ss.str();
    }
    double density() const override
    {
        return show_mean(T00);
    }
    double sum_phi() const override
    {
        // const double N = size();
        // return show_mean(phi) * N * N *N;
        return show_mean(phi);
    }
    virtual void save_to_file(std::string prefix) const 
    {
        // save the energy-momentum tensor
        T00.saveHDF5 (prefix + "_T00.h5");
        T0i.saveHDF5 (prefix + "_T0i.h5");
        Tij.saveHDF5 (prefix + "_Tij.h5");
        
        // save the potentials
        phi.saveHDF5 (prefix + "_phi.h5");
        chi.saveHDF5 (prefix + "_chi.h5");
        Bi.saveHDF5 (prefix + "_B.h5");
    }
};
} // namespace gevolution
