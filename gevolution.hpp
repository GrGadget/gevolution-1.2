//////////////////////////
// gevolution.hpp
//////////////////////////
//
// Geneva algorithms for evolution of metric perturbations
// and relativistic free-streaming particles (gevolution)
//
// 1. Suite of Fourier-based methods for the computation of the
//    relativistic scalar (Phi, Phi-Psi) and vector modes [see J. Adamek,
//    R. Durrer, and M. Kunz, Class. Quant. Grav. 31, 234006 (2014)]
//
// 2. Collection of "update position" and "update velocity/momentum" methods
//    [see J. Adamek, D. Daverio, R. Durrer, and M. Kunz, JCAP 1607, 053
//    (2016)]
//
// 3. Collection of projection methods for the construction of the
//    stress-energy-tensor
//
// 4. Fourier-space projection methods for the computation of the
//    curl and divergence of the velocity field
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen
// Mary University of London)
//
// Last modified: April 2019
//
//////////////////////////

#ifndef GEVOLUTION_HEADER
#define GEVOLUTION_HEADER

#include "LATfield2.hpp"
#include "real_type.hpp"
#include "Particles_gevolution.hpp"
#include <cstdlib>
#include <iostream>

namespace gevolution
{
using LATfield2::Field;
using LATfield2::part_simple_info;
using LATfield2::Particles;
using LATfield2::rKSite;
using LATfield2::Site;

//////////////////////////
// prepareFTsource (1)
//////////////////////////
// Description:
//   construction of real-space source tensor for Fourier-based solvers
//
// Arguments:
//   phi        reference to field configuration
//   Tij        reference to symmetric tensor field containing the space-space
//              components of the stress-energy tensor (rescaled by a^3)
//   Sij        reference to allocated symmetric tensor field which will
//   contain
//              the source tensor (may be identical to Tji)
//   coeff      scaling coefficient for Tij ("8 pi G dx^2 / a")
//
// Returns:
//
//////////////////////////

template <class FieldType>
void prepareFTsource (Field<FieldType> &phi, Field<FieldType> &Tij,
                      Field<FieldType> &Sij, const double coeff)
{
    Site x (phi.lattice ());

    for (x.first (); x.test (); x.next ())
    {
        // 0-0-component:
        Sij (x, 0, 0) = coeff * Tij (x, 0, 0);
#ifdef PHINONLINEAR
#ifdef ORIGINALMETRIC
        Sij (x, 0, 0)
            -= 4. * phi (x) * (phi (x - 0) + phi (x + 0) - 2. * phi (x));
        Sij (x, 0, 0)
            -= 0.5 * (phi (x + 0) - phi (x - 0)) * (phi (x + 0) - phi (x - 0));
#else
        Sij (x, 0, 0)
            += 0.5 * (phi (x + 0) - phi (x - 0)) * (phi (x + 0) - phi (x - 0));
#endif
#endif

        // 1-1-component:
        Sij (x, 1, 1) = coeff * Tij (x, 1, 1);
#ifdef PHINONLINEAR
#ifdef ORIGINALMETRIC
        Sij (x, 1, 1)
            -= 4. * phi (x) * (phi (x - 1) + phi (x + 1) - 2. * phi (x));
        Sij (x, 1, 1)
            -= 0.5 * (phi (x + 1) - phi (x - 1)) * (phi (x + 1) - phi (x - 1));
#else
        Sij (x, 1, 1)
            += 0.5 * (phi (x + 1) - phi (x - 1)) * (phi (x + 1) - phi (x - 1));
#endif
#endif

        // 2-2-component:
        Sij (x, 2, 2) = coeff * Tij (x, 2, 2);
#ifdef PHINONLINEAR
#ifdef ORIGINALMETRIC
        Sij (x, 2, 2)
            -= 4. * phi (x) * (phi (x - 2) + phi (x + 2) - 2. * phi (x));
        Sij (x, 2, 2)
            -= 0.5 * (phi (x + 2) - phi (x - 2)) * (phi (x + 2) - phi (x - 2));
#else
        Sij (x, 2, 2)
            += 0.5 * (phi (x + 2) - phi (x - 2)) * (phi (x + 2) - phi (x - 2));
#endif
#endif

        // 0-1-component:
        Sij (x, 0, 1) = coeff * Tij (x, 0, 1);
#ifdef PHINONLINEAR
        Sij (x, 0, 1) += phi (x + 0) * phi (x + 1) - phi (x) * phi (x + 0 + 1);
#ifdef ORIGINALMETRIC
        Sij (x, 0, 1) -= 1.5 * phi (x) * phi (x);
        Sij (x, 0, 1) += 1.5 * phi (x + 0) * phi (x + 0);
        Sij (x, 0, 1) += 1.5 * phi (x + 1) * phi (x + 1);
        Sij (x, 0, 1) -= 1.5 * phi (x + 0 + 1) * phi (x + 0 + 1);
#else
        Sij (x, 0, 1) += 0.5 * phi (x) * phi (x);
        Sij (x, 0, 1) -= 0.5 * phi (x + 0) * phi (x + 0);
        Sij (x, 0, 1) -= 0.5 * phi (x + 1) * phi (x + 1);
        Sij (x, 0, 1) += 0.5 * phi (x + 0 + 1) * phi (x + 0 + 1);
#endif
#endif

        // 0-2-component:
        Sij (x, 0, 2) = coeff * Tij (x, 0, 2);
#ifdef PHINONLINEAR
        Sij (x, 0, 2) += phi (x + 0) * phi (x + 2) - phi (x) * phi (x + 0 + 2);
#ifdef ORIGINALMETRIC
        Sij (x, 0, 2) -= 1.5 * phi (x) * phi (x);
        Sij (x, 0, 2) += 1.5 * phi (x + 0) * phi (x + 0);
        Sij (x, 0, 2) += 1.5 * phi (x + 2) * phi (x + 2);
        Sij (x, 0, 2) -= 1.5 * phi (x + 0 + 2) * phi (x + 0 + 2);
#else
        Sij (x, 0, 2) += 0.5 * phi (x) * phi (x);
        Sij (x, 0, 2) -= 0.5 * phi (x + 0) * phi (x + 0);
        Sij (x, 0, 2) -= 0.5 * phi (x + 2) * phi (x + 2);
        Sij (x, 0, 2) += 0.5 * phi (x + 0 + 2) * phi (x + 0 + 2);
#endif
#endif

        // 1-2-component:
        Sij (x, 1, 2) = coeff * Tij (x, 1, 2);
#ifdef PHINONLINEAR
        Sij (x, 1, 2) += phi (x + 1) * phi (x + 2) - phi (x) * phi (x + 1 + 2);
#ifdef ORIGINALMETRIC
        Sij (x, 1, 2) -= 1.5 * phi (x) * phi (x);
        Sij (x, 1, 2) += 1.5 * phi (x + 1) * phi (x + 1);
        Sij (x, 1, 2) += 1.5 * phi (x + 2) * phi (x + 2);
        Sij (x, 1, 2) -= 1.5 * phi (x + 1 + 2) * phi (x + 1 + 2);
#else
        Sij (x, 1, 2) += 0.5 * phi (x) * phi (x);
        Sij (x, 1, 2) -= 0.5 * phi (x + 1) * phi (x + 1);
        Sij (x, 1, 2) -= 0.5 * phi (x + 2) * phi (x + 2);
        Sij (x, 1, 2) += 0.5 * phi (x + 1 + 2) * phi (x + 1 + 2);
#endif
#endif
    }
}

//////////////////////////
// prepareFTsource (2)
//////////////////////////
// Description:
//   construction of real-space source field for Fourier-based solvers
//
// Arguments:
//   phi        reference to field configuration (first Bardeen potential)
//   chi        reference to field configuration (difference between Bardeen
//   potentials, phi-psi) source     reference to fully dressed source field
//   (rescaled by a^3) bgmodel    background model of the source (rescaled by
//   a^3) to be subtracted result     reference to allocated field which will
//   contain the result (may be identical to source) coeff      diffusion
//   coefficient ("3 H_conformal dx^2 / dtau") coeff2     scaling coefficient
//   for the source ("4 pi G dx^2 / a") coeff3     scaling coefficient for the
//   psi-term ("3 H_conformal^2 dx^2")
//
// Returns:
//
//////////////////////////

template <class FieldType>
void prepareFTsource (Field<FieldType> &phi, Field<FieldType> &chi,
                      Field<FieldType> &source, const FieldType bgmodel,
                      Field<FieldType> &result, const double coeff,
                      const double coeff2, const double coeff3)
{
    Site x (phi.lattice ());

    for (x.first (); x.test (); x.next ())
    {
        result (x) = coeff2 * (source (x) - bgmodel);
#ifdef PHINONLINEAR
#ifdef ORIGINALMETRIC
        result (x) *= 1. - 4. * phi (x);
        result (x) -= 0.375 * (phi (x - 0) - phi (x + 0))
                      * (phi (x - 0) - phi (x + 0));
        result (x) -= 0.375 * (phi (x - 1) - phi (x + 1))
                      * (phi (x - 1) - phi (x + 1));
        result (x) -= 0.375 * (phi (x - 2) - phi (x + 2))
                      * (phi (x - 2) - phi (x + 2));
#else
        result (x) *= 1. - 2. * phi (x);
        result (x) += 0.125 * (phi (x - 0) - phi (x + 0))
                      * (phi (x - 0) - phi (x + 0));
        result (x) += 0.125 * (phi (x - 1) - phi (x + 1))
                      * (phi (x - 1) - phi (x + 1));
        result (x) += 0.125 * (phi (x - 2) - phi (x + 2))
                      * (phi (x - 2) - phi (x + 2));
#endif
#endif
        result (x) += (coeff3 - coeff) * phi (x) - coeff3 * chi (x);
    }
}

#ifdef FFT3D
//////////////////////////
// projectFTscalar
//////////////////////////
// Description:
//   projection of the Fourier image of a tensor field on the trace-free
//   longitudinal (scalar) component
//
// Arguments:
//   SijFT      reference to the Fourier image of the input tensor field
//   chiFT      reference to allocated field which will contain the Fourier
//              image of the trace-free longitudinal (scalar) component
//
// Returns:
//
//////////////////////////

void projectFTscalar (Field<Cplx> &SijFT, Field<Cplx> &chiFT,
                      const int add = 0);

//////////////////////////
// evolveFTvector
//////////////////////////
// Description:
//   projects the Fourier image of a tensor field on the spin-1 component
//   used as a source for the evolution of the vector perturbation
//
// Arguments:
//   SijFT      reference to the Fourier image of the input tensor field
//   BiFT       reference to the Fourier image of the vector perturbation
//   a2dtau     conformal time step times scale factor squared (a^2 * dtau)
//
// Returns:
//
//////////////////////////

void evolveFTvector (Field<Cplx> &SijFT, Field<Cplx> &BiFT, const Real a2dtau);

//////////////////////////
// projectFTvector
//////////////////////////
// Description:
//   projects the Fourier image of a vector field on the transverse component
//   and solves the constraint equation for the vector perturbation
//
// Arguments:
//   SiFT       reference to the Fourier image of the input vector field
//   BiFT       reference to the Fourier image of the vector perturbation (can
//   be identical to input) coeff      rescaling coefficient (default 1) modif
//   modification k^2 -> k^2 + modif (default 0)
//
// Returns:
//
//////////////////////////

void projectFTvector (Field<Cplx> &SiFT, Field<Cplx> &BiFT,
                      const Real coeff = 1., const Real modif = 0.);

//////////////////////////
// projectFTtensor
//////////////////////////
// Description:
//   projection of the Fourier image of a tensor field on the transverse
//   trace-free tensor component
//
// Arguments:
//   SijFT      reference to the Fourier image of the input tensor field
//   hijFT      reference to allocated field which will contain the Fourier
//              image of the transverse trace-free tensor component
//
// Returns:
//
//////////////////////////

void projectFTtensor (Field<Cplx> &SijFT, Field<Cplx> &hijFT);

//////////////////////////
// solveModifiedPoissonFT
//////////////////////////
// Description:
//   Modified Poisson solver using the standard Fourier method
//
// Arguments:
//   sourceFT   reference to the Fourier image of the source field
//   potFT      reference to the Fourier image of the potential
//   coeff      coefficient applied to the source ("4 pi G / a")
//   modif      modification k^2 -> k^2 + modif (default 0 gives standard
//   Poisson equation)
//
// Returns:
//
//////////////////////////

void solveModifiedPoissonFT (Field<Cplx> &sourceFT, Field<Cplx> &potFT,
                             Real coeff, const Real modif = 0.);
#endif

//////////////////////////
// update_q
//////////////////////////
// Description:
//   Update momentum method (arbitrary momentum)
//   Note that vel[3] in the particle structure is used to store q[3] in units
//   of the particle mass, such that as q^2 << m^2 a^2 the meaning of vel[3]
//   is ~ v*a.
//
// Arguments:
//   dtau       time step
//   dx         lattice unit
//   part       pointer to particle structure
//   ref_dist   distance vector to reference point
//   partInfo   global particle properties (unused)
//   fields     array of pointers to fields appearing in geodesic equation
//              fields[0] = phi
//              fields[1] = chi
//              fields[2] = Bi
//   sites      array of sites on the respective lattices
//   nfield     number of fields
//   params     array of additional parameters
//              params[0] = a
//              params[1] = scaling coefficient for Bi
//   outputs    array of reduction variables
//   noutputs   number of reduction variables
//
// Returns: squared velocity of particle after update
//
//////////////////////////

Real update_q (double dtau, double dx, particle *part, double *ref_dist,
               part_simple_info partInfo, Field<Real> **fields, Site *sites,
               int nfield, double *params, double *outputs, int noutputs);

//////////////////////////
// update_q_Newton
//////////////////////////
// Description:
//   Update momentum method (Newtonian version)
//   Note that vel[3] in the particle structure is used to store q[3] in units
//   of the particle mass, such that the meaning of vel[3] is v*a.
//
// Arguments:
//   dtau       time step
//   dx         lattice unit
//   part       pointer to particle structure
//   ref_dist   distance vector to reference point
//   partInfo   global particle properties (unused)
//   fields     array of pointers to fields appearing in geodesic equation
//              fields[0] = psi
//              fields[1] = chi
//   sites      array of sites on the respective lattices
//   nfield     number of fields (should be 1)
//   params     array of additional parameters
//              params[0] = a
//   outputs    array of reduction variables
//   noutputs   number of reduction variables
//
// Returns: squared velocity of particle after update
//
//////////////////////////

//Real update_q_Newton (double dtau, double dx, part_simple *part,
//                      double *ref_dist, part_simple_info partInfo,
//                      Field<Real> **fields, Site *sites, int nfield,
//                      double *params, double *outputs, int noutputs);
// Real update_q_Newton (double dtau, double dx, part_simple *part,
//                       double *ref_dist, part_simple_info partInfo,
//                       Field<Real> **fields, Site *sites, int nfield,
//                       double *params, double *outputs, int noutputs);
Real update_q_Newton ( 
                      particle& part,
                      const Field<Real>& psi, 
                      const Site& xpart,
                      double dtau,
                      double dx,
                      double a);
Real update_q_Newton ( 
                      particle& part,
                      const double dtau);

//////////////////////////
// update_pos
//////////////////////////
// Description:
//   Update position method (arbitrary momentum)
//   Note that vel[3] in the particle structure is used to store q[3] in units
//   of the particle mass, such that as q^2 << m^2 a^2 the meaning of vel[3]
//   is ~ v*a.
//
// Arguments:
//   dtau       time step
//   dx         lattice unit
//   part       pointer to particle structure
//   ref_dist   distance vector to reference point
//   partInfo   global particle properties (unused)
//   fields     array of pointers to fields appearing in geodesic equation
//              fields[0] = phi
//              fields[1] = chi
//              fields[2] = Bi
//   sites      array of sites on the respective lattices
//   nfield     number of fields
//   params     array of additional parameters
//              params[0] = a
//              params[1] = scaling coefficient for Bi
//   outputs    array of reduction variables
//   noutputs   number of reduction variables
//
// Returns:
//
//////////////////////////

void update_pos (double dtau, double dx, particle *part, double *ref_dist,
                 part_simple_info partInfo, Field<Real> **fields, Site *sites,
                 int nfield, double *params, double *outputs, int noutputs);

//////////////////////////
// update_pos_Newton
//////////////////////////
// Description:
//   Update position method (Newtonian version)
//   Note that vel[3] in the particle structure is used to store q[3] in units
//   of the particle mass, such that the meaning of vel[3] is v*a.
//
// Arguments:
//   dtau       time step
//   dx         lattice unit (unused)
//   part       pointer to particle structure
//   ref_dist   distance vector to reference point (unused)
//   partInfo   global particle properties (unused)
//   fields     array of pointers to fields appearing in geodesic equation
//   (unused) sites      array of sites on the respective lattices (unused)
//   nfield     number of fields (unused)
//   params     array of additional parameters
//              params[0] = a
//   outputs    array of reduction variables (unused)
//   noutputs   number of reduction variables (unused)
//
// Returns:
//
//////////////////////////

void update_pos_Newton (double dtau, double dx, particle *part,
                        double *ref_dist, part_simple_info partInfo,
                        Field<Real> **fields, Site *sites, int nfield,
                        double *params, double *outputs, int noutputs);

//////////////////////////
// projection_T00_project
//////////////////////////
// Description:
//   Particle-mesh projection for T00, including geometric corrections
//
// Arguments:
//   pcls       pointer to particle handler
//   T00        pointer to target field
//   a          scale factor at projection (needed in order to convert
//              canonical momenta to energies)
//   phi        pointer to Bardeen potential which characterizes the
//              geometric corrections (volume distortion); can be set to
//              NULL which will result in no corrections applied
//   coeff      coefficient applied to the projection operation (default 1)
//
// Returns:
//
//////////////////////////

template <typename part, typename part_info, typename part_dataType>
void projection_T00_project (Particles<part, part_info, part_dataType> *pcls,
                             Field<Real> *T00, double a = 1.,
                             Field<Real> *phi = NULL, double coeff = 1.)
{
    if (T00->lattice ().halo () == 0)
    {
        std::cout << "projection_T00_project: target field needs halo > 0"
                  << std::endl;
        std::exit (-1);
    }

    Site xPart (pcls->lattice ());
    Site xField (T00->lattice ());

    Real referPos[3];
    Real weightScalarGridUp[3];
    Real weightScalarGridDown[3];
    Real dx = pcls->res ();

    double mass = coeff / (dx * dx * dx);
    mass *= pcls->parts_info () -> mass;
    mass /= a;

    Real e = a, f = 0.;

    Real localCube[8]; // XYZ = 000 | 001 | 010 | 011 | 100 | 101 | 110 | 111
    Real localCubePhi[8];

    for (int i = 0; i < 8; i++)
        localCubePhi[i] = 0.0;

    for (xPart.first (), xField.first (); xPart.test ();
         xPart.next (), xField.next ())
    {
        if (pcls->field () (xPart).size != 0)
        {
            for (int i = 0; i < 3; i++)
                referPos[i] = xPart.coord (i) * dx;
            for (int i = 0; i < 8; i++)
                localCube[i] = 0.0;

            if (phi != NULL)
            {
                localCubePhi[0] = (*phi) (xField);
                localCubePhi[1] = (*phi) (xField + 2);
                localCubePhi[2] = (*phi) (xField + 1);
                localCubePhi[3] = (*phi) (xField + 1 + 2);
                localCubePhi[4] = (*phi) (xField + 0);
                localCubePhi[5] = (*phi) (xField + 0 + 2);
                localCubePhi[6] = (*phi) (xField + 0 + 1);
                localCubePhi[7] = (*phi) (xField + 0 + 1 + 2);
            }

            for (const auto& p : pcls->field ()(xPart).parts)
            {
                for (int i = 0; i < 3; i++)
                {
                    weightScalarGridUp[i] = (p.pos[i] - referPos[i]) / dx;
                    weightScalarGridDown[i] = 1.0l - weightScalarGridUp[i];
                }

                if (phi != NULL)
                {
                    const auto &q = p.vel;
                    f = q[0] * q[0] + q[1] * q[1] + q[2] * q[2];
                    e = sqrt (f + a * a);
                    f = 3. * e + f / e;
                }

                // 000
                localCube[0]
                    += weightScalarGridDown[0] * weightScalarGridDown[1]
                       * weightScalarGridDown[2] * (e + f * localCubePhi[0]);
                // 001
                localCube[1]
                    += weightScalarGridDown[0] * weightScalarGridDown[1]
                       * weightScalarGridUp[2] * (e + f * localCubePhi[1]);
                // 010
                localCube[2] += weightScalarGridDown[0] * weightScalarGridUp[1]
                                * weightScalarGridDown[2]
                                * (e + f * localCubePhi[2]);
                // 011
                localCube[3] += weightScalarGridDown[0] * weightScalarGridUp[1]
                                * weightScalarGridUp[2]
                                * (e + f * localCubePhi[3]);
                // 100
                localCube[4] += weightScalarGridUp[0] * weightScalarGridDown[1]
                                * weightScalarGridDown[2]
                                * (e + f * localCubePhi[4]);
                // 101
                localCube[5] += weightScalarGridUp[0] * weightScalarGridDown[1]
                                * weightScalarGridUp[2]
                                * (e + f * localCubePhi[5]);
                // 110
                localCube[6] += weightScalarGridUp[0] * weightScalarGridUp[1]
                                * weightScalarGridDown[2]
                                * (e + f * localCubePhi[6]);
                // 111
                localCube[7] += weightScalarGridUp[0] * weightScalarGridUp[1]
                                * weightScalarGridUp[2]
                                * (e + f * localCubePhi[7]);
            }

            (*T00) (xField) += localCube[0] * mass;
            (*T00) (xField + 2) += localCube[1] * mass;
            (*T00) (xField + 1) += localCube[2] * mass;
            (*T00) (xField + 1 + 2) += localCube[3] * mass;
            (*T00) (xField + 0) += localCube[4] * mass;
            (*T00) (xField + 0 + 2) += localCube[5] * mass;
            (*T00) (xField + 0 + 1) += localCube[6] * mass;
            (*T00) (xField + 0 + 1 + 2) += localCube[7] * mass;
        }
    }
}

#define projection_T00_comm scalarProjectionCIC_comm

//////////////////////////
// projection_T0i_project
//////////////////////////
// Description:
//   Particle-mesh projection for T0i, including geometric corrections
//
// Arguments:
//   pcls       pointer to particle handler
//   T0i        pointer to target field
//   phi        pointer to Bardeen potential which characterizes the
//              geometric corrections (volume distortion); can be set to
//              NULL which will result in no corrections applied
//   coeff      coefficient applied to the projection operation (default 1)
//
// Returns:
//
//////////////////////////

template <typename part, typename part_info, typename part_dataType>
void projection_T0i_project (Particles<part, part_info, part_dataType> *pcls,
                             Field<Real> *T0i, Field<Real> *phi = NULL,
                             double coeff = 1.)
{
    if (T0i->lattice ().halo () == 0)
    {
        std::cout << "projection_T0i_project: target field needs halo > 0"
                  << std::endl;
        std::exit (-1);
    }

    Site xPart (pcls->lattice ());
    Site xT0i (T0i->lattice ());

    Real referPos[3];
    Real weightScalarGridDown[3];
    Real weightScalarGridUp[3];
    Real dx = pcls->res ();

    double mass = coeff / (dx * dx * dx);
    mass *= pcls->parts_info () -> mass;

    Real w;

    Real qi[12];
    Real localCubePhi[8];

    for (int i = 0; i < 8; i++)
        localCubePhi[i] = 0;

    for (xPart.first (), xT0i.first (); xPart.test ();
         xPart.next (), xT0i.next ())
    {
        if (pcls->field () (xPart).size != 0)
        {
            for (int i = 0; i < 3; i++)
                referPos[i] = xPart.coord (i) * dx;

            for (int i = 0; i < 12; i++)
                qi[i] = 0.0;

            if (phi != NULL)
            {
                localCubePhi[0] = (*phi) (xT0i);
                localCubePhi[1] = (*phi) (xT0i + 2);
                localCubePhi[2] = (*phi) (xT0i + 1);
                localCubePhi[3] = (*phi) (xT0i + 1 + 2);
                localCubePhi[4] = (*phi) (xT0i + 0);
                localCubePhi[5] = (*phi) (xT0i + 0 + 2);
                localCubePhi[6] = (*phi) (xT0i + 0 + 1);
                localCubePhi[7] = (*phi) (xT0i + 0 + 1 + 2);
            }

            for (const auto& p: pcls->field ()(xPart).parts)
            {
                for (int i = 0; i < 3; i++)
                {
                    weightScalarGridUp[i] = (p.pos[i] - referPos[i]) / dx;
                    weightScalarGridDown[i] = 1.0l - weightScalarGridUp[i];
                }

                const auto &q = p.vel;

                w = mass * q[0];

                qi[0] += w * weightScalarGridDown[1] * weightScalarGridDown[2];
                qi[1] += w * weightScalarGridUp[1] * weightScalarGridDown[2];
                qi[2] += w * weightScalarGridDown[1] * weightScalarGridUp[2];
                qi[3] += w * weightScalarGridUp[1] * weightScalarGridUp[2];

                w = mass * q[1];

                qi[4] += w * weightScalarGridDown[0] * weightScalarGridDown[2];
                qi[5] += w * weightScalarGridUp[0] * weightScalarGridDown[2];
                qi[6] += w * weightScalarGridDown[0] * weightScalarGridUp[2];
                qi[7] += w * weightScalarGridUp[0] * weightScalarGridUp[2];

                w = mass * q[2];

                qi[8] += w * weightScalarGridDown[0] * weightScalarGridDown[1];
                qi[9] += w * weightScalarGridUp[0] * weightScalarGridDown[1];
                qi[10] += w * weightScalarGridDown[0] * weightScalarGridUp[1];
                qi[11] += w * weightScalarGridUp[0] * weightScalarGridUp[1];
            }

            (*T0i) (xT0i, 0)
                += qi[0] * (1. + localCubePhi[0] + localCubePhi[4]);
            (*T0i) (xT0i, 1)
                += qi[4] * (1. + localCubePhi[0] + localCubePhi[2]);
            (*T0i) (xT0i, 2)
                += qi[8] * (1. + localCubePhi[0] + localCubePhi[1]);

            (*T0i) (xT0i + 0, 1)
                += qi[5] * (1. + localCubePhi[4] + localCubePhi[6]);
            (*T0i) (xT0i + 0, 2)
                += qi[9] * (1. + localCubePhi[4] + localCubePhi[5]);

            (*T0i) (xT0i + 1, 0)
                += qi[1] * (1. + localCubePhi[2] + localCubePhi[6]);
            (*T0i) (xT0i + 1, 2)
                += qi[10] * (1. + localCubePhi[2] + localCubePhi[3]);

            (*T0i) (xT0i + 2, 0)
                += qi[2] * (1. + localCubePhi[1] + localCubePhi[5]);
            (*T0i) (xT0i + 2, 1)
                += qi[6] * (1. + localCubePhi[1] + localCubePhi[3]);

            (*T0i) (xT0i + 1 + 2, 0)
                += qi[3] * (1. + localCubePhi[3] + localCubePhi[7]);
            (*T0i) (xT0i + 0 + 2, 1)
                += qi[7] * (1. + localCubePhi[5] + localCubePhi[7]);
            (*T0i) (xT0i + 0 + 1, 2)
                += qi[11] * (1. + localCubePhi[6] + localCubePhi[7]);
        }
    }
}

#define projection_T0i_comm vectorProjectionCICNGP_comm

//////////////////////////
// projection_Tij_project
//////////////////////////
// Description:
//   Particle-mesh projection for Tij, including geometric corrections
//
// Arguments:
//   pcls       pointer to particle handler
//   Tij        pointer to target field
//   a          scale factor at projection (needed in order to convert
//              canonical momenta to energies)
//   phi        pointer to Bardeen potential which characterizes the
//              geometric corrections (volume distortion); can be set to
//              NULL which will result in no corrections applied
//   coeff      coefficient applied to the projection operation (default 1)
//
// Returns:
//
//////////////////////////

template <typename part, typename part_info, typename part_dataType>
void projection_Tij_project (Particles<part, part_info, part_dataType> *pcls,
                             Field<Real> *Tij, double a = 1.,
                             Field<Real> *phi = NULL, double coeff = 1.)
{
    if (Tij->lattice ().halo () == 0)
    {
        std::cout << "projection_Tij_project: target field needs halo > 0"
                  << std::endl;
        std::exit (-1);
    }

    Site xPart (pcls->lattice ());
    Site xTij (Tij->lattice ());

    Real referPos[3];
    Real weightScalarGridDown[3];
    Real weightScalarGridUp[3];
    Real dx = pcls->res ();

    double mass = coeff / (dx * dx * dx);
    mass *= pcls->parts_info () -> mass;
    mass /= a;

    Real e, f, w;

    Real tij[6];  // local cube
    Real tii[24]; // local cube
    Real localCubePhi[8];

    for (int i = 0; i < 8; i++)
        localCubePhi[i] = 0;

    for (xPart.first (), xTij.first (); xPart.test ();
         xPart.next (), xTij.next ())
    {
        if (pcls->field () (xPart).size != 0)
        {
            for (int i = 0; i < 3; i++)
                referPos[i] = (double)xPart.coord (i) * dx;

            for (int i = 0; i < 6; i++)
                tij[i] = 0.0;
            for (int i = 0; i < 24; i++)
                tii[i] = 0.0;

            if (phi != NULL)
            {
                localCubePhi[0] = (*phi) (xTij);
                localCubePhi[1] = (*phi) (xTij + 2);
                localCubePhi[2] = (*phi) (xTij + 1);
                localCubePhi[3] = (*phi) (xTij + 1 + 2);
                localCubePhi[4] = (*phi) (xTij + 0);
                localCubePhi[5] = (*phi) (xTij + 0 + 2);
                localCubePhi[6] = (*phi) (xTij + 0 + 1);
                localCubePhi[7] = (*phi) (xTij + 0 + 1 + 2);
            }

            for (const auto& p : pcls->field () (xPart).parts)
            {
                for (int i = 0; i < 3; i++)
                {
                    weightScalarGridUp[i] = (p.pos[i] - referPos[i]) / dx;
                    weightScalarGridDown[i] = 1.0l - weightScalarGridUp[i];
                }

                const auto &q = p.vel;
                f = q[0] * q[0] + q[1] * q[1] + q[2] * q[2];
                e = sqrt (f + a * a);
                f = 4. + a * a / (f + a * a);

                // diagonal components
                for (int i = 0; i < 3; i++)
                {
                    w = mass * q[i] * q[i] / e;
                    // 000
                    tii[0 + i * 8] += w * weightScalarGridDown[0]
                                      * weightScalarGridDown[1]
                                      * weightScalarGridDown[2]
                                      * (1. + f * localCubePhi[0]);
                    // 001
                    tii[1 + i * 8]
                        += w * weightScalarGridDown[0] * weightScalarGridDown[1]
                           * weightScalarGridUp[2] * (1. + f * localCubePhi[1]);
                    // 010
                    tii[2 + i * 8] += w * weightScalarGridDown[0]
                                      * weightScalarGridUp[1]
                                      * weightScalarGridDown[2]
                                      * (1. + f * localCubePhi[2]);
                    // 011
                    tii[3 + i * 8]
                        += w * weightScalarGridDown[0] * weightScalarGridUp[1]
                           * weightScalarGridUp[2] * (1. + f * localCubePhi[3]);
                    // 100
                    tii[4 + i * 8] += w * weightScalarGridUp[0]
                                      * weightScalarGridDown[1]
                                      * weightScalarGridDown[2]
                                      * (1. + f * localCubePhi[4]);
                    // 101
                    tii[5 + i * 8]
                        += w * weightScalarGridUp[0] * weightScalarGridDown[1]
                           * weightScalarGridUp[2] * (1. + f * localCubePhi[5]);
                    // 110
                    tii[6 + i * 8] += w * weightScalarGridUp[0]
                                      * weightScalarGridUp[1]
                                      * weightScalarGridDown[2]
                                      * (1. + f * localCubePhi[6]);
                    // 111
                    tii[7 + i * 8]
                        += w * weightScalarGridUp[0] * weightScalarGridUp[1]
                           * weightScalarGridUp[2] * (1. + f * localCubePhi[7]);
                }

                w = mass * q[0] * q[1] / e;
                tij[0] += w * weightScalarGridDown[2]
                          * (1.
                             + f * 0.25
                                   * (localCubePhi[0] + localCubePhi[2]
                                      + localCubePhi[4] + localCubePhi[6]));
                tij[1] += w * weightScalarGridUp[2]
                          * (1.
                             + f * 0.25
                                   * (localCubePhi[1] + localCubePhi[3]
                                      + localCubePhi[5] + localCubePhi[7]));

                w = mass * q[0] * q[2] / e;
                tij[2] += w * weightScalarGridDown[1]
                          * (1.
                             + f * 0.25
                                   * (localCubePhi[0] + localCubePhi[1]
                                      + localCubePhi[4] + localCubePhi[5]));
                tij[3] += w * weightScalarGridUp[1]
                          * (1.
                             + f * 0.25
                                   * (localCubePhi[2] + localCubePhi[3]
                                      + localCubePhi[6] + localCubePhi[7]));

                w = mass * q[1] * q[2] / e;
                tij[4] += w * weightScalarGridDown[0]
                          * (1.
                             + f * 0.25
                                   * (localCubePhi[0] + localCubePhi[1]
                                      + localCubePhi[2] + localCubePhi[3]));
                tij[5] += w * weightScalarGridUp[0]
                          * (1.
                             + f * 0.25
                                   * (localCubePhi[4] + localCubePhi[5]
                                      + localCubePhi[6] + localCubePhi[7]));
            }

            for (int i = 0; i < 3; i++)
                (*Tij) (xTij, i, i) += tii[8 * i];
            (*Tij) (xTij, 0, 1) += tij[0];
            (*Tij) (xTij, 0, 2) += tij[2];
            (*Tij) (xTij, 1, 2) += tij[4];

            for (int i = 0; i < 3; i++)
                (*Tij) (xTij + 0, i, i) += tii[4 + 8 * i];
            (*Tij) (xTij + 0, 1, 2) += tij[5];

            for (int i = 0; i < 3; i++)
                (*Tij) (xTij + 1, i, i) += tii[2 + 8 * i];
            (*Tij) (xTij + 1, 0, 2) += tij[3];

            for (int i = 0; i < 3; i++)
                (*Tij) (xTij + 2, i, i) += tii[1 + 8 * i];
            (*Tij) (xTij + 2, 0, 1) += tij[1];

            for (int i = 0; i < 3; i++)
                (*Tij) (xTij + 0 + 1, i, i) += tii[6 + 8 * i];
            for (int i = 0; i < 3; i++)
                (*Tij) (xTij + 0 + 2, i, i) += tii[5 + 8 * i];
            for (int i = 0; i < 3; i++)
                (*Tij) (xTij + 1 + 2, i, i) += tii[3 + 8 * i];
            for (int i = 0; i < 3; i++)
                (*Tij) (xTij + 0 + 1 + 2, i, i) += tii[7 + 8 * i];
        }
    }
}

#ifndef projection_Tij_comm
#define projection_Tij_comm symtensorProjectionCICNGP_comm
#endif

//////////////////////////
// projection_Ti0_project
//////////////////////////
// Description:
//   Particle-mesh projection for Ti0, including geometric corrections
//
// Arguments:
//   pcls       pointer to particle handler
//   Ti0        pointer to target field
//   phi        pointer to Bardeen potential which characterizes the
//              geometric corrections (volume distortion); can be set to
//              NULL which will result in no corrections applied
//   chi        pointer to difference between the Bardeen potentials which
//              characterizes additional corrections; can be set to
//              NULL which will result in no corrections applied
//   coeff      coefficient applied to the projection operation (default 1)
//
// Returns:
//
//////////////////////////

template <typename part, typename part_info, typename part_dataType>
void projection_Ti0_project (Particles<part, part_info, part_dataType> *pcls,
                             Field<Real> *Ti0, Field<Real> *phi = NULL,
                             Field<Real> *chi = NULL, double coeff = 1.)
{
    if (Ti0->lattice ().halo () == 0)
    {
        std::cout << "projection_Ti0_project: target field needs halo > 0"
                  << std::endl;
        std::exit (-1);
    }

    Site xPart (pcls->lattice ());
    Site xField (Ti0->lattice ());

    Real referPos[3];
    Real weightScalarGridUp[3];
    Real weightScalarGridDown[3];
    Real dx = pcls->res ();

    double mass = coeff / (dx * dx * dx);
    mass *= pcls->parts_info ()->mass;

    Real localCube[24]; // XYZ = 000 | 001 | 010 | 011 | 100 | 101 | 110 | 111
    Real localCubePhi[8];
    Real localCubeChi[8];

    for (int i = 0; i < 8; i++)
        localCubePhi[i] = 0.0;
    for (int i = 0; i < 8; i++)
        localCubeChi[i] = 0.0;

    for (xPart.first (), xField.first (); xPart.test ();
         xPart.next (), xField.next ())
    {
        if (pcls->field () (xPart).size != 0)
        {
            for (int i = 0; i < 3; i++)
                referPos[i] = xPart.coord (i) * dx;
            for (int i = 0; i < 24; i++)
                localCube[i] = 0.0;

            if (phi != NULL)
            {
                localCubePhi[0] = (*phi) (xField);
                localCubePhi[1] = (*phi) (xField + 2);
                localCubePhi[2] = (*phi) (xField + 1);
                localCubePhi[3] = (*phi) (xField + 1 + 2);
                localCubePhi[4] = (*phi) (xField + 0);
                localCubePhi[5] = (*phi) (xField + 0 + 2);
                localCubePhi[6] = (*phi) (xField + 0 + 1);
                localCubePhi[7] = (*phi) (xField + 0 + 1 + 2);
            }
            if (chi != NULL)
            {
                localCubeChi[0] = (*chi) (xField);
                localCubeChi[1] = (*chi) (xField + 2);
                localCubeChi[2] = (*chi) (xField + 1);
                localCubeChi[3] = (*chi) (xField + 1 + 2);
                localCubeChi[4] = (*chi) (xField + 0);
                localCubeChi[5] = (*chi) (xField + 0 + 2);
                localCubeChi[6] = (*chi) (xField + 0 + 1);
                localCubeChi[7] = (*chi) (xField + 0 + 1 + 2);
            }

            for (const auto& p : pcls->field () (xPart).parts)
            {
                for (int i = 0; i < 3; i++)
                {
                    weightScalarGridUp[i] = (p.pos[i] - referPos[i]) / dx;
                    weightScalarGridDown[i] = 1.0l - weightScalarGridUp[i];
                }

                const auto &q = p.vel;

                for (int i = 0; i < 3; i++)
                {
                    // 000
                    localCube[8 * i]
                        += weightScalarGridDown[0] * weightScalarGridDown[1]
                           * weightScalarGridDown[2] * q[i]
                           * (1. + 6 * localCubePhi[0] - localCubeChi[0]);
                    // 001
                    localCube[8 * i + 1]
                        += weightScalarGridDown[0] * weightScalarGridDown[1]
                           * weightScalarGridUp[2] * q[i]
                           * (1. + 6 * localCubePhi[1] - localCubeChi[1]);
                    // 010
                    localCube[8 * i + 2]
                        += weightScalarGridDown[0] * weightScalarGridUp[1]
                           * weightScalarGridDown[2] * q[i]
                           * (1. + 6 * localCubePhi[2] - localCubeChi[2]);
                    // 011
                    localCube[8 * i + 3]
                        += weightScalarGridDown[0] * weightScalarGridUp[1]
                           * weightScalarGridUp[2] * q[i]
                           * (1. + 6 * localCubePhi[3] - localCubeChi[3]);
                    // 100
                    localCube[8 * i + 4]
                        += weightScalarGridUp[0] * weightScalarGridDown[1]
                           * weightScalarGridDown[2] * q[i]
                           * (1. + 6 * localCubePhi[4] - localCubeChi[4]);
                    // 101
                    localCube[8 * i + 5]
                        += weightScalarGridUp[0] * weightScalarGridDown[1]
                           * weightScalarGridUp[2] * q[i]
                           * (1. + 6 * localCubePhi[5] - localCubeChi[5]);
                    // 110
                    localCube[8 * i + 6]
                        += weightScalarGridUp[0] * weightScalarGridUp[1]
                           * weightScalarGridDown[2] * q[i]
                           * (1. + 6 * localCubePhi[6] - localCubeChi[6]);
                    // 111
                    localCube[8 * i + 7]
                        += weightScalarGridUp[0] * weightScalarGridUp[1]
                           * weightScalarGridUp[2] * q[i]
                           * (1. + 6 * localCubePhi[7] - localCubeChi[7]);
                }
            }
            for (int i = 0; i < 3; i++)
            {
                (*Ti0) (xField, i) += localCube[8 * i] * mass;
                (*Ti0) (xField + 2, i) += localCube[8 * i + 1] * mass;
                (*Ti0) (xField + 1, i) += localCube[8 * i + 2] * mass;
                (*Ti0) (xField + 1 + 2, i) += localCube[8 * i + 3] * mass;
                (*Ti0) (xField + 0, i) += localCube[8 * i + 4] * mass;
                (*Ti0) (xField + 0 + 2, i) += localCube[8 * i + 5] * mass;
                (*Ti0) (xField + 0 + 1, i) += localCube[8 * i + 6] * mass;
                (*Ti0) (xField + 0 + 1 + 2, i) += localCube[8 * i + 7] * mass;
            }
        }
    }
}

//////////////////////////
// projectFTtheta
//////////////////////////
// Description:
//   Compute the diverge of the velocity in Fourier space
//
// Arguments:
//   thFT       reference to the Fourier image of the divergence of the
//   velocity field viFT       reference to the Fourier image of the velocity
//   field
//
// Returns:
//
//////////////////////////

void projectFTtheta (Field<Cplx> &thFT, Field<Cplx> &viFT);

//////////////////////////
// projectFTomega
//////////////////////////
// Description:
//   Compute the curl part of the velocity field in Fourier space
//
// Arguments:
//   viFT      reference to the input Fourier image of the velocity field
//             the divergence part will be projected out
//
// Returns:
//
//////////////////////////

void projectFTomega (Field<Cplx> &viFT);
}
#endif
