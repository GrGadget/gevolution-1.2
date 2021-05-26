//////////////////////////
// ic_prevolution.hpp
//////////////////////////
//
// initial condition generator for gevolution consistently
// generating the full phase space for relativistic particles
// [see C.-P. Ma and E. Bertschinger, Astrophys. J. 429, 22 (1994)]
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen
// Mary University of London)
//
// Last modified: November 2019
//
//////////////////////////

#ifndef IC_PREVOLUTION_HEADER
#define IC_PREVOLUTION_HEADER

#ifdef HAVE_CLASS

#include "LATfield2.hpp"
#include "tools.hpp"
namespace gevolution
{
//////////////////////////
// generateIC_prevolution
//////////////////////////
// Description:
//   initial condition generator consistently generating the full phase
//   space for relativistic particles
//
// Arguments:
//   sim            simulation metadata structure
//   ic             settings for IC generation
//   cosmo          cosmological parameter structure
//   a              reference to scale factor
//   tau            reference to conformal coordinate time
//   dtau           time step
//   dtau_old       previous time step (will be passed back)
//   pcls_cdm       pointer to (uninitialized) particle handler for CDM
//   pcls_b         pointer to (uninitialized) particle handler for baryons
//   pcls_ncdm      array of (uninitialized) particle handlers for
//                  non-cold DM (may be set to NULL)
//   maxvel         array that will contain the maximum q/m/a (max. velocity)
//   phi            pointer to allocated field
//   chi            pointer to allocated field
//   Bi             pointer to allocated field
//   source         pointer to allocated field
//   Sij            pointer to allocated field
//   scalarFT       pointer to allocated field
//   BiFT           pointer to allocated field
//   SijFT          pointer to allocated field
//   plan_phi       pointer to FFT planner
//   plan_chi       pointer to FFT planner
//   plan_Bi        pointer to FFT planner
//   plan_source    pointer to FFT planner
//   plan_Sij       pointer to FFT planner
//   params         pointer to array of precision settings for CLASS (can be
//   NULL) numparam       number of precision settings for CLASS (can be 0)
//
// Returns:
//
//////////////////////////

void generateIC_prevolution (
    metadata &sim, icsettings &ic, const cosmology cosmo,
    double &a, double &tau, double &dtau, double &dtau_old,
    Particles_gevolution *pcls_cdm,
    Particles_gevolution *pcls_b,
    Particles_gevolution *pcls_ncdm,
    double *maxvel, Field<Real> *phi, Field<Real> *chi, Field<Real> *Bi,
    Field<Real> *source, Field<Real> *Sij, Field<Cplx> *scalarFT,
    Field<Cplx> *BiFT, Field<Cplx> *SijFT, PlanFFT<Cplx> *plan_phi,
    PlanFFT<Cplx> *plan_chi, PlanFFT<Cplx> *plan_Bi, PlanFFT<Cplx> *plan_source,
    PlanFFT<Cplx> *plan_Sij, parameter *params, int &numparam);
}
#endif
#endif
