//////////////////////////
// radiation.hpp
//////////////////////////
//
// code components related to radiation and linear relativistic species
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen
// Mary University of London)
//
// Last modified: November 2019
//
//////////////////////////

#ifndef RADIATION_HEADER
#define RADIATION_HEADER

#include "gevolution/config.h"
#include "gevolution/real_type.hpp"

#ifdef HAVE_CLASS

namespace gevolution
{
//////////////////////////
// projection_T00_project (radiation module)
//////////////////////////
// Description:
//   provides a realization of the linear density field of radiation and
//   non-cold species using linear transfer functions precomputed with CLASS;
//   the contributions for the various species are included only until some
//   individual redshift values are reached (after which no linear treatment
//   is requested)
//
// Arguments:
//   class_background  CLASS structure that contains the background
//   class_perturbs    CLASS structure that contains the perturbations
//   source            reference to field that will contain the realization
//   scalarFT          reference to Fourier image of that field
//   plan_source       pointer to FFT planner
//   sim               simulation metadata structure
//   ic                settings for IC generation (contains the random seed)
//   cosmo             cosmological parameter structure
//   fourpiG           4 pi G (in code units)
//   a                 scale factor
//   coeff             multiplicative coefficient (default 1)
//
// Returns:
//
//////////////////////////

void projection_T00_project (background &class_background,
                             perturbs &class_perturbs, Field<Real> &source,
                             Field<Cplx> &scalarFT, PlanFFT<Cplx> *plan_source,
                             metadata &sim, icsettings &ic, cosmology &cosmo,
                             const double fourpiG, double a, double coeff = 1.);

//////////////////////////
// prepareFTchiLinear
//////////////////////////
// Description:
//   provides a (Fourier-space) realization of chi (generated by radiation and
//   non-cold species) from the linear transfer functions precomputed with
//   CLASS
//
// Arguments:
//   class_background  CLASS structure that contains the background
//   class_perturbs    CLASS structure that contains the perturbations
//   scalarFT          reference to Fourier image of field; will contain the
//                     (Fourier image) of the realization
//   sim               simulation metadata structure
//   ic                settings for IC generation (contains the random seed)
//   cosmo             cosmological parameter structure
//   fourpiG           4 pi G (in code units)
//   a                 scale factor
//   coeff             multiplicative coefficient (default 1)
//
// Returns:
//
//////////////////////////

void prepareFTchiLinear (background &class_background, perturbs &class_perturbs,
                         Field<Cplx> &scalarFT, metadata &sim, icsettings &ic,
                         cosmology &cosmo, const double fourpiG, double a,
                         double coeff = 1.);

}
#endif

#endif