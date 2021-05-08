//////////////////////////
// background.hpp
//////////////////////////
//
// code components related to background evolution
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen
// Mary University of London)
//
// Last modified: September 2018
//
//////////////////////////

#ifndef BACKGROUND_HEADER
#define BACKGROUND_HEADER

#include "metadata.hpp"
#include <gsl/gsl_integration.h>

namespace gevolution
{

// double FermiDiracIntegrand (double q, void *w);

//////////////////////////
// FermiDiracIntegral
//////////////////////////
// Description:
//   computes the integral of the relativistic Fermi-Dirac distribution
//
// Arguments:
//   w          parameter in the F-D distribution, "(m a / kB T)^2"
//
// Returns: value for the integral
//
//////////////////////////

double FermiDiracIntegral (double &w);

//////////////////////////
// bg_ncdm (1)
//////////////////////////
// Description:
//   computes the background model for one ncdm species by integrating the
//   relativistic Fermi-Dirac distribution
//
// Arguments:
//   a          scale factor at which to compute the background model
//   cosmo      structure containing the cosmological parameters
//   p          index of the ncdm species
//
// Returns: value for the background model
//
//////////////////////////

double bg_ncdm (const double a, const cosmology cosmo, const int p);

//////////////////////////
// bg_ncdm (2)
//////////////////////////
// Description:
//   computes the background model for all ncdm species by integrating the
//   relativistic Fermi-Dirac distribution
//
// Arguments:
//   a          scale factor at which to compute the background model
//   cosmo      structure containing the cosmological parameters
//
// Note:
//   For optimization, the last value of a is stored in a static variable such
//   that multiple calls at the same value of a will not result in multiple
//   integrations being carried out. This assumes that the cosmological model
//   should not change!
//
// Returns: value for the background model
//
//////////////////////////

double bg_ncdm (const double a, const cosmology cosmo);

//////////////////////////
// Hconf
//////////////////////////
// Description:
//   computes the conformal Hubble rate at given scale factor
//
// Arguments:
//   a          scale factor
//   fourpiG    "4 pi G"
//   cosmo      structure containing the cosmological parameters
//
// Returns: conformal Hubble rate
//
//////////////////////////

double Hconf (const double a, const double fourpiG, const cosmology cosmo);

double Omega_m (const double a, const cosmology cosmo);

double Omega_rad (const double a, const cosmology cosmo);

double Omega_Lambda (const double a, const cosmology cosmo);

//////////////////////////
// rungekutta4bg
//////////////////////////
// Description:
//   integrates the Friedmann equation for the background model using a
//   fourth-order Runge-Kutta method
//
// Arguments:
//   a          scale factor (will be advanced by dtau)
//   fourpiG    "4 pi G"
//   cosmo      structure containing the cosmological parameters
//   dtau       time step by which the scale factor should be advanced
//
// Returns:
//
//////////////////////////

void rungekutta4bg (double &a, const double fourpiG, const cosmology cosmo,
                    const double dtau);

// double particleHorizonIntegrand (double sqrta, void *cosmo);

//////////////////////////
// particleHorizon
//////////////////////////
// Description:
//   computes the particle horizon (tau) at given scale factor
//
// Arguments:
//   a          scale factor
//   fourpiG    "4 pi G"
//   cosmo      structure containing the cosmological parameters
//
// Returns: particle horizon (tau)
//
//////////////////////////

double particleHorizon (const double a, const double fourpiG, cosmology &cosmo);
}
#endif
