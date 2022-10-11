//////////////////////////
// tools.hpp
//////////////////////////
//
// Collection of analysis tools for gevolution
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen
// Mary University of London)
//
// Last modified: April 2019
//
//////////////////////////

#ifndef TOOLS_HEADER
#define TOOLS_HEADER

#include "gevolution/config.h"
#include "LATfield2.hpp"
#include "gevolution/basic_types.hpp"
#include "gevolution/metadata.hpp"
#include <iostream>
#include <string>

namespace gevolution
{
using LATfield2::Field;
using LATfield2::parallel;
using LATfield2::rKSite;
using LATfield2::Site;


#define KTYPE_GRID 0
#define KTYPE_LINEAR 1

#ifdef FFT3D
//////////////////////////
// extractCrossSpectrum
//////////////////////////
// Description:
//   generates the cross spectrum for two Fourier images
//
// Arguments:
//   fld1FT     reference to the first Fourier image for which the cross
//   spectrum should be extracted fld2FT     reference to the second Fourier
//   image for which the cross spectrum should be extracted kbin allocated
//   array that will contain the central k-value for the bins power allocated
//   array that will contain the average power in each bin kscatter   allocated
//   array that will contain the k-scatter for each bin pscatter   allocated
//   array that will contain the scatter in power for each bin occupation
//   allocated array that will count the number of grid points contributing to
//   each bin numbin     number of bins (minimum size of all arrays) ktype flag
//   indicating which definition of momentum to be used
//                  0: grid momentum
//                  1: linear (default)
//   comp1      for component-wise cross spectra, the component for the first
//   field (ignored if negative) comp2      for component-wise cross spectra,
//   the component for the second field (ignored if negative)
//
// Returns:
//
//////////////////////////

void extractCrossSpectrum (Field<Cplx> &fld1FT, Field<Cplx> &fld2FT, Real *kbin,
                           Real *power, Real *kscatter, Real *pscatter,
                           int *occupation, const int numbins,
                           const bool deconvolve = true,
                           const int ktype = KTYPE_LINEAR, const int comp1 = -1,
                           const int comp2 = -1);

//////////////////////////
// extractPowerSpectrum
//////////////////////////
// Description:
//   generates the power spectrum for a Fourier image
//
// Arguments:
//   fldFT      reference to the Fourier image for which the power spectrum
//   should be extracted kbin       allocated array that will contain the
//   central k-value for the bins power      allocated array that will contain
//   the average power in each bin kscatter   allocated array that will contain
//   the k-scatter for each bin pscatter   allocated array that will contain
//   the scatter in power for each bin occupation allocated array that will
//   count the number of grid points contributing to each bin numbin     number
//   of bins (minimum size of all arrays) ktype      flag indicating which
//   definition of momentum to be used
//                  0: grid momentum
//                  1: linear (default)
//
// Returns:
//
//////////////////////////

void extractPowerSpectrum (Field<Cplx> &fldFT, Real *kbin, Real *power,
                           Real *kscatter, Real *pscatter, int *occupation,
                           const int numbins, const bool deconvolve = true,
                           const int ktype = KTYPE_LINEAR);
#endif

//////////////////////////
// writePowerSpectrum
//////////////////////////
// Description:
//   writes power spectra as tabulated data into ASCII file
//
// Arguments:
//   kbin           array containing the central values of k for each bin
//   power          array containing the central values of P(k) for each bin
//   kscatter       array containing the statistical error on k for each bin
//   pscatter       array containing the statistical error on P(k) for each bin
//   occupation     array containing the number of k-modes contributing to each
//   bin numbins        total number of bins (length of the arrays) rescalek
//   unit conversion factor for k rescalep       unit conversion factor for
//   P(k) filename       output file name description    descriptive header a
//   scale factor for this spectrum z_target       target redshift for this
//   output (used only if EXACT_OUTPUT_REDSHIFTS is defined)
//
// Returns:
//
//////////////////////////

void writePowerSpectrum (Real *kbin, Real *power, Real *kscatter,
                         Real *pscatter, int *occupation, const int numbins,
                         const Real rescalek, const Real rescalep,
                         const char *filename, const char *description,
                         double a, const double z_target = -1);

//////////////////////////
// computeVectorDiagnostics
//////////////////////////
// Description:
//   computes some diagnostics for the spin-1 perturbation
//
// Arguments:
//   Bi         reference to the real-space vector field to analyze
//   mdivB      will contain the maximum value of the divergence of Bi
//   mcurlB     will contain the maximum value of the curl of Bi
//
// Returns:
//
//////////////////////////

void computeVectorDiagnostics (Field<Real> &Bi, Real &mdivB, Real &mcurlB);

//////////////////////////
// computeTensorDiagnostics
//////////////////////////
// Description:
//   computes some diagnostics for the spin-2 perturbation
//
// Arguments:
//   hij        reference to the real-space tensor field to analyze
//   mdivh      will contain the maximum value of the divergence of hij
//   mtraceh    will contain the maximum value of the trace of hij
//   mnormh     will contain the maximum value of the norm of hij
//
// Returns:
//
//////////////////////////

void computeTensorDiagnostics (Field<Real> &hij, Real &mdivh, Real &mtraceh,
                               Real &mnormh);

//////////////////////////
// findIntersectingLightcones
//////////////////////////
// Description:
//   determines periodic copies of light cone vertex for which the present
//   look-back interval may overlap with a given spatial domain
//
// Arguments:
//   lightcone  reference to structure describing light cone geometry
//   outer      outer (far) limit of look-back interval
//   inner      inner (close) limit of look-back interval
//   domain     array of domain boundaries
//   vertex     will contain array of relevant vertex locations
//
// Returns:
//   number of vertices found
//
//////////////////////////

int findIntersectingLightcones (lightcone_geometry &lightcone, double outer,
                                double inner, double *domain,
                                double vertex[MAX_INTERSECTS][3]);

//////////////////////////
// hourMinSec
//////////////////////////
// Description:
//   generates formatted output for cpu-time: hh..h:mm:ss.s
//
// Arguments:
//   seconds    number of seconds
//
// Returns:
//   formatted string
//
//////////////////////////

std::string hourMinSec (double seconds);
}
#endif
