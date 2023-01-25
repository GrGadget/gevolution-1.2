/////////////////////////
// velocity.hpp
//////////////////////////
//
// Utilities for the velocity field
//
// Author: Goran Jelic-Cizmek (Université de Genève)
// Author: Francesca Lepori (SISSA Trieste & INFN Trieste & Université de
// Genève) Author: Julian Adamek (Queen Mary University of London)
//
// Last modified: April 2019
//
//////////////////////////

#ifndef VELOCITY_HEADER
#define VELOCITY_HEADER

#include "gevolution/config.h"
#include "LATfield2.hpp"
#include "gevolution/basic_types.hpp"
#include "gevolution/metadata.hpp"

namespace gevolution
{
using LATfield2::Field;
using LATfield2::Site;
/**
    computes and stores all the background functions
**/

//////////////////////////
// compute_vi_rescaled
//////////////////////////
// Description:
//   Compute the velocity field as v^i = T^i_0/T^0_0, if a = 1 then vi = a v^i
//   If T^0_0 = 0 the velocity field is set to be the one at the previous time
//   step, rescaled as v^i(a) = v^i(a_past) a*Hconf(a) dD1/da (velocity method
//   = rescaled) [see G. Jelic-Cizmek, F. Lepori, J. Adamek, and R. Durrer,
//   JCAP 1809, 006 (2018)]
//
// Arguments:
//   cosmo      structure containing the cosmological parameters
//   vi         reference to the velocity field, contains (a^3 T^i_0)
//   source     reference to the field source (a^4 T^0_0)
//   Ti0        reference to the field Ti0 (a^3 T^i_0)
//   a          scale factor at current time step
//   a_old      scale factor at previous time step
//
// Returns:
//
//////////////////////////

void compute_vi_rescaled (cosmology &cosmo, Field<Real> *vi,
                          Field<Real> *source, Field<Real> *Ti0, double a = 1.,
                          double a_old = 1.);
}
#endif
