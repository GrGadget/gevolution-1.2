//////////////////////////
// hibernation.hpp
//////////////////////////
//
// Auxiliary functions for hibernation
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen
// Mary University of London)
//
// Last modified: February 2019
//
//////////////////////////

#ifndef HIBERNATION_HEADER
#define HIBERNATION_HEADER

#include "gevolution/config.h"
#include "LATfield2.hpp"
#include "gevolution/real_type.hpp"
#include "gevolution/metadata.hpp"
#include "gevolution/Particles_gevolution.hpp"

namespace gevolution
{
using LATfield2::Field;
using LATfield2::Particles;

//////////////////////////
// writeRestartSettings
//////////////////////////
// Description:
//   writes a settings file containing all the relevant metadata for restarting
//   a run from a hibernation point
//
// Arguments:
//   sim            simulation metadata structure
//   ic             settings for IC generation
//   cosmo          cosmological parameter structure
//   a              scale factor
//   tau            conformal coordinate time
//   dtau           time step
//   cycle          current main control loop cycle count
//   restartcount   restart counter aka number of hibernation point (default
//   -1)
//                  if < 0 no number is associated to the hibernation point
//
// Returns:
//
//////////////////////////

void writeRestartSettings (metadata &sim, icsettings &ic, cosmology &cosmo,
                           const double a, const double tau, const double dtau,
                           const int cycle, const int restartcount = -1);

//////////////////////////
// hibernate
//////////////////////////
// Description:
//   creates a hibernation point by writing snapshots of the simulation data
//   and metadata
//
// Arguments:
//   sim            simulation metadata structure
//   ic             settings for IC generation
//   cosmo          cosmological parameter structure
//   pcls_cdm       pointer to particle handler for CDM
//   pcls_b         pointer to particle handler for baryons
//   pcls_ncdm      array of particle handlers for non-cold DM
//   phi            reference to field containing first Bardeen potential
//   chi            reference to field containing difference of Bardeen
//   potentials Bi             reference to vector field containing
//   frame-dragging potential a              scale factor tau conformal
//   coordinate time dtau           time step cycle          current main
//   control loop cycle count restartcount   restart counter aka number of
//   hibernation point (default -1)
//                  if < 0 no number is associated to the hibernation point
//
// Returns:
//
//////////////////////////

void hibernate (
    metadata &sim, icsettings &ic, cosmology &cosmo,
    Particles_gevolution *pcls_cdm,
    Particles_gevolution *pcls_b,
    Particles_gevolution *pcls_ncdm,
    Field<Real> &phi, Field<Real> &chi, Field<Real> &Bi, const double a,
    const double tau, const double dtau, const int cycle,
    const int restartcount = -1);
}
#endif
