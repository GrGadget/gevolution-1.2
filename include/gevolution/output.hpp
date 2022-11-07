//////////////////////////
// output.hpp
//////////////////////////
//
// Output of snapshots, light cones and spectra
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen
// Mary University of London)
//
// Last modified: March 2020
//
//////////////////////////

#ifndef OUTPUT_HEADER
#define OUTPUT_HEADER

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <set>
#include <string>

#include "gevolution/config.h"
#include "LATfield2.hpp"
#include "gevolution/basic_types.hpp"
#include "gevolution/Particles_gevolution.hpp"
#include "gevolution/background.hpp"
#include "gevolution/gevolution.hpp"
#include "gevolution/metadata.hpp"
#include "gevolution/tools.hpp"

namespace gevolution
{
using LATfield2::FFT_BACKWARD;
using LATfield2::FFT_FORWARD;
using LATfield2::parallel;
using LATfield2::PlanFFT;
//////////////////////////
// writeSnapshots
//////////////////////////
// Description:
//   output of snapshots
//
// Arguments:
//   sim            simulation metadata structure
//   cosmo          cosmological parameter structure
//   a              scale factor
//   snapcount      snapshot index
//   h5filename     base name for HDF5 output file
//   pcls_cdm       pointer to particle handler for CDM
//   pcls_b         pointer to particle handler for baryons
//   pcls_ncdm      array of particle handlers for
//                  non-cold DM (may be set to NULL)
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
//   Bi_check       pointer to allocated field
//   BiFT_check     pointer to allocated field
//   plan_Bi_check  pointer to FFT planner
//   vi             pointer to allocated field
//
// Returns:
//
//////////////////////////

void write_snapshot (
    const metadata & sim, 
    const cosmology& cosmo,
    const double a /* scale factor */, 
    const Particles_gevolution & pcls,
    std::string filename);

// void writeSnapshots (
//     metadata &sim, const cosmology cosmo, const double a,
//     const double dtau_old, const int done_hij, const int snapcount,
//     std::string h5filename,
//     Particles_gevolution *pcls_cdm,
//     Particles_gevolution *pcls_b,
//     Particles_gevolution *pcls_ncdm,
//     Field<Real> *phi, Field<Real> *chi, Field<Real> *Bi, Field<Real> *source,
//     Field<Real> *Sij, Field<Cplx> *scalarFT, Field<Cplx> *BiFT,
//     Field<Cplx> *SijFT, PlanFFT<Cplx> *plan_phi, PlanFFT<Cplx> *plan_chi,
//     PlanFFT<Cplx> *plan_Bi, PlanFFT<Cplx> *plan_source, PlanFFT<Cplx> *plan_Sij
// #ifdef CHECK_B
//     ,
//     Field<Real> *Bi_check, Field<Cplx> *BiFT_check, PlanFFT<Cplx> *plan_Bi_check
// #endif
// #ifdef VELOCITY
//     ,
//     Field<Real> *vi
// #endif
// );

//////////////////////////
// writeLightcones
//////////////////////////
// Description:
//   output of light cones
//
// Arguments:
//   sim            simulation metadata structure
//   cosmo          cosmological parameter structure
//   a              scale factor
//   tau            conformal time
//   dtau           conformal time step
//   dtau_old       conformal time step of previous cycle
//   maxvel         maximum cdm velocity
//   cycle          current simulation cycle
//   h5filename     base name for HDF5 output file
//   pcls_cdm       pointer to particle handler for CDM
//   pcls_b         pointer to particle handler for baryons
//   pcls_ncdm      array of particle handlers for
//                  non-cold DM (may be set to NULL)
//   phi            pointer to allocated field
//   chi            pointer to allocated field
//   Bi             pointer to allocated field
//   Sij            pointer to allocated field
//   BiFT           pointer to allocated field
//   SijFT          pointer to allocated field
//   plan_Bi        pointer to FFT planner
//   plan_Sij       pointer to FFT planner
//   done_hij       reference to tensor projection flag
//   IDbacklog      IDs of particles written in previous cycle
//
// Returns:
//
//////////////////////////

void writeLightcones (
    metadata &sim, const cosmology cosmo, const double a,
    const double tau, const double dtau, const double dtau_old,
    const double maxvel, const int cycle, std::string h5filename,
    Particles_gevolution *pcls_cdm,
    Particles_gevolution *pcls_b,
    Particles_gevolution *pcls_ncdm,
    Field<Real> *phi, Field<Real> *chi, Field<Real> *Bi, Field<Real> *Sij,
    Field<Cplx> *BiFT, Field<Cplx> *SijFT, PlanFFT<Cplx> *plan_Bi,
    PlanFFT<Cplx> *plan_Sij, int &done_hij, std::set<long> *IDbacklog);

void writeLightcone_shell(
    const std::string filename,
    const gadget2_header hdr,
    const lightcone_geometry lc,
    const Particles_gevolution &Pcdm)
// Notice that lc describes the lightcone shell we want to save. The caller of this function is the
// one responsible to split a full lightcone into thin lightcone shells as well as construct the
// snapshot header template.
{
    auto particles_in_lightcone = [&](const gevolution::particle&, const LATfield2::Site&)
    {
        // FIXME: select the particles that lie in the lightcone
        return false;
    };
    
    Pcdm.saveGadget2 (
        filename, 
        hdr,
        particles_in_lightcone);
}

//////////////////////////
// writeSpectra
//////////////////////////
// Description:
//   output of spectra
//
// Arguments:
//   sim            simulation metadata structure
//   cosmo          cosmological parameter structure
//   a              scale factor
//   pkcount        spectrum output index
//   pcls_cdm       pointer to particle handler for CDM
//   pcls_b         pointer to particle handler for baryons
//   pcls_ncdm      array of particle handlers for
//                  non-cold DM (may be set to NULL)
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
//   Bi_check       pointer to allocated field (or NULL)
//   BiFT_check     pointer to allocated field (or NULL)
//   plan_Bi_check  pointer to FFT planner (or NULL)
//
// Returns:
//
//////////////////////////

void writeSpectra (
    metadata &sim, const cosmology cosmo, const double a,
    const int pkcount,
#ifdef HAVE_CLASS
    background &class_background, perturbs &class_perturbs, icsettings &ic,
#endif
    Particles_gevolution *pcls_cdm,
    Particles_gevolution *pcls_b,
    Particles_gevolution *pcls_ncdm,
    Field<Real> *phi, Field<Real> *chi, Field<Real> *Bi, Field<Real> *source,
    Field<Real> *Sij, Field<Cplx> *scalarFT, Field<Cplx> *BiFT,
    Field<Cplx> *SijFT, PlanFFT<Cplx> *plan_phi, PlanFFT<Cplx> *plan_chi,
    PlanFFT<Cplx> *plan_Bi, PlanFFT<Cplx> *plan_source, PlanFFT<Cplx> *plan_Sij
#ifdef CHECK_B
    ,
    Field<Real> *Bi_check, Field<Cplx> *BiFT_check, PlanFFT<Cplx> *plan_Bi_check
#endif
#ifdef VELOCITY
    ,
    Field<Real> *vi, Field<Cplx> *viFT, PlanFFT<Cplx> *plan_vi
#endif
);
}
#endif
