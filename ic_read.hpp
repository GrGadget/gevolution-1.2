//////////////////////////
// ic_read.hpp
//////////////////////////
//
// read initial conditions from disk
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen
// Mary University of London)
//
// Last modified: August 2019
//
//////////////////////////

#ifndef IC_READ_HEADER
#define IC_READ_HEADER

#include "LATfield2.hpp"
#include "Particles_gevolution.hpp" // Particles_gevolution
#include "metadata.hpp"             // metadata
#include "tools.hpp"                // Cplx
#include <iostream>
#include <set>
#include <string>

namespace gevolution
{
using LATfield2::FFT_BACKWARD;
using LATfield2::FFT_FORWARD;
using LATfield2::fileDsc;
using LATfield2::get_fileDsc_global;
using LATfield2::get_fileDsc_local;
using LATfield2::parallel;
using LATfield2::part_simple;
using LATfield2::part_simple_dataType;
using LATfield2::part_simple_info;
using LATfield2::PlanFFT;
using LATfield2::Real;

//////////////////////////
// readIC
//////////////////////////
// Description:
//   reads initial conditions from disk
//
// Arguments:
//   sim            simulation metadata structure
//   ic             settings for IC generation
//   cosmo          cosmological parameter structure
//   fourpiG        4 pi G (in code units)
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
//   cycle          reference to cycle counter (for restart after hibernation)
//   snapcount      reference to snapshot counter (for restart after
//   hibernation) pkcount        reference to spectra counter (for restart
//   after hibernation) restartcount   reference to restart counter (for
//   restart after hibernation)
//
// Returns:
//
//////////////////////////

void readIC (
    metadata &sim, icsettings &ic, cosmology &cosmo, const double fourpiG,
    double &a, double &tau, double &dtau, double &dtau_old,
    Particles_gevolution<part_simple, part_simple_info, part_simple_dataType>
        *pcls_cdm,
    Particles_gevolution<part_simple, part_simple_info, part_simple_dataType>
        *pcls_b,
    Particles_gevolution<part_simple, part_simple_info, part_simple_dataType>
        *pcls_ncdm,
    double *maxvel, Field<Real> *phi, Field<Real> *chi, Field<Real> *Bi,
    Field<Real> *source, Field<Real> *Sij, Field<Cplx> *scalarFT,
    Field<Cplx> *BiFT, Field<Cplx> *SijFT, PlanFFT<Cplx> *plan_phi,
    PlanFFT<Cplx> *plan_chi, PlanFFT<Cplx> *plan_Bi, PlanFFT<Cplx> *plan_source,
    PlanFFT<Cplx> *plan_Sij, int &cycle, int &snapcount, int &pkcount,
    int &restartcount, std::set<long> *IDbacklog);
}
#endif
