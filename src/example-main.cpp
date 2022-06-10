//////////////////////////
// Copyright (c) 2022 Eduardo Quintana Miranda
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESSED OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
//////////////////////////

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
namespace mpi = boost::mpi;

#include "gevolution/gevolution.hpp"
#include "gevolution/newtonian_pm.hpp"
#include "gevolution/gr_pm.hpp"

#include "gevolution/Particles_gevolution.hpp"
#include "gevolution/background.hpp"
#include "gevolution/class_tools.hpp"
#include "gevolution/ic_basic.hpp"
#include "gevolution/ic_read.hpp"
#include "gevolution/metadata.hpp"
#include "gevolution/tools.hpp"
#ifdef ICGEN_PREVOLUTION
#include "gevolution/ic_prevolution.hpp"
#endif
#include "gevolution/hibernation.hpp"
#include "gevolution/output.hpp"
#include "gevolution/parser.hpp"
#include "gevolution/radiation.hpp"
#ifdef VELOCITY
#include "gevolution/velocity.hpp"
#endif

#include <filesystem>
namespace fs = std::filesystem;
using namespace std;
using namespace LATfield2;
using namespace gevolution;

int main (int argc, char **argv)
{
    mpi::environment env;
    mpi::communicator com_world;
    
    int n = 0, m = 0;

    int cycle = 0, snapcount = 0, pkcount = 0, restartcount = 0,
              usedparams, numparam = 0;
    
    int box[3];
    double dtau, dtau_old, dx, tau, a, tmp;
    char filename[2 * PARAM_MAX_LENGTH + 24];
    char *settingsfile = NULL;

    parameter *params = NULL;
    metadata sim;
    cosmology cosmo;
    icsettings ic;
    
    double maxvel[MAX_PCL_SPECIES];

    for (int i = 1; i < argc; i++)
    {
        if (argv[i][0] != '-')
            continue;
        switch (argv[i][1])
        {
        case 's':
            settingsfile = argv[++i]; // settings file name
            break;
        case 'n':
            n = atoi (argv[++i]); // size of the dim 1 of the processor grid
            break;
        case 'm':
            m = atoi (argv[++i]); // size of the dim 2 of the processor grid
            break;
        }
    }

    parallel.initialize (com_world, n, m);

    COUT << "  _   _      _         __ ,  _" << endl;
    COUT << " (_| (-' \\/ (_) (_ (_| (  ( (_) /\\/	version 1.2         "
            "running on "
         << n * m << " cores." << endl;
    COUT << "  -'" << endl << endl;

    if (settingsfile == NULL)
    {
        COUT << " error"
             << ": no settings file specified!" << endl;
        parallel.abortForce ();
    }

    COUT << " initializing..." << endl;

    numparam = loadParameterFile (settingsfile, params);

    usedparams = parseMetadata (params, numparam, sim, cosmo, ic);

    COUT << " parsing of settings file completed. " << numparam
         << " parameters found, " << usedparams << " were used." << endl;

    sprintf (filename, "%s%s_settings_used.ini", sim.output_path,
             sim.basename_generic);
    saveParameterFile (filename, params, numparam);

    free (params);

        numparam = 0;

    box[0] = sim.numpts;
    box[1] = sim.numpts;
    box[2] = sim.numpts;

    Lattice lat (3, box, 2);
    Lattice latFT;
    latFT.initializeRealFFT (lat, 0);

    Particles_gevolution
        pcls_cdm,pcls_b,pcls_ncdm[MAX_PCL_SPECIES-2];
    set<long> IDbacklog[MAX_PCL_SPECIES];



    Site x (lat);
    rKSite kFT (latFT);

    dx = 1.0 / (double)sim.numpts;
    // numpts3d = (long)sim.numpts * (long)sim.numpts * (long)sim.numpts;

    for (int i = 0; i < 3;
         i++) // particles may never move farther than to the adjacent domain
    {
        if (lat.sizeLocal (i) - 1 < sim.movelimit)
            sim.movelimit = lat.sizeLocal (i) - 1;
    }
    parallel.min (sim.movelimit);

    cosmo.fourpiG
        = 1.5 * sim.boxsize * sim.boxsize / cosmo.C_SPEED_OF_LIGHT / cosmo.C_SPEED_OF_LIGHT;
    a = 1. / (1. + sim.z_in);
    tau = particleHorizon (a, cosmo);
    unique_ptr<debugger_t> Debugger_ptr{
        Debugger = new
        debugger_t(
            com_world,
            "forcetest.bin",
            1000 * sim.boxsize,
            10 * cosmo.C_SPEED_OF_LIGHT * cosmo.C_SPEED_OF_LIGHT/sim.boxsize)};

    dtau = std::min(sim.Cf * dx, sim.steplimit/Hconf(a,cosmo));
    dtau_old = dtau;

    {
    Field<Real> phi;
    Field<Real> source;
    Field<Real> chi;
    Field<Real> Sij;
    Field<Real> Bi;
    Field<Cplx> scalarFT;
    Field<Cplx> SijFT;
    Field<Cplx> BiFT;
    source.initialize (lat, 1);
    phi.initialize (lat, 1);
    chi.initialize (lat, 1);
    scalarFT.initialize (latFT, 1);
    PlanFFT<Cplx> plan_source (&source, &scalarFT);
    PlanFFT<Cplx> plan_phi (&phi, &scalarFT);
    PlanFFT<Cplx> plan_chi (&chi, &scalarFT);
    Sij.initialize (lat, 3, 3, matrix_symmetry::symmetric);
    SijFT.initialize (latFT, 3, 3, matrix_symmetry::symmetric);
    PlanFFT<Cplx> plan_Sij (&Sij, &SijFT);
    Bi.initialize (lat, 3);
    BiFT.initialize (latFT, 3);
    PlanFFT<Cplx> plan_Bi (&Bi, &BiFT);
    if (ic.generator == ICGEN_BASIC)
        generateIC_basic (sim, ic, cosmo, &pcls_cdm, &pcls_b,
                          pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij,
                          &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi,
                          &plan_Bi, &plan_source, &plan_Sij, params,
                          numparam); // generates ICs on the fly
    else if (ic.generator == ICGEN_READ_FROM_DISK)
        readIC (sim, ic, cosmo, a, tau, dtau, dtau_old, &pcls_cdm,
                &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij,
                &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi,
                &plan_source, &plan_Sij, cycle, snapcount, pkcount,
                restartcount, IDbacklog);
#ifdef ICGEN_PREVOLUTION
    else if (ic.generator == ICGEN_PREVOLUTION)
        generateIC_prevolution (sim, ic, cosmo, a, tau, dtau, dtau_old,
                                &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi,
                                &chi, &Bi, &source, &Sij, &scalarFT, &BiFT,
                                &SijFT, &plan_phi, &plan_chi, &plan_Bi,
                                &plan_source, &plan_Sij, params, numparam);
#endif
    else
    {
        COUT << " error: IC generator not implemented!" << endl;
        parallel.abortForce ();
    }
    }
    
    relativistic_pm<Cplx,Particles_gevolution> PM(sim.numpts,com_world);

    do // main loop
    {
        COUT << "Cycle: " << cycle << std::endl;
        // PM step 1. construction of the energy momentum tensor
        PM.clear_sources();
        PM.sample(pcls_cdm,a);
        
        // PM step 2. compute the potentials
        PM.compute_potential(
            cosmo.fourpiG, 
            a,
            Hconf(a,cosmo),
            cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm (a, cosmo));
        
        PM.save_power_spectrum("power"+std::to_string(pkcount));
        
        PM.compute_forces(pcls_cdm,1.0,a);
        
        // Kick
        pcls_cdm.for_each(
            [&]
            (particle& part, const Site& /*xpart*/)
            {
               const double dtau_eff =  
                               (dtau + dtau_old) * 0.5 ;
               for(int i=0;i<3;++i)
               {
                   part.momentum[i] += dtau_eff * part.force[i];
               }
            }
            );
        
        rungekutta4bg (a, cosmo,
                       0.5 * dtau); // evolve background by half a time step
        
        PM.compute_velocities(pcls_cdm,a);
        
        // Drift
        pcls_cdm.for_each(
            [&](particle& part, const Site& /* xpart */)
            {
                for(int i=0;i<3;++i)
                {
                    part.pos[i] += dtau * part.vel[i];
                }
            }
        );
        
        // re-arrange particles in domains
        pcls_cdm.moveParticles();
        
        rungekutta4bg (a, cosmo,
                       0.5 * dtau); // evolve background by half a time step
        
        tau += dtau;
        dtau_old = dtau;
        dtau = std::min(sim.Cf,sim.steplimit/Hconf(a,cosmo));
        cycle++;

    }while(cycle<5);
    
    COUT << "Simulation complete" << std::endl;

    return 0;
}
