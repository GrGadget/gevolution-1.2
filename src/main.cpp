//////////////////////////
// Copyright (c) 2015-2019 Julian Adamek
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

//////////////////////////
// main.cpp
//////////////////////////
//
// main control sequence of Geneva N-body code with evolution of metric
// perturbations (gevolution)
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen
// Mary University of London)
//
// Last modified: November 2019
//
//////////////////////////

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
namespace mpi = boost::mpi;

#include "gevolution/gevolution.hpp"
#include "gevolution/newtonian_pm.hpp"
#include "gevolution/gr_pm.hpp"
#include "gevolution/debugger.hpp"
#include "version.h"
#include <set>
#include <stdlib.h>
#include <vector>
#ifdef HAVE_CLASS
#include "class.h"
#undef MAX // due to macro collision this has to be done BEFORE including
           // LATfield2 headers!
           // then do not use macros
#undef MIN
#endif
#include "LATfield2.hpp"
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
#ifdef ICGEN_FALCONIC
#include "fcn/togevolution.hpp"
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

// stop condition, by external file 'stop'
// in the execution directory
bool stop()
{
    fs::path p{"stop"};
    bool ret = false;
    
    if(fs::exists(p))
    {
        ret = true;
        fs::remove(p);
    }
    return ret;
}

void show_constant(double v, std::string name)
{
    COUT << "Constant " << name << " = " << v << "\n";
}


int main (int argc, char **argv)
{
    mpi::environment env;
    mpi::communicator com_world;
    
#ifdef BENCHMARK
    // benchmarking variables
    double ref_time, ref2_time, cycle_start_time;
    double initialization_time;
    double run_time;
    double cycle_time = 0;
    double projection_time = 0;
    double snapshot_output_time = 0;
    double spectra_output_time = 0;
    double lightcone_output_time = 0;
    double gravity_solver_time = 0;
    double fft_time = 0;
    int fft_count = 0;
    double update_q_time = 0;
    int update_q_count = 0;
    double moveParts_time = 0;
    int moveParts_count = 0;
#endif // BENCHMARK

    int n = 0, m = 0;
#ifdef EXTERNAL_IO
    int io_size = 0;
    int io_group_size = 0;
#endif

    int i, j, cycle = 0, snapcount = 0, pkcount = 0, restartcount = 0,
              usedparams, numparam = 0, numspecies, done_hij;
    int numsteps_ncdm[MAX_PCL_SPECIES - 2];
    long numpts3d;
    int box[3];
    double dtau, dtau_old, dx, tau, a, tmp, start_time;
    double maxvel[MAX_PCL_SPECIES];
    FILE *outfile;
    char filename[2 * PARAM_MAX_LENGTH + 24];
    string h5filename;
    char *settingsfile = NULL;

#ifdef HAVE_CLASS
    char *precisionfile = NULL;
#endif
    parameter *params = NULL;
    metadata sim;
    cosmology cosmo;
    icsettings ic;
    double T00hom;

#ifndef H5_DEBUG
    H5Eset_auto2 (H5E_DEFAULT, NULL, NULL);
#endif

    for (i = 1; i < argc; i++)
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
        case 'p':
            cout << "HAVE_CLASS needs to be set at compilation to use "
                    "CLASS "
                    "precision files"
                 << endl;
            exit (-100);
            break;
        case 'i':
            cout << "EXTERNAL_IO needs to be set at compilation to use "
                    "the I/O "
                    "server"
                 << endl;
            exit (-1000);
            break;
        case 'g':
            cout << "EXTERNAL_IO needs to be set at compilation to use "
                    "the I/O "
                    "server"
                 << endl;
            exit (-1000);
        }
    }

    parallel.initialize (com_world, n, m);

    COUT << COLORTEXT_WHITE << endl;
    COUT << "  _   _      _         __ ,  _" << endl;
    COUT << " (_| (-' \\/ (_) (_ (_| (  ( (_) /\\/	version 1.2         "
            "running on "
         << n * m << " cores." << endl;
    COUT << "  -'" << endl << COLORTEXT_RESET << endl;
    COUT << "Version date: " GIT_DATE "\n"
            "Commit: " GIT_COMMIT "\n\n";

    if (settingsfile == NULL)
    {
        COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET
             << ": no settings file specified!" << endl;
        parallel.abortForce ();
    }

    COUT << " initializing..." << endl;

    start_time = MPI_Wtime ();

    numparam = loadParameterFile (settingsfile, params);

    usedparams = parseMetadata (params, numparam, sim, cosmo, ic);

    COUT << " parsing of settings file completed. " << numparam
         << " parameters found, " << usedparams << " were used." << endl;

    sprintf (filename, "%s%s_settings_used.ini", sim.output_path,
             sim.basename_generic);
    saveParameterFile (filename, params, numparam);

    free (params);

        numparam = 0;

    h5filename.reserve (2 * PARAM_MAX_LENGTH);
    h5filename.assign (sim.output_path);

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
    numpts3d = (long)sim.numpts * (long)sim.numpts * (long)sim.numpts;

    for (i = 0; i < 3;
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
#ifdef ICGEN_FALCONIC
    else if (ic.generator == ICGEN_FALCONIC)
        maxvel[0] = generateIC_FalconIC (
            sim, ic, cosmo, dtau, &pcls_cdm, pcls_ncdm, maxvel + 1,
            &phi, &source, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT,
            &plan_phi, &plan_source, &plan_chi, &plan_Bi, &plan_source,
            &plan_Sij);
#endif
    else
    {
        COUT << " error: IC generator not implemented!" << endl;
        parallel.abortForce ();
    }
    }
    
    if (sim.baryon_flag > 1)
    {
        COUT << " error: baryon_flag > 1 after IC generation, something "
                "went "
                "wrong in IC generator!"
             << endl;
        parallel.abortForce ();
    }

    numspecies = 1 + sim.baryon_flag + cosmo.num_ncdm;
    parallel.max<double> (maxvel, numspecies);

    if (sim.gr_flag == gravity_theory::GR)
    {
        for (i = 0; i < numspecies; i++)
            maxvel[i] /= sqrt (maxvel[i] * maxvel[i] + 1.0);
    }


        COUT << COLORTEXT_GREEN << " initialization complete."
             << COLORTEXT_RESET << endl
             << endl;

    
    std::unique_ptr< particle_mesh<Cplx,Particles_gevolution> > PM;
    
    if(sim.gr_flag==gravity_theory::GR)
    {
        PM.reset(
            new relativistic_pm<Cplx,Particles_gevolution>(sim.numpts)
        );
    }else
    {
        PM.reset(
            new newtonian_pm<Cplx,Particles_gevolution>(sim.numpts)
        );
    }
    
    pcls_cdm.update_mass(); // fix the mass legacy problem
    
    do // main loop
    {
        COUT << "Starting cycle: " << cycle << '\n';        
        
        // PM step 1. construction of the energy momentum tensor
        PM->clear_sources();
        PM->sample(pcls_cdm,a);
        
        for(auto report = PM->report();;)
        {
            COUT << report << "\n";
            break;
        }
        
        // EM tensor
        T00hom = PM->density();
        if (cycle % CYCLE_INFO_INTERVAL == 0)
        {
            COUT << " cycle " << cycle
                 << ", background information: z = " << (1. / a) - 1.
                 << ", average T00 = " << T00hom << ", background model = "
                 << cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm (a, cosmo)
                 << endl;
            
            std::array<double,3> mpv;
            mpv = PM->test_velocities(pcls_cdm);
            
            COUT << " mean     mass: " << mpv[0] << "\n";
            COUT << " mean sqr(pos): " << mpv[1] << "\n";
            COUT << " mean sqr(vel): " << mpv[2] << "\n";
        }
         
        // PM step 2. compute the potentials
        PM->compute_potential(cosmo.fourpiG, a,Hconf(a,cosmo),dtau_old,
            cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm (a, cosmo));
        
        // Sources
        for(auto report = PM->report();;)
        {
            COUT << report << "\n";
            break;
        }

        // record some background data
        // if (kFT.setCoord (0, 0, 0))
        // {
        //     sprintf (filename, "%s%s_background.dat", sim.output_path,
        //              sim.basename_generic);
        //     outfile = fopen (filename, "a");
        //     if (outfile == NULL)
        //     {
        //         cout << " error opening file for background output!" << endl;
        //     }
        //     else
        //     {
        //         if (cycle == 0)
        //             fprintf (outfile, "# background statistics\n# cycle   "
        //                               "tau/boxsize    a      "
        //                               "       conformal H/H0  phi(k=0)      "
        //                               " T00(k=0)\n");
        //         fprintf (outfile, " %6d   %e   %e   %e   %e   %e\n", cycle, tau,
        //                  a,
        //                  Hconf (a, cosmo) / Hconf (1., cosmo),
        //                  PM->phi_FT (kFT).real (), T00hom);
        //         fclose (outfile);
        //     }
        // }
        // done recording background data

        // lightcone output
        if (sim.num_lightcone > 0)
            //writeLightcones (sim, cosmo, a, tau, dtau, dtau_old,
            //                 maxvel[0], cycle,
            //                 h5filename + sim.basename_lightcone, &pcls_cdm,
            //                 &pcls_b, pcls_ncdm, &PM->phi, &PM->chi, &PM->Bi,
            //                 &PM->Tij, &PM->Bi_FT,
            //                 &PM->Tij_FT, &PM->plan_Bi, &PM->plan_Tij, done_hij, IDbacklog);
            ;
        else
            done_hij = 0;


        // TODO: snapshot output
        if (snapcount < sim.num_snapshot
            && 1. / a < sim.z_snapshot[snapcount] + 1.)
        {
            COUT << COLORTEXT_CYAN << " writing snapshot" << COLORTEXT_RESET
                 << " at z = " << ((1. / a) - 1.) << " (cycle " << cycle
                 << "), tau/boxsize = " << tau << endl;

            // writeSnapshots (sim, cosmo, a, dtau_old, done_hij,
            //                 snapcount, h5filename + sim.basename_snapshot,
            //                 &pcls_cdm, &pcls_b, pcls_ncdm, &PM->phi, &PM->chi,
            //                 &PM->Bi,
            //                 &PM->T00, &PM->Tij, &PM->T00_FT, &PM->Bi_FT,
            //                 &PM->Tij_FT, &PM->plan_phi,
            //                 &PM->plan_chi, &PM->plan_Bi, &PM->plan_T00,
            //                 &PM->plan_Tij
            // );

            snapcount++;
        }


        // TODO: power spectra output
        if (pkcount < sim.num_pk && 1. / a < sim.z_pk[pkcount] + 1.)
        {
            COUT << COLORTEXT_CYAN << " writing power spectra"
                 << COLORTEXT_RESET << " at z = " << ((1. / a) - 1.)
                 << " (cycle " << cycle << "), tau/boxsize = " << tau << endl;

            // writeSpectra (sim, cosmo, a, pkcount,
            //               &pcls_cdm, &pcls_b, pcls_ncdm, 
            //               &PM->phi,     &PM->chi,    &PM->Bi,
            //               &PM->T00,     &PM->Tij, 
            //               &PM->T00_FT,  &PM->Bi_FT,  &PM->Tij_FT, 
            //               &PM->plan_phi,&PM->plan_chi, 
            //               &PM->plan_Bi, &PM->plan_T00, &PM->plan_Tij
            // );

            pkcount++;
        }

#ifdef EXACT_OUTPUT_REDSHIFTS
        tmp = a;
        rungekutta4bg (tmp, cosmo, 0.5 * dtau);
        rungekutta4bg (tmp, cosmo, 0.5 * dtau);

        if (pkcount < sim.num_pk && 1. / tmp < sim.z_pk[pkcount] + 1.)
        {
            // writeSpectra (sim, cosmo, a, pkcount,
            //               &pcls_cdm, &pcls_b, pcls_ncdm, 
            //               &PM->phi, &PM->chi, &PM->Bi,
            //               &PM->T00, &PM->Tij, 
            //               &PM->T00_FT,   &PM->Bi_FT, &PM->Tij_FT, 
            //               &PM->plan_phi, &PM->plan_chi, 
            //               &PM->plan_Bi,  &PM->plan_T00, &PM->plan_Tij
            // );
        }
#endif // EXACT_OUTPUT_REDSHIFTS


        if (pkcount >= sim.num_pk && snapcount >= sim.num_snapshot)
        {
            for (i = 0; i < sim.num_lightcone; i++)
            {
                if (sim.lightcone[i].z + 1. < 1. / a)
                    i = sim.num_lightcone + 1;
            }
            if (i == sim.num_lightcone)
                break; // simulation complete
        }

        if (cycle % CYCLE_INFO_INTERVAL == 0)
        {
            COUT << " cycle " << cycle
                 << ", time integration information: max |v| = " << maxvel[0]
                 << " (cdm Courant factor = " << maxvel[0] * dtau / dx;
            
            COUT << "), time step / Hubble time = "
                 << Hconf (a, cosmo) * dtau;

            COUT << endl;
        }

        // cdm and baryon particle update
        PM->compute_forces(pcls_cdm,1.0,a);
        maxvel[0]=0;
        pcls_cdm.for_each(
            [&]
            (particle& part, const Site& /*xpart*/)
            {
               const double dtau_eff =  
                               (dtau + dtau_old) * 0.5 ;
               double v2 = 0;
               for(int i=0;i<3;++i)
               {
                   part.vel[i] += dtau_eff * part.acc[i];
                   v2 += part.vel[i]*part.vel[i];
               }
               maxvel[0]=std::max(maxvel[0],v2);
            }
            );
        maxvel[0] = std::sqrt(maxvel[0])/a;              
        
        Debugger_ptr -> flush();

        rungekutta4bg (a, cosmo,
                       0.5 * dtau); // evolve background by half a time step

        pcls_cdm.for_each(
            [&](particle& part, const Site& xpart)
            {
                std::array<Real,3> velocity =
                PM->momentum_to_velocity(
                    {part.vel[0],part.vel[1],part.vel[2]},
                    {part.pos[0],part.pos[1],part.pos[2]},
                    xpart,
                    a);
                for(int i=0;i<3;++i)
                {
                    part.pos[i] += dtau * velocity[i];
                }
            }
        );
            
        rungekutta4bg (a, cosmo,
                       0.5 * dtau); // evolve background by half a time step

        parallel.max<double> (maxvel, numspecies);

        if (sim.gr_flag == gravity_theory::GR)
        {
            for (i = 0; i < numspecies; i++)
                maxvel[i] /= sqrt (maxvel[i] * maxvel[i] + 1.0);
        }
        // done particle update

        tau += dtau;

        if (sim.wallclocklimit > 0.) // check for wallclock time limit
        {
            tmp = MPI_Wtime () - start_time;
            parallel.max (tmp);
            if (tmp > sim.wallclocklimit) // hibernate
            {
                COUT << COLORTEXT_YELLOW
                     << " reaching hibernation wallclock limit, "
                        "hibernating..."
                     << COLORTEXT_RESET << endl;
                COUT << COLORTEXT_CYAN << " writing hibernation point"
                     << COLORTEXT_RESET << " at z = " << ((1. / a) - 1.)
                     << " (cycle " << cycle << "), tau/boxsize = " << tau
                     << endl;
                if (sim.vector_flag == VECTOR_PARABOLIC 
                    && sim.gr_flag == gravity_theory::Newtonian)
                {
                    // PM->plan_Bi.execute (FFT_BACKWARD);
                }
                // hibernate (sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm,
                //            PM->phi, PM->chi, PM->Bi, a, tau, dtau, cycle);
                break;
            }
        }

        if (restartcount < sim.num_restart
            && 1. / a < sim.z_restart[restartcount] + 1.)
        {
            COUT << COLORTEXT_CYAN << " writing hibernation point"
                 << COLORTEXT_RESET << " at z = " << ((1. / a) - 1.)
                 << " (cycle " << cycle << "), tau/boxsize = " << tau << endl;
            if (sim.vector_flag == VECTOR_PARABOLIC 
                && sim.gr_flag ==gravity_theory::Newtonian)
            {
                // PM->plan_Bi.execute (FFT_BACKWARD);
            }
            // hibernate (sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, PM->phi,
            //               PM->chi, PM->Bi, a, tau, dtau, cycle, restartcount);
            restartcount++;
        }

        dtau_old = dtau;
        dtau = std::min(sim.Cf,sim.steplimit/Hconf(a,cosmo));
        cycle++;
    }while( not stop() );

    COUT << COLORTEXT_GREEN << " simulation complete." << COLORTEXT_RESET
         << endl;
    return 0;
}
