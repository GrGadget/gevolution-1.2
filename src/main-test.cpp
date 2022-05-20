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
bool stop(mpi::communicator &com_world)
{
    bool ret = false;
    if(com_world.rank()==0)
    {
        fs::path p{"stop"};
        if(fs::exists(p))
        {
            ret = true;
            fs::remove(p);
        }
    }
    mpi::all_reduce(com_world,ret,ret,std::logical_or());
    return ret;
}

void show_constant(double v, std::string name)
{
    COUT << "Constant " << name << " = " << v << "\n";
}

auto string_fill(std::string s, int n , char c = ' ')
{
    std::string prefix;
    n -= s.size();
    if(n>0)
    {
        prefix = std::string(n,c);
    }
    return prefix + s;
}


int main (int argc, char **argv)
{
    mpi::environment env;
    mpi::communicator com_world;
    
    int n = 0, m = 0;

    int cycle = 0, snapcount = 0, pkcount = 0, restartcount = 0,
              usedparams, numparam = 0, numspecies;
    // int done_hij;
    // int numsteps_ncdm[MAX_PCL_SPECIES - 2];
    // long numpts3d;
    int box[3];
    double dtau, dtau_old, dx, tau, a, tmp, start_time;
    double maxvel[MAX_PCL_SPECIES];
    // FILE *outfile;
    char filename[2 * PARAM_MAX_LENGTH + 24];
    string h5filename;
    char *settingsfile = NULL;

    parameter *params = NULL;
    metadata sim;
    cosmology cosmo;
    icsettings ic;

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

    numspecies = 1 + sim.baryon_flag + cosmo.num_ncdm;
    parallel.max<double> (maxvel, numspecies);

    if (sim.gr_flag == gravity_theory::GR)
    {
        for (int i = 0; i < numspecies; i++)
            maxvel[i] /= sqrt (maxvel[i] * maxvel[i] + 1.0);
    }


        COUT << COLORTEXT_GREEN << " initialization complete."
             << COLORTEXT_RESET << endl
             << endl;

    
    std::unique_ptr< particle_mesh<Cplx,Particles_gevolution> > PM;
    
    if(sim.gr_flag==gravity_theory::GR)
    {
        PM.reset(
            new relativistic_pm<Cplx,Particles_gevolution>(sim.numpts,com_world)
        );
    }else
    {
        PM.reset(
            new newtonian_pm<Cplx,Particles_gevolution>(sim.numpts,com_world)
        );
    }
    
    pcls_cdm.update_mass(); // fix the mass legacy problem
    
    // background file initialization
    fs::path BackgroundPath{sim.output_path};
    BackgroundPath /= std::string(sim.basename_generic) + "_background.dat";
    constexpr int tabwidth = 16;
    
    if(com_world.rank()==0)
    // select main process
    {
        std::ofstream f(BackgroundPath);
        f << "# background statistics\n"
          << "#" 
          << std::setw(tabwidth-1) << "cycle"
          << std::setw(tabwidth) << "tau"
          << std::setw(tabwidth) << "a"
          << std::setw(tabwidth) << "conformal H/H0"
          << std::setw(tabwidth) << "sum(phi)"
          << std::setw(tabwidth) << "<T00>"
          << "\n";
    }
    
    do // main loop
    {
        COUT << "Starting cycle: " << cycle << '\n';        
        
        for(int i=0;i<10;++i)
        {
            // PM step 1. construction of the energy momentum tensor
            PM->clear_sources();
            PM->sample(pcls_cdm,a);
            // PM step 2. compute the potentials
            PM->compute_potential(
                cosmo.fourpiG, 
                a,
                Hconf(a,cosmo),
                // dtau_old,
                cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm (a, cosmo));
        }
        COUT << "fourpiG = " << cosmo.fourpiG << '\n';
        COUT << "a       = " << a << '\n';
        
        // TODO: power spectra output
        {
            COUT << COLORTEXT_CYAN << " writing power spectra"
                 << COLORTEXT_RESET << " at z = " << ((1. / a) - 1.)
                 << " (cycle " << cycle << "), tau/boxsize = " << tau << endl;
            
            PM->save_power_spectrum("power"+std::to_string(pkcount));
        }

    }while( false );

    COUT << COLORTEXT_GREEN << " simulation complete." << COLORTEXT_RESET
         << endl;
    return 0;
}
