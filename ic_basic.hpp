//////////////////////////
// ic_basic.hpp
//////////////////////////
//
// basic initial condition generator for gevolution
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen
// Mary University of London)
//
// Last modified: November 2019
//
//////////////////////////

#ifndef IC_BASIC_HEADER
#define IC_BASIC_HEADER

#include "LATfield2.hpp"
#include "real_type.hpp"
#include "gevolution.hpp"  // solveModifiedPoissonFT
#include "parser.hpp"      // parameter
#include "prng_engine.hpp" // sitmo
#include <gsl/gsl_spline.h>

#define MAX_LINESIZE 2048

// should be larger than maximum Ngrid
#ifndef HUGE_SKIP
#define HUGE_SKIP 65536
#endif

namespace gevolution
{
using LATfield2::FFT_BACKWARD;
using LATfield2::FFT_FORWARD;
using LATfield2::Field;
using LATfield2::parallel;
using LATfield2::part_simple_dataType;
using LATfield2::part_simple_info;
using LATfield2::Particles;
using LATfield2::PlanFFT;
using LATfield2::rKSite;
using LATfield2::Site;

//////////////////////////
// displace_pcls_ic_basic
//////////////////////////
// Description:
//   displaces particles according to gradient of displacement field
//   (accepts two displacement fields for baryon treatment = hybrid)
//
// Arguments:
//   coeff             coefficient to be applied to displacement
//   lat_resolution    1 / Ngrid
//   part              particle to be displaced
//   ref_dist          distance vector (in lattice units) to reference lattice
//   point partInfo          particle metadata fields            array of
//   pointers to displacement fields sites             array of respective
//   sites nfields           number of fields passed params additional
//   parameters (ignored) outputs           array for reduction variables;
//   first entry will contain the displacement noutputs          number of
//   reduction variables
//
// Returns:
//
//////////////////////////

void displace_pcls_ic_basic (double coeff, double lat_resolution,
                             particle *part, double *ref_dist,
                             part_simple_info partInfo, Field<Real> **fields,
                             Site *sites, int nfield, double *params,
                             double *outputs, int noutputs);

//////////////////////////
// initialize_q_ic_basic
//////////////////////////
// Description:
//   initializes velocities proportional to gradient of potential
//   (accepts two potentials for baryon treatment = hybrid)
//
// Arguments:
//   coeff             coefficient to be applied to velocity
//   lat_resolution    1 / Ngrid
//   part              particle to be manipulated
//   ref_dist          distance vector (in lattice units) to reference lattice
//   point partInfo          particle metadata fields            array of
//   pointers to velocity potentials sites             array of respective
//   sites nfields           number of fields passed params additional
//   parameters (ignored) outputs           array for reduction variables
//   (ignored) noutputs          number of reduction variables
//
// Returns: square of the velocity, (q/m)^2
//
//////////////////////////

Real initialize_q_ic_basic (double coeff, double lat_resolution,
                            particle *part, double *ref_dist,
                            part_simple_info partInfo, Field<Real> **fields,
                            Site *sites, int nfield, double *params,
                            double *outputs, int noutputs);

//////////////////////////
// loadHomogeneousTemplate
//////////////////////////
// Description:
//   loads a homogeneous template from a GADGET-2 file
//
// Arguments:
//   filename   string containing the path to the template file
//   numpart    will contain the number of particles of the template
//   partdata   will contain the particle positions (memory will be allocated)
//
// Returns:
//
//////////////////////////

void loadHomogeneousTemplate (const char *filename, long &numpart,
                              float *&partdata);

//////////////////////////
// loadPowerSpectrum
//////////////////////////
// Description:
//   loads a tabulated matter power spectrum from a file
//
// Arguments:
//   filename   string containing the path to the template file
//   pkspline   will point to the gsl_spline which holds the tabulated
//              power spectrum (memory will be allocated)
//   boxsize    comoving box size (in the same units as used in the file)
//
// Returns:
//
//////////////////////////

void loadPowerSpectrum (const char *filename, gsl_spline *&pkspline,
                        const double boxsize);

//////////////////////////
// loadTransferFunctions (1)
//////////////////////////
// Description:
//   loads a set of tabulated transfer functions from a file
//
// Arguments:
//   filename   string containing the path to the template file
//   tk_delta   will point to the gsl_spline which holds the tabulated
//              transfer function for delta (memory will be allocated)
//   tk_theta   will point to the gsl_spline which holds the tabulated
//              transfer function for theta (memory will be allocated)
//   qname      string containing the name of the component (e.g. "cdm")
//   boxsize    comoving box size (in the same units as used in the file)
//   h          conversion factor between 1/Mpc and h/Mpc (theta is in units of
//   1/Mpc)
//
// Returns:
//
//////////////////////////

void loadTransferFunctions (const char *filename, gsl_spline *&tk_delta,
                            gsl_spline *&tk_theta, const char *qname,
                            const double boxsize, const double h);

//////////////////////////
// generateCICKernel
//////////////////////////
// Description:
//   generates convolution kernel for CIC projection
//
// Arguments:
//   ker        reference to allocated field that will contain the convolution
//   kernel numpcl     number of particles in the pcldata pointer pcldata raw
//   particle data from which the kernel will be constructed;
//              a standard kernel will be provided in case no particles are
//              specified
//   numtile    tiling factor used for the particle template
//
// Returns:
//
//////////////////////////

void generateCICKernel (Field<Real> &ker, const long numpcl = 0,
                        float *pcldata = nullptr, const int numtile = 1);

#ifdef FFT3D

//////////////////////////
// generateDisplacementField (generateRealization)
//////////////////////////
// Description:
//   generates particle displacement field
//
// Non-type template parameters:
//   ignorekernel  this is effectively an optimization flag defaulted to 0;
//   instantiating with 1 instead will cause
//                 the function to ignore the convolution kernel, allowing the
//                 function to be used for generating realizations
//                 (generateRealization is simply an alias for
//                 generateDisplacementField<1>)
//
// Arguments:
//   potFT         reference to allocated field that contains the convolution
//   kernel relating the potential
//                 (generating the displacement field) with the bare density
//                 perturbation; will contain the Fourier image of the
//                 potential generating the displacement field
//   coeff         gauge correction coefficient "H_conformal^2"
//   pkspline      pointer to a gsl_spline which holds a tabulated power
//   spectrum seed          initial seed for random number generator ksphere
//   flag to indicate that only a sphere in k-space should be initialized
//                 (default = 0: full k-space cube is initialized)
//   deconvolve_f  flag to indicate deconvolution function
//                 0: no deconvolution
//                 1: sinc (default)
//
// Returns:
//
//////////////////////////

#ifndef generateRealization
#define generateRealization generateDisplacementField<1>
#endif

template <int ignorekernel = 0>
void generateDisplacementField (Field<Cplx> &potFT, const Real coeff,
                                const gsl_spline *pkspline,
                                const unsigned int seed, const int ksphere = 0,
                                const int deconvolve_f = 1)
{
    const int linesize = potFT.lattice ().size (1);
    const int kmax = (linesize / 2) - 1;
    rKSite k (potFT.lattice ());
    int kx, ky, kz, i, j;
    int kymin, kymax, kzmin, kzmax;
    // long jumpy, jumpz;
    float r1, r2, k2, s;
    float *sinc;
    sitmo::prng_engine prng;
    uint64_t huge_skip = HUGE_SKIP;
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();

    sinc = (float *)malloc (linesize * sizeof (float));

    sinc[0] = 1.;
    if (deconvolve_f == 1)
    {
        for (i = 1; i < linesize; i++)
            sinc[i] = sin (M_PI * (float)i / (float)linesize) * (float)linesize
                      / (M_PI * (float)i);
    }
    else
    {
        for (i = 1; i < linesize; i++)
            sinc[i] = 1.;
    }

    k.initialize (potFT.lattice (), potFT.lattice ().siteLast ());
    kymax = k.coord (1);
    kzmax = k.coord (2);
    k.initialize (potFT.lattice (), potFT.lattice ().siteFirst ());
    kymin = k.coord (1);
    kzmin = k.coord (2);

    if (kymin < (linesize / 2) + 1 && kzmin < (linesize / 2) + 1)
    {
        prng.seed (seed);

        if (kymin == 0 && kzmin == 0)
        {
            k.setCoord (0, 0, 0);
            potFT (k) = Cplx (0., 0.);
            kx = 1;
        }
        else
        {
            kx = 0;
            prng.discard (((uint64_t)kzmin * huge_skip + (uint64_t)kymin)
                          * huge_skip);
        }

        for (kz = kzmin; kz < (linesize / 2) + 1 && kz <= kzmax; kz++)
        {
            for (ky = kymin, j = 0; ky < (linesize / 2) + 1 && ky <= kymax;
                 ky++, j++)
            {
                for (i = 0; kx < (linesize / 2) + 1; kx++)
                {
                    k.setCoord (kx, ky, kz);

                    k2 = (float)(kx * kx) + (float)(ky * ky) + (float)(kz * kz);

                    if (kx >= kmax || ky >= kmax || kz >= kmax
                        || (k2 >= kmax * kmax && ksphere > 0))
                    {
                        potFT (k) = Cplx (0., 0.);
                    }
                    else
                    {
                        s = sinc[kx] * sinc[ky] * sinc[kz];
                        k2 *= 4. * M_PI * M_PI;
                        do
                        {
                            r1 = (float)prng ()
                                 / (float)sitmo::prng_engine::max ();
                            i++;
                        } while (r1 == 0.);
                        r2 = (float)prng () / (float)sitmo::prng_engine::max ();
                        i++;

                        potFT (k)
                            = (ignorekernel
                                   ? Cplx (cos (2. * M_PI * r2),
                                           sin (2. * M_PI * r2))
                                   : Cplx (cos (2. * M_PI * r2),
                                           sin (2. * M_PI * r2))
                                         * (1. + 7.5 * coeff / k2) / potFT (k))
                              * sqrt (-2. * log (r1))
                              * gsl_spline_eval (pkspline, sqrt (k2), acc) * s;
                    }
                }
                prng.discard (huge_skip - (uint64_t)i);
                kx = 0;
            }
            prng.discard (huge_skip * (huge_skip - (uint64_t)j));
        }
    }

    if (kymax >= (linesize / 2) + 1 && kzmin < (linesize / 2) + 1)
    {
        prng.seed (seed);
        prng.discard (((huge_skip + (uint64_t)kzmin) * huge_skip
                       + (uint64_t) (linesize - kymax))
                      * huge_skip);

        for (kz = kzmin; kz < (linesize / 2) + 1 && kz <= kzmax; kz++)
        {
            for (ky = kymax, j = 0; ky >= (linesize / 2) + 1 && ky >= kymin;
                 ky--, j++)
            {
                for (kx = 0, i = 0; kx < (linesize / 2) + 1; kx++)
                {
                    k.setCoord (kx, ky, kz);

                    k2 = (float)(kx * kx)
                         + (float)((linesize - ky) * (linesize - ky))
                         + (float)(kz * kz);

                    if (kx >= kmax || (linesize - ky) >= kmax || kz >= kmax
                        || (k2 >= kmax * kmax && ksphere > 0))
                    {
                        potFT (k) = Cplx (0., 0.);
                    }
                    else
                    {
                        s = sinc[kx] * sinc[linesize - ky] * sinc[kz];
                        k2 *= 4. * M_PI * M_PI;
                        do
                        {
                            r1 = (float)prng ()
                                 / (float)sitmo::prng_engine::max ();
                            i++;
                        } while (r1 == 0.);
                        r2 = (float)prng () / (float)sitmo::prng_engine::max ();
                        i++;

                        potFT (k)
                            = (ignorekernel
                                   ? Cplx (cos (2. * M_PI * r2),
                                           sin (2. * M_PI * r2))
                                   : Cplx (cos (2. * M_PI * r2),
                                           sin (2. * M_PI * r2))
                                         * (1. + 7.5 * coeff / k2) / potFT (k))
                              * sqrt (-2. * log (r1))
                              * gsl_spline_eval (pkspline, sqrt (k2), acc) * s;
                    }
                }
                prng.discard (huge_skip - (uint64_t)i);
            }
            prng.discard (huge_skip * (huge_skip - (uint64_t)j));
        }
    }

    if (kymin < (linesize / 2) + 1 && kzmax >= (linesize / 2) + 1)
    {
        prng.seed (seed);
        prng.discard (
            ((huge_skip + huge_skip + (uint64_t) (linesize - kzmax)) * huge_skip
             + (uint64_t)kymin)
            * huge_skip);

        for (kz = kzmax; kz >= (linesize / 2) + 1 && kz >= kzmin; kz--)
        {
            for (ky = kymin, j = 0; ky < (linesize / 2) + 1 && ky <= kymax;
                 ky++, j++)
            {
                for (kx = 1, i = 0; kx < (linesize / 2) + 1; kx++)
                {
                    k.setCoord (kx, ky, kz);

                    k2 = (float)(kx * kx) + (float)(ky * ky)
                         + (float)((linesize - kz) * (linesize - kz));

                    if (kx >= kmax || ky >= kmax || (linesize - kz) >= kmax
                        || (k2 >= kmax * kmax && ksphere > 0))
                    {
                        potFT (k) = Cplx (0., 0.);
                    }
                    else
                    {
                        s = sinc[kx] * sinc[ky] * sinc[linesize - kz];
                        k2 *= 4. * M_PI * M_PI;
                        do
                        {
                            r1 = (float)prng ()
                                 / (float)sitmo::prng_engine::max ();
                            i++;
                        } while (r1 == 0.);
                        r2 = (float)prng () / (float)sitmo::prng_engine::max ();
                        i++;

                        potFT (k)
                            = (ignorekernel
                                   ? Cplx (cos (2. * M_PI * r2),
                                           sin (2. * M_PI * r2))
                                   : Cplx (cos (2. * M_PI * r2),
                                           sin (2. * M_PI * r2))
                                         * (1. + 7.5 * coeff / k2) / potFT (k))
                              * sqrt (-2. * log (r1))
                              * gsl_spline_eval (pkspline, sqrt (k2), acc) * s;
                    }
                }
                prng.discard (huge_skip - (uint64_t)i);
            }
            prng.discard (huge_skip * (huge_skip - (uint64_t)j));
        }

        prng.seed (seed);
        prng.discard (
            ((uint64_t) (linesize - kzmax) * huge_skip + (uint64_t)kymin)
            * huge_skip);
        kx = 0;

        for (kz = kzmax; kz >= (linesize / 2) + 1 && kz >= kzmin; kz--)
        {
            for (ky = kymin, j = 0; ky < (linesize / 2) + 1 && ky <= kymax;
                 ky++, j++)
            {
                k.setCoord (kx, ky, kz);

                k2 = (float)(ky * ky)
                     + (float)((linesize - kz) * (linesize - kz));
                i = 0;

                if (ky >= kmax || (linesize - kz) >= kmax
                    || (k2 >= kmax * kmax && ksphere > 0))
                {
                    potFT (k) = Cplx (0., 0.);
                }
                else
                {
                    s = sinc[ky] * sinc[linesize - kz];
                    k2 *= 4. * M_PI * M_PI;
                    do
                    {
                        r1 = (float)prng () / (float)sitmo::prng_engine::max ();
                        i++;
                    } while (r1 == 0.);
                    r2 = (float)prng () / (float)sitmo::prng_engine::max ();
                    i++;

                    potFT (k)
                        = (ignorekernel
                               ? Cplx (cos (2. * M_PI * r2),
                                       -sin (2. * M_PI * r2))
                               : Cplx (cos (2. * M_PI * r2),
                                       -sin (2. * M_PI * r2))
                                     * (1. + 7.5 * coeff / k2) / potFT (k))
                          * sqrt (-2. * log (r1))
                          * gsl_spline_eval (pkspline, sqrt (k2), acc) * s;
                }

                prng.discard (huge_skip - (uint64_t)i);
            }
            prng.discard (huge_skip * (huge_skip - (uint64_t)j));
        }
    }

    if (kymax >= (linesize / 2) + 1 && kzmax >= (linesize / 2) + 1)
    {
        prng.seed (seed);
        prng.discard (
            ((huge_skip + huge_skip + huge_skip + (uint64_t) (linesize - kzmax))
                 * huge_skip
             + (uint64_t) (linesize - kymax))
            * huge_skip);

        for (kz = kzmax; kz >= (linesize / 2) + 1 && kz >= kzmin; kz--)
        {
            for (ky = kymax, j = 0; ky >= (linesize / 2) + 1 && ky >= kymin;
                 ky--, j++)
            {
                for (kx = 1, i = 0; kx < (linesize / 2) + 1; kx++)
                {
                    k.setCoord (kx, ky, kz);

                    k2 = (float)(kx * kx)
                         + (float)((linesize - ky) * (linesize - ky))
                         + (float)((linesize - kz) * (linesize - kz));

                    if (kx >= kmax || (linesize - ky) >= kmax
                        || (linesize - kz) >= kmax
                        || (k2 >= kmax * kmax && ksphere > 0))
                    {
                        potFT (k) = Cplx (0., 0.);
                    }
                    else
                    {
                        s = sinc[kx] * sinc[linesize - ky]
                            * sinc[linesize - kz];
                        k2 *= 4. * M_PI * M_PI;
                        do
                        {
                            r1 = (float)prng ()
                                 / (float)sitmo::prng_engine::max ();
                            i++;
                        } while (r1 == 0.);
                        r2 = (float)prng () / (float)sitmo::prng_engine::max ();
                        i++;

                        potFT (k)
                            = (ignorekernel
                                   ? Cplx (cos (2. * M_PI * r2),
                                           sin (2. * M_PI * r2))
                                   : Cplx (cos (2. * M_PI * r2),
                                           sin (2. * M_PI * r2))
                                         * (1. + 7.5 * coeff / k2) / potFT (k))
                              * sqrt (-2. * log (r1))
                              * gsl_spline_eval (pkspline, sqrt (k2), acc) * s;
                    }
                }
                prng.discard (huge_skip - (uint64_t)i);
            }
            prng.discard (huge_skip * (huge_skip - (uint64_t)j));
        }

        prng.seed (seed);
        prng.discard (
            ((huge_skip + huge_skip + (uint64_t) (linesize - kzmax)) * huge_skip
             + (uint64_t) (linesize - kymax))
            * huge_skip);
        kx = 0;

        for (kz = kzmax; kz >= (linesize / 2) + 1 && kz >= kzmin; kz--)
        {
            for (ky = kymax, j = 0; ky >= (linesize / 2) + 1 && ky >= kymin;
                 ky--, j++)
            {
                k.setCoord (kx, ky, kz);

                k2 = (float)((linesize - ky) * (linesize - ky))
                     + (float)((linesize - kz) * (linesize - kz));
                i = 0;

                if ((linesize - ky) >= kmax || (linesize - kz) >= kmax
                    || (k2 >= kmax * kmax && ksphere > 0))
                {
                    potFT (k) = Cplx (0., 0.);
                }
                else
                {
                    s = sinc[linesize - ky] * sinc[linesize - kz];
                    k2 *= 4. * M_PI * M_PI;
                    do
                    {
                        r1 = (float)prng () / (float)sitmo::prng_engine::max ();
                        i++;
                    } while (r1 == 0.);
                    r2 = (float)prng () / (float)sitmo::prng_engine::max ();
                    i++;

                    potFT (k)
                        = (ignorekernel
                               ? Cplx (cos (2. * M_PI * r2),
                                       -sin (2. * M_PI * r2))
                               : Cplx (cos (2. * M_PI * r2),
                                       -sin (2. * M_PI * r2))
                                     * (1. + 7.5 * coeff / k2) / potFT (k))
                          * sqrt (-2. * log (r1))
                          * gsl_spline_eval (pkspline, sqrt (k2), acc) * s;
                }

                prng.discard (huge_skip - (uint64_t)i);
            }
            prng.discard (huge_skip * (huge_skip - (uint64_t)j));
        }
    }

    gsl_interp_accel_free (acc);
    free (sinc);
}
#endif

//////////////////////////
// initializeParticlePositions
//////////////////////////
// Description:
//   initializes particle positions using a homogeneous template
//
// Arguments:
//   numpart    number of particles of the template
//   partdata   particle positions in the template
//   numtile    tiling factor for homogeneous template - total particle number
//   will be
//              numpart * numtile^3
//   pcls       reference to (empty) particle object which will contain the new
//   particle ensemble
//
// Returns:
//
//////////////////////////

void initializeParticlePositions (
    const long numpart, const float *partdata, const int numtile,
    Particles_gevolution &pcls);

//////////////////////////
// applyMomentumDistribution
//////////////////////////
// Description:
//   adds a random momentum vector drawn from a Fermi-Dirac distribution to
//   each particle. The current implementation uses a simple rejection-sampling
//   method based on the ziggurat algorithm [G. Marsaglia and W.W. Tsang,
//   J. Stat. Softw. 5 (2000) 1]. The "ziggurat" is hard-coded and was
//   precomputed for the ultrarelativistic limit of a Fermi-Dirac distribution
//   with zero chemical potential. A method to construct the ziggurat
//   on-the-fly for more general distribution functions could be implemented in
//   the future.
//
// Arguments:
//   pcls   pointer to particle handler
//   seed   seed for random number generator
//   T_m    dimensionless parameter in the distribution function
//          in most cases the ratio of temperature and fundamental mass
//   delta  Field containing a local dT/T (optional)
//
// Returns: sum of momenta over all particles (for reduction)
//
//////////////////////////

double applyMomentumDistribution (
    Particles_gevolution *pcls,
    unsigned int seed, float T_m = 0., Field<Real> *delta = nullptr);

#ifdef FFT3D

//////////////////////////
// generateIC_basic
//////////////////////////
// Description:
//   basic initial condition generator
//
// Arguments:
//   sim            simulation metadata structure
//   ic             settings for IC generation
//   cosmo          cosmological parameter structure
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
//   params         pointer to array of precision settings for CLASS (can be
//   NULL) numparam       number of precision settings for CLASS (can be 0)
//
// Returns:
//
//////////////////////////

void generateIC_basic (
    metadata &sim, icsettings &ic, const cosmology cosmo,
    Particles_gevolution *pcls_cdm,
    Particles_gevolution *pcls_b,
    Particles_gevolution *pcls_ncdm,
    double *maxvel, Field<Real> *phi, Field<Real> *chi, Field<Real> *Bi,
    Field<Real> *source, Field<Real> *Sij, Field<Cplx> *scalarFT,
    Field<Cplx> *BiFT, Field<Cplx> *SijFT, PlanFFT<Cplx> *plan_phi,
    PlanFFT<Cplx> *plan_chi, PlanFFT<Cplx> *plan_Bi, PlanFFT<Cplx> *plan_source,
    PlanFFT<Cplx> *plan_Sij, parameter *params, int &numparam);

#endif
}
#endif
