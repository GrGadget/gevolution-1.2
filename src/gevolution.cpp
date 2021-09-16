//////////////////////////
//
// Geneva algorithms for evolution of metric perturbations
// and relativistic free-streaming particles (gevolution)
//
// 1. Suite of Fourier-based methods for the computation of the
//    relativistic scalar (Phi, Phi-Psi) and vector modes [see J. Adamek,
//    R. Durrer, and M. Kunz, Class. Quant. Grav. 31, 234006 (2014)]
//
// 2. Collection of "update position" and "update velocity/momentum" methods
//    [see J. Adamek, D. Daverio, R. Durrer, and M. Kunz, JCAP 1607, 053
//    (2016)]
//
// 3. Collection of projection methods for the construction of the
//    stress-energy-tensor
//
// 4. Fourier-space projection methods for the computation of the
//    curl and divergence of the velocity field
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen
// Mary University of London)
//
// Last modified: April 2019
//
//////////////////////////

#include "gevolution/debugger.hpp"
#include "gevolution/gevolution.hpp"
#include <cstdlib>
#include <iostream>

namespace gevolution
{
using LATfield2::Field;
using LATfield2::Particles;
using LATfield2::rKSite;
using LATfield2::Site;

#ifdef FFT3D
//////////////////////////
// projectFTscalar
//////////////////////////
// Description:
//   projection of the Fourier image of a tensor field on the trace-free
//   longitudinal (scalar) component
//
// Arguments:
//   SijFT      reference to the Fourier image of the input tensor field
//   chiFT      reference to allocated field which will contain the Fourier
//              image of the trace-free longitudinal (scalar) component
//
// Returns:
//
//////////////////////////

void projectFTscalar (Field<Cplx> &SijFT, Field<Cplx> &chiFT, const int add)
{
    const int linesize = chiFT.lattice ().size (1);
    int i;
    Real *gridk2;
    Cplx *kshift;
    rKSite k (chiFT.lattice ());

    gridk2 = (Real *)malloc (linesize * sizeof (Real));
    kshift = (Cplx *)malloc (linesize * sizeof (Cplx));

    for (i = 0; i < linesize; i++)
    {
        gridk2[i] = 2. * (Real)linesize * sin (M_PI * (Real)i / (Real)linesize);
        kshift[i] = gridk2[i]
                    * Cplx (cos (M_PI * (Real)i / (Real)linesize),
                            -sin (M_PI * (Real)i / (Real)linesize));
        gridk2[i] *= gridk2[i];
    }

    k.first ();
    if (k.coord (0) == 0 && k.coord (1) == 0 && k.coord (2) == 0)
    {
        chiFT (k) = Cplx (0., 0.);
        k.next ();
    }

    if (add)
    {
        for (; k.test (); k.next ())
        {
            chiFT (k) += ((gridk2[k.coord (1)] + gridk2[k.coord (2)]
                           - 2. * gridk2[k.coord (0)])
                              * SijFT (k, 0, 0)
                          + (gridk2[k.coord (0)] + gridk2[k.coord (2)]
                             - 2. * gridk2[k.coord (1)])
                                * SijFT (k, 1, 1)
                          + (gridk2[k.coord (0)] + gridk2[k.coord (1)]
                             - 2. * gridk2[k.coord (2)])
                                * SijFT (k, 2, 2)
                          - 6. * kshift[k.coord (0)] * kshift[k.coord (1)]
                                * SijFT (k, 0, 1)
                          - 6. * kshift[k.coord (0)] * kshift[k.coord (2)]
                                * SijFT (k, 0, 2)
                          - 6. * kshift[k.coord (1)] * kshift[k.coord (2)]
                                * SijFT (k, 1, 2))
                         / (2.
                            * (gridk2[k.coord (0)] + gridk2[k.coord (1)]
                               + gridk2[k.coord (2)])
                            * (gridk2[k.coord (0)] + gridk2[k.coord (1)]
                               + gridk2[k.coord (2)])
                            * linesize);
        }
    }
    else
    {
        for (; k.test (); k.next ())
        {
            chiFT (k) = ((gridk2[k.coord (1)] + gridk2[k.coord (2)]
                          - 2. * gridk2[k.coord (0)])
                             * SijFT (k, 0, 0)
                         + (gridk2[k.coord (0)] + gridk2[k.coord (2)]
                            - 2. * gridk2[k.coord (1)])
                               * SijFT (k, 1, 1)
                         + (gridk2[k.coord (0)] + gridk2[k.coord (1)]
                            - 2. * gridk2[k.coord (2)])
                               * SijFT (k, 2, 2)
                         - 6. * kshift[k.coord (0)] * kshift[k.coord (1)]
                               * SijFT (k, 0, 1)
                         - 6. * kshift[k.coord (0)] * kshift[k.coord (2)]
                               * SijFT (k, 0, 2)
                         - 6. * kshift[k.coord (1)] * kshift[k.coord (2)]
                               * SijFT (k, 1, 2))
                        / (2.
                           * (gridk2[k.coord (0)] + gridk2[k.coord (1)]
                              + gridk2[k.coord (2)])
                           * (gridk2[k.coord (0)] + gridk2[k.coord (1)]
                              + gridk2[k.coord (2)])
                           * linesize);
        }
    }

    free (gridk2);
    free (kshift);
}

//////////////////////////
// evolveFTvector
//////////////////////////
// Description:
//   projects the Fourier image of a tensor field on the spin-1 component
//   used as a source for the evolution of the vector perturbation
//
// Arguments:
//   SijFT      reference to the Fourier image of the input tensor field
//   BiFT       reference to the Fourier image of the vector perturbation
//   a2dtau     conformal time step times scale factor squared (a^2 * dtau)
//
// Returns:
//
//////////////////////////

void evolveFTvector (Field<Cplx> &SijFT, Field<Cplx> &BiFT, const Real a2dtau)
{
    const int linesize = BiFT.lattice ().size (1);
    int i;
    Real *gridk2;
    Cplx *kshift;
    rKSite k (BiFT.lattice ());
    Real k4;

    gridk2 = (Real *)malloc (linesize * sizeof (Real));
    kshift = (Cplx *)malloc (linesize * sizeof (Cplx));

    for (i = 0; i < linesize; i++)
    {
        gridk2[i] = 2. * (Real)linesize * sin (M_PI * (Real)i / (Real)linesize);
        kshift[i] = gridk2[i]
                    * Cplx (cos (M_PI * (Real)i / (Real)linesize),
                            -sin (M_PI * (Real)i / (Real)linesize));
        gridk2[i] *= gridk2[i];
    }

    k.first ();
    if (k.coord (0) == 0 && k.coord (1) == 0 && k.coord (2) == 0)
    {
        BiFT (k, 0) = Cplx (0., 0.);
        BiFT (k, 1) = Cplx (0., 0.);
        BiFT (k, 2) = Cplx (0., 0.);
        k.next ();
    }

    for (; k.test (); k.next ())
    {
        k4 = gridk2[k.coord (0)] + gridk2[k.coord (1)] + gridk2[k.coord (2)];
        k4 *= k4;

        BiFT (k, 0) += Cplx (0., -2. * a2dtau / k4)
                       * (kshift[k.coord (0)].conj ()
                              * ((gridk2[k.coord (1)] + gridk2[k.coord (2)])
                                     * SijFT (k, 0, 0)
                                 - gridk2[k.coord (1)] * SijFT (k, 1, 1)
                                 - gridk2[k.coord (2)] * SijFT (k, 2, 2)
                                 - 2. * kshift[k.coord (1)]
                                       * kshift[k.coord (2)] * SijFT (k, 1, 2))
                          + (gridk2[k.coord (1)] + gridk2[k.coord (2)]
                             - gridk2[k.coord (0)])
                                * (kshift[k.coord (1)] * SijFT (k, 0, 1)
                                   + kshift[k.coord (2)] * SijFT (k, 0, 2)));
        BiFT (k, 1) += Cplx (0., -2. * a2dtau / k4)
                       * (kshift[k.coord (1)].conj ()
                              * ((gridk2[k.coord (0)] + gridk2[k.coord (2)])
                                     * SijFT (k, 1, 1)
                                 - gridk2[k.coord (0)] * SijFT (k, 0, 0)
                                 - gridk2[k.coord (2)] * SijFT (k, 2, 2)
                                 - 2. * kshift[k.coord (0)]
                                       * kshift[k.coord (2)] * SijFT (k, 0, 2))
                          + (gridk2[k.coord (0)] + gridk2[k.coord (2)]
                             - gridk2[k.coord (1)])
                                * (kshift[k.coord (0)] * SijFT (k, 0, 1)
                                   + kshift[k.coord (2)] * SijFT (k, 1, 2)));
        BiFT (k, 2) += Cplx (0., -2. * a2dtau / k4)
                       * (kshift[k.coord (2)].conj ()
                              * ((gridk2[k.coord (0)] + gridk2[k.coord (1)])
                                     * SijFT (k, 2, 2)
                                 - gridk2[k.coord (0)] * SijFT (k, 0, 0)
                                 - gridk2[k.coord (1)] * SijFT (k, 1, 1)
                                 - 2. * kshift[k.coord (0)]
                                       * kshift[k.coord (1)] * SijFT (k, 0, 1))
                          + (gridk2[k.coord (0)] + gridk2[k.coord (1)]
                             - gridk2[k.coord (2)])
                                * (kshift[k.coord (0)] * SijFT (k, 0, 2)
                                   + kshift[k.coord (1)] * SijFT (k, 1, 2)));
    }

    free (gridk2);
    free (kshift);
}

//////////////////////////
// projectFTvector
//////////////////////////
// Description:
//   projects the Fourier image of a vector field on the transverse component
//   and solves the constraint equation for the vector perturbation
//
// Arguments:
//   SiFT       reference to the Fourier image of the input vector field
//   BiFT       reference to the Fourier image of the vector perturbation (can
//   be identical to input) coeff      rescaling coefficient (default 1) modif
//   modification k^2 -> k^2 + modif (default 0)
//
// Returns:
//
//////////////////////////

void projectFTvector (Field<Cplx> &SiFT, Field<Cplx> &BiFT, const Real coeff,
                      const Real modif)
{
    const int linesize = BiFT.lattice ().size (1);
    int i;
    Real *gridk2;
    Cplx *kshift;
    rKSite k (BiFT.lattice ());
    Real k2;
    Cplx tmp (0., 0.);

    gridk2 = (Real *)malloc (linesize * sizeof (Real));
    kshift = (Cplx *)malloc (linesize * sizeof (Cplx));

    for (i = 0; i < linesize; i++)
    {
        gridk2[i] = 2. * (Real)linesize * sin (M_PI * (Real)i / (Real)linesize);
        kshift[i] = gridk2[i]
                    * Cplx (cos (M_PI * (Real)i / (Real)linesize),
                            -sin (M_PI * (Real)i / (Real)linesize));
        gridk2[i] *= gridk2[i];
    }

    k.first ();
    if (k.coord (0) == 0 && k.coord (1) == 0 && k.coord (2) == 0)
    {
        BiFT (k, 0) = Cplx (0., 0.);
        BiFT (k, 1) = Cplx (0., 0.);
        BiFT (k, 2) = Cplx (0., 0.);
        k.next ();
    }

    for (; k.test (); k.next ())
    {
        k2 = gridk2[k.coord (0)] + gridk2[k.coord (1)] + gridk2[k.coord (2)];

        tmp = (kshift[k.coord (0)] * SiFT (k, 0)
               + kshift[k.coord (1)] * SiFT (k, 1)
               + kshift[k.coord (2)] * SiFT (k, 2))
              / k2;

        BiFT (k, 0) = (SiFT (k, 0) - kshift[k.coord (0)].conj () * tmp) * 4.
                      * coeff / (k2 + modif);
        BiFT (k, 1) = (SiFT (k, 1) - kshift[k.coord (1)].conj () * tmp) * 4.
                      * coeff / (k2 + modif);
        BiFT (k, 2) = (SiFT (k, 2) - kshift[k.coord (2)].conj () * tmp) * 4.
                      * coeff / (k2 + modif);
    }

    free (gridk2);
    free (kshift);
}

//////////////////////////
// projectFTtensor
//////////////////////////
// Description:
//   projection of the Fourier image of a tensor field on the transverse
//   trace-free tensor component
//
// Arguments:
//   SijFT      reference to the Fourier image of the input tensor field
//   hijFT      reference to allocated field which will contain the Fourier
//              image of the transverse trace-free tensor component
//
// Returns:
//
//////////////////////////

void projectFTtensor (Field<Cplx> &SijFT, Field<Cplx> &hijFT)
{
    const int linesize = hijFT.lattice ().size (1);
    int i;
    Real *gridk2;
    Cplx *kshift;
    rKSite k (hijFT.lattice ());
    Cplx SxxFT, SxyFT, SxzFT, SyyFT, SyzFT, SzzFT;
    Real k2, k6;

    gridk2 = (Real *)malloc (linesize * sizeof (Real));
    kshift = (Cplx *)malloc (linesize * sizeof (Cplx));

    for (i = 0; i < linesize; i++)
    {
        gridk2[i] = 2. * (Real)linesize * sin (M_PI * (Real)i / (Real)linesize);
        kshift[i] = gridk2[i]
                    * Cplx (cos (M_PI * (Real)i / (Real)linesize),
                            -sin (M_PI * (Real)i / (Real)linesize));
        gridk2[i] *= gridk2[i];
    }

    k.first ();
    if (k.coord (0) == 0 && k.coord (1) == 0 && k.coord (2) == 0)
    {
        for (i = 0; i < hijFT.components (); i++)
            hijFT (k, i) = Cplx (0., 0.);

        k.next ();
    }

    for (; k.test (); k.next ())
    {
        SxxFT = SijFT (k, 0, 0);
        SxyFT = SijFT (k, 0, 1);
        SxzFT = SijFT (k, 0, 2);
        SyyFT = SijFT (k, 1, 1);
        SyzFT = SijFT (k, 1, 2);
        SzzFT = SijFT (k, 2, 2);

        k2 = gridk2[k.coord (0)] + gridk2[k.coord (1)] + gridk2[k.coord (2)];
        k6 = k2 * k2 * k2 * linesize;

        hijFT (k, 0, 0)
            = ((gridk2[k.coord (0)] - k2)
                   * ((gridk2[k.coord (0)] - k2) * SxxFT
                      + 2. * kshift[k.coord (0)]
                            * (kshift[k.coord (1)] * SxyFT
                               + kshift[k.coord (2)] * SxzFT))
               + ((gridk2[k.coord (0)] + k2) * (gridk2[k.coord (1)] + k2)
                  - 2. * k2 * k2)
                     * SyyFT
               + ((gridk2[k.coord (0)] + k2) * (gridk2[k.coord (2)] + k2)
                  - 2. * k2 * k2)
                     * SzzFT
               + 2. * (gridk2[k.coord (0)] + k2) * kshift[k.coord (1)]
                     * kshift[k.coord (2)] * SyzFT)
              / k6;

        hijFT (k, 0, 1)
            = (2. * (gridk2[k.coord (0)] - k2) * (gridk2[k.coord (1)] - k2)
                   * SxyFT
               + (gridk2[k.coord (2)] + k2) * kshift[k.coord (0)].conj ()
                     * kshift[k.coord (1)].conj () * SzzFT
               + (gridk2[k.coord (0)] - k2) * kshift[k.coord (1)].conj ()
                     * (kshift[k.coord (0)].conj () * SxxFT
                        + 2. * kshift[k.coord (2)] * SxzFT)
               + (gridk2[k.coord (1)] - k2) * kshift[k.coord (0)].conj ()
                     * (kshift[k.coord (1)].conj () * SyyFT
                        + 2. * kshift[k.coord (2)] * SyzFT))
              / k6;

        hijFT (k, 0, 2)
            = (2. * (gridk2[k.coord (0)] - k2) * (gridk2[k.coord (2)] - k2)
                   * SxzFT
               + (gridk2[k.coord (1)] + k2) * kshift[k.coord (0)].conj ()
                     * kshift[k.coord (2)].conj () * SyyFT
               + (gridk2[k.coord (0)] - k2) * kshift[k.coord (2)].conj ()
                     * (kshift[k.coord (0)].conj () * SxxFT
                        + 2. * kshift[k.coord (1)] * SxyFT)
               + (gridk2[k.coord (2)] - k2) * kshift[k.coord (0)].conj ()
                     * (kshift[k.coord (2)].conj () * SzzFT
                        + 2. * kshift[k.coord (1)] * SyzFT))
              / k6;

        hijFT (k, 1, 1)
            = ((gridk2[k.coord (1)] - k2)
                   * ((gridk2[k.coord (1)] - k2) * SyyFT
                      + 2. * kshift[k.coord (1)]
                            * (kshift[k.coord (0)] * SxyFT
                               + kshift[k.coord (2)] * SyzFT))
               + ((gridk2[k.coord (1)] + k2) * (gridk2[k.coord (0)] + k2)
                  - 2. * k2 * k2)
                     * SxxFT
               + ((gridk2[k.coord (1)] + k2) * (gridk2[k.coord (2)] + k2)
                  - 2. * k2 * k2)
                     * SzzFT
               + 2. * (gridk2[k.coord (1)] + k2) * kshift[k.coord (0)]
                     * kshift[k.coord (2)] * SxzFT)
              / k6;

        hijFT (k, 1, 2)
            = (2. * (gridk2[k.coord (1)] - k2) * (gridk2[k.coord (2)] - k2)
                   * SyzFT
               + (gridk2[k.coord (0)] + k2) * kshift[k.coord (1)].conj ()
                     * kshift[k.coord (2)].conj () * SxxFT
               + (gridk2[k.coord (1)] - k2) * kshift[k.coord (2)].conj ()
                     * (kshift[k.coord (1)].conj () * SyyFT
                        + 2. * kshift[k.coord (0)] * SxyFT)
               + (gridk2[k.coord (2)] - k2) * kshift[k.coord (1)].conj ()
                     * (kshift[k.coord (2)].conj () * SzzFT
                        + 2. * kshift[k.coord (0)] * SxzFT))
              / k6;

        hijFT (k, 2, 2)
            = ((gridk2[k.coord (2)] - k2)
                   * ((gridk2[k.coord (2)] - k2) * SzzFT
                      + 2. * kshift[k.coord (2)]
                            * (kshift[k.coord (0)] * SxzFT
                               + kshift[k.coord (1)] * SyzFT))
               + ((gridk2[k.coord (2)] + k2) * (gridk2[k.coord (0)] + k2)
                  - 2. * k2 * k2)
                     * SxxFT
               + ((gridk2[k.coord (2)] + k2) * (gridk2[k.coord (1)] + k2)
                  - 2. * k2 * k2)
                     * SyyFT
               + 2. * (gridk2[k.coord (2)] + k2) * kshift[k.coord (0)]
                     * kshift[k.coord (1)] * SxyFT)
              / k6;
    }

    free (gridk2);
    free (kshift);
}

#endif

//////////////////////////
// update_q
//////////////////////////
// Description:
//   Update momentum method (arbitrary momentum)
//   Note that vel[3] in the particle structure is used to store q[3] in units
//   of the particle mass, such that as q^2 << m^2 a^2 the meaning of vel[3]
//   is ~ v*a.
//
// Arguments:
//   dtau       time step
//   dx         lattice unit
//   part       pointer to particle structure
//   ref_dist   distance vector to reference point
//   partInfo   global particle properties (unused)
//   fields     array of pointers to fields appearing in geodesic equation
//              fields[0] = phi
//              fields[1] = chi
//              fields[2] = Bi
//   sites      array of sites on the respective lattices
//   nfield     number of fields
//   params     array of additional parameters
//              params[0] = a
//              params[1] = scaling coefficient for Bi
//   outputs    array of reduction variables
//   noutputs   number of reduction variables
//
// Returns: squared velocity of particle after update
//
//////////////////////////

Real update_q (double dtau, double dx, particle *part, double *ref_dist,
               particle_info partInfo, Field<Real> **fields, Site *sites,
               int nfield, double *params, double *outputs, int noutputs)
{
#define phi (*fields[0])
#define chi (*fields[1])
#define Bi (*fields[2])
#define xphi (sites[0])
#define xchi (sites[1])
#define xB (sites[2])

    Real gradphi[3] = { 0, 0, 0 };
    Real pgradB[3] = { 0, 0, 0 };
    Real v2 = (*part).vel[0] * (*part).vel[0] + (*part).vel[1] * (*part).vel[1]
              + (*part).vel[2] * (*part).vel[2];
    Real e2 = v2 + params[0] * params[0];

    gradphi[0] = (1. - ref_dist[1]) * (1. - ref_dist[2])
                 * (phi (xphi + 0) - phi (xphi));
    gradphi[1] = (1. - ref_dist[0]) * (1. - ref_dist[2])
                 * (phi (xphi + 1) - phi (xphi));
    gradphi[2] = (1. - ref_dist[0]) * (1. - ref_dist[1])
                 * (phi (xphi + 2) - phi (xphi));
    gradphi[0] += ref_dist[1] * (1. - ref_dist[2])
                  * (phi (xphi + 1 + 0) - phi (xphi + 1));
    gradphi[1] += ref_dist[0] * (1. - ref_dist[2])
                  * (phi (xphi + 1 + 0) - phi (xphi + 0));
    gradphi[2] += ref_dist[0] * (1. - ref_dist[1])
                  * (phi (xphi + 2 + 0) - phi (xphi + 0));
    gradphi[0] += (1. - ref_dist[1]) * ref_dist[2]
                  * (phi (xphi + 2 + 0) - phi (xphi + 2));
    gradphi[1] += (1. - ref_dist[0]) * ref_dist[2]
                  * (phi (xphi + 2 + 1) - phi (xphi + 2));
    gradphi[2] += (1. - ref_dist[0]) * ref_dist[1]
                  * (phi (xphi + 2 + 1) - phi (xphi + 1));
    gradphi[0] += ref_dist[1] * ref_dist[2]
                  * (phi (xphi + 2 + 1 + 0) - phi (xphi + 2 + 1));
    gradphi[1] += ref_dist[0] * ref_dist[2]
                  * (phi (xphi + 2 + 1 + 0) - phi (xphi + 2 + 0));
    gradphi[2] += ref_dist[0] * ref_dist[1]
                  * (phi (xphi + 2 + 1 + 0) - phi (xphi + 1 + 0));

    gradphi[0] *= (v2 + e2) / e2;
    gradphi[1] *= (v2 + e2) / e2;
    gradphi[2] *= (v2 + e2) / e2;

    if (nfield >= 2 && fields[1] != NULL)
    {
        gradphi[0] -= (1. - ref_dist[1]) * (1. - ref_dist[2])
                      * (chi (xchi + 0) - chi (xchi));
        gradphi[1] -= (1. - ref_dist[0]) * (1. - ref_dist[2])
                      * (chi (xchi + 1) - chi (xchi));
        gradphi[2] -= (1. - ref_dist[0]) * (1. - ref_dist[1])
                      * (chi (xchi + 2) - chi (xchi));
        gradphi[0] -= ref_dist[1] * (1. - ref_dist[2])
                      * (chi (xchi + 1 + 0) - chi (xchi + 1));
        gradphi[1] -= ref_dist[0] * (1. - ref_dist[2])
                      * (chi (xchi + 1 + 0) - chi (xchi + 0));
        gradphi[2] -= ref_dist[0] * (1. - ref_dist[1])
                      * (chi (xchi + 2 + 0) - chi (xchi + 0));
        gradphi[0] -= (1. - ref_dist[1]) * ref_dist[2]
                      * (chi (xchi + 2 + 0) - chi (xchi + 2));
        gradphi[1] -= (1. - ref_dist[0]) * ref_dist[2]
                      * (chi (xchi + 2 + 1) - chi (xchi + 2));
        gradphi[2] -= (1. - ref_dist[0]) * ref_dist[1]
                      * (chi (xchi + 2 + 1) - chi (xchi + 1));
        gradphi[0] -= ref_dist[1] * ref_dist[2]
                      * (chi (xchi + 2 + 1 + 0) - chi (xchi + 2 + 1));
        gradphi[1] -= ref_dist[0] * ref_dist[2]
                      * (chi (xchi + 2 + 1 + 0) - chi (xchi + 2 + 0));
        gradphi[2] -= ref_dist[0] * ref_dist[1]
                      * (chi (xchi + 2 + 1 + 0) - chi (xchi + 1 + 0));
    }

    e2 = sqrt (e2);

    if (nfield >= 3 && fields[2] != NULL)
    {
        pgradB[0] = ((1. - ref_dist[2]) * (Bi (xB + 0, 1) - Bi (xB, 1))
                     + ref_dist[2] * (Bi (xB + 2 + 0, 1) - Bi (xB + 2, 1)))
                    * (*part).vel[1];
        pgradB[0] += ((1. - ref_dist[1]) * (Bi (xB + 0, 2) - Bi (xB, 2))
                      + ref_dist[1] * (Bi (xB + 1 + 0, 2) - Bi (xB + 1, 2)))
                     * (*part).vel[2];
        pgradB[0] += (1. - ref_dist[1]) * (1. - ref_dist[2])
                     * ((ref_dist[0] - 1.) * Bi (xB - 0, 0)
                        + (1. - 2. * ref_dist[0]) * Bi (xB, 0)
                        + ref_dist[0] * Bi (xB + 0, 0))
                     * (*part).vel[0];
        pgradB[0] += ref_dist[1] * (1. - ref_dist[2])
                     * ((ref_dist[0] - 1.) * Bi (xB + 1 - 0, 0)
                        + (1. - 2. * ref_dist[0]) * Bi (xB + 1, 0)
                        + ref_dist[0] * Bi (xB + 1 + 0, 0))
                     * (*part).vel[0];
        pgradB[0] += (1. - ref_dist[1]) * ref_dist[2]
                     * ((ref_dist[0] - 1.) * Bi (xB + 2 - 0, 0)
                        + (1. - 2. * ref_dist[0]) * Bi (xB + 2, 0)
                        + ref_dist[0] * Bi (xB + 2 + 0, 0))
                     * (*part).vel[0];
        pgradB[0] += ref_dist[1] * ref_dist[2]
                     * ((ref_dist[0] - 1.) * Bi (xB + 2 + 1 - 0, 0)
                        + (1. - 2. * ref_dist[0]) * Bi (xB + 2 + 1, 0)
                        + ref_dist[0] * Bi (xB + 2 + 1 + 0, 0))
                     * (*part).vel[0];

        pgradB[1] = ((1. - ref_dist[0]) * (Bi (xB + 1, 2) - Bi (xB, 2))
                     + ref_dist[0] * (Bi (xB + 1 + 0, 2) - Bi (xB + 0, 2)))
                    * (*part).vel[2];
        pgradB[1] += ((1. - ref_dist[2]) * (Bi (xB + 1, 0) - Bi (xB, 0))
                      + ref_dist[2] * (Bi (xB + 1 + 2, 0) - Bi (xB + 2, 0)))
                     * (*part).vel[0];
        pgradB[1] += (1. - ref_dist[0]) * (1. - ref_dist[2])
                     * ((ref_dist[1] - 1.) * Bi (xB - 1, 1)
                        + (1. - 2. * ref_dist[1]) * Bi (xB, 1)
                        + ref_dist[1] * Bi (xB + 1, 1))
                     * (*part).vel[1];
        pgradB[1] += ref_dist[0] * (1. - ref_dist[2])
                     * ((ref_dist[1] - 1.) * Bi (xB + 0 - 1, 1)
                        + (1. - 2. * ref_dist[1]) * Bi (xB + 0, 1)
                        + ref_dist[1] * Bi (xB + 0 + 1, 1))
                     * (*part).vel[1];
        pgradB[1] += (1. - ref_dist[0]) * ref_dist[2]
                     * ((ref_dist[1] - 1.) * Bi (xB + 2 - 1, 1)
                        + (1. - 2. * ref_dist[1]) * Bi (xB + 2, 1)
                        + ref_dist[1] * Bi (xB + 2 + 1, 1))
                     * (*part).vel[1];
        pgradB[1] += ref_dist[0] * ref_dist[2]
                     * ((ref_dist[1] - 1.) * Bi (xB + 2 + 0 - 1, 1)
                        + (1. - 2. * ref_dist[1]) * Bi (xB + 2 + 0, 1)
                        + ref_dist[1] * Bi (xB + 2 + 0 + 1, 1))
                     * (*part).vel[1];

        pgradB[2] = ((1. - ref_dist[1]) * (Bi (xB + 2, 0) - Bi (xB, 0))
                     + ref_dist[1] * (Bi (xB + 2 + 1, 0) - Bi (xB + 1, 0)))
                    * (*part).vel[0];
        pgradB[2] += ((1. - ref_dist[0]) * (Bi (xB + 2, 1) - Bi (xB, 1))
                      + ref_dist[0] * (Bi (xB + 2 + 0, 1) - Bi (xB + 0, 1)))
                     * (*part).vel[1];
        pgradB[2] += (1. - ref_dist[0]) * (1. - ref_dist[1])
                     * ((ref_dist[2] - 1.) * Bi (xB - 2, 2)
                        + (1. - 2. * ref_dist[2]) * Bi (xB, 2)
                        + ref_dist[2] * Bi (xB + 2, 2))
                     * (*part).vel[2];
        pgradB[2] += ref_dist[0] * (1. - ref_dist[1])
                     * ((ref_dist[2] - 1.) * Bi (xB + 0 - 2, 2)
                        + (1. - 2. * ref_dist[2]) * Bi (xB + 0, 2)
                        + ref_dist[2] * Bi (xB + 0 + 2, 2))
                     * (*part).vel[2];
        pgradB[2] += (1. - ref_dist[0]) * ref_dist[1]
                     * ((ref_dist[2] - 1.) * Bi (xB + 1 - 2, 2)
                        + (1. - 2. * ref_dist[2]) * Bi (xB + 1, 2)
                        + ref_dist[2] * Bi (xB + 2 + 1, 2))
                     * (*part).vel[2];
        pgradB[2] += ref_dist[0] * ref_dist[1]
                     * ((ref_dist[2] - 1.) * Bi (xB + 1 + 0 - 2, 2)
                        + (1. - 2. * ref_dist[2]) * Bi (xB + 1 + 0, 2)
                        + ref_dist[2] * Bi (xB + 1 + 0 + 2, 2))
                     * (*part).vel[2];

        gradphi[0] += pgradB[0] / params[1] / e2;
        gradphi[1] += pgradB[1] / params[1] / e2;
        gradphi[2] += pgradB[2] / params[1] / e2;
    }

    v2 = 0.;
    for (int i = 0; i < 3; i++)
    {
        (*part).vel[i] -= dtau * e2 * gradphi[i] / dx;
        v2 += (*part).vel[i] * (*part).vel[i];
    }

    return v2 / params[0] / params[0];

#undef phi
#undef chi
#undef Bi
#undef xphi
#undef xchi
#undef xB
}

//////////////////////////
// update_q_Newton
//////////////////////////
// Description:
//   Update momentum method (Newtonian version)
//   Note that vel[3] in the particle structure is used to store q[3] in units
//   of the particle mass, such that the meaning of vel[3] is v*a.
//
// Arguments:
//   dtau       time step
//   dx         lattice unit
//   part       pointer to particle structure
//   ref_dist   distance vector to reference point
//   partInfo   global particle properties (unused)
//   fields     array of pointers to fields appearing in geodesic equation
//              fields[0] = psi
//              fields[1] = chi
//   sites      array of sites on the respective lattices
//   nfield     number of fields (should be 1)
//   params     array of additional parameters
//              params[0] = a
//   outputs    array of reduction variables
//   noutputs   number of reduction variables
//
// Returns: squared velocity of particle after update
//
//////////////////////////
Real update_q_Newton ( 
                      particle& part,
                      const Field<Real>& psi, 
                      const Site& xpart,
                      double dtau,
                      double dx,
                      double a)
{
    std::array<double,3> ref_dist;
    for(int l=0;l<3;++l)
        ref_dist[l] = part.pos[l]/dx - xpart.coord(l);
    
    std::array<Real,3> gradpsi{ 0, 0, 0 };

    gradpsi[0] = (1. - ref_dist[1]) * (1. - ref_dist[2])
                 * (psi (xpart + 0) - psi (xpart));
    gradpsi[1] = (1. - ref_dist[0]) * (1. - ref_dist[2])
                 * (psi (xpart + 1) - psi (xpart));
    gradpsi[2] = (1. - ref_dist[0]) * (1. - ref_dist[1])
                 * (psi (xpart + 2) - psi (xpart));
    gradpsi[0] += ref_dist[1] * (1. - ref_dist[2])
                  * (psi (xpart + 1 + 0) - psi (xpart + 1));
    gradpsi[1] += ref_dist[0] * (1. - ref_dist[2])
                  * (psi (xpart + 1 + 0) - psi (xpart + 0));
    gradpsi[2] += ref_dist[0] * (1. - ref_dist[1])
                  * (psi (xpart + 2 + 0) - psi (xpart + 0));
    gradpsi[0] += (1. - ref_dist[1]) * ref_dist[2]
                  * (psi (xpart + 2 + 0) - psi (xpart + 2));
    gradpsi[1] += (1. - ref_dist[0]) * ref_dist[2]
                  * (psi (xpart + 2 + 1) - psi (xpart + 2));
    gradpsi[2] += (1. - ref_dist[0]) * ref_dist[1]
                  * (psi (xpart + 2 + 1) - psi (xpart + 1));
    gradpsi[0] += ref_dist[1] * ref_dist[2]
                  * (psi (xpart + 2 + 1 + 0) - psi (xpart + 2 + 1));
    gradpsi[1] += ref_dist[0] * ref_dist[2]
                  * (psi (xpart + 2 + 1 + 0) - psi (xpart + 2 + 0));
    gradpsi[2] += ref_dist[0] * ref_dist[1]
                  * (psi (xpart + 2 + 1 + 0) - psi (xpart + 1 + 0));

    Real v2 = 0.;
    std::array<double,3> acc;
    for (int i = 0; i < 3; i++)
    {
        acc[i] =  (-1)*a * gradpsi[i] / dx;
        part.vel[i] += dtau * acc[i];
        v2 += part.vel[i] * part.vel[i];
    }
    return v2 / a / a;
}

//////////////////////////
// update_pos
//////////////////////////
// Description:
//   Update position method (arbitrary momentum)
//   Note that vel[3] in the particle structure is used to store q[3] in units
//   of the particle mass, such that as q^2 << m^2 a^2 the meaning of vel[3]
//   is ~ v*a.
//
// Arguments:
//   dtau       time step
//   dx         lattice unit
//   part       pointer to particle structure
//   ref_dist   distance vector to reference point
//   partInfo   global particle properties (unused)
//   fields     array of pointers to fields appearing in geodesic equation
//              fields[0] = phi
//              fields[1] = chi
//              fields[2] = Bi
//   sites      array of sites on the respective lattices
//   nfield     number of fields
//   params     array of additional parameters
//              params[0] = a
//              params[1] = scaling coefficient for Bi
//   outputs    array of reduction variables
//   noutputs   number of reduction variables
//
// Returns:
//
//////////////////////////

void update_pos (double dtau, double dx, particle *part, double *ref_dist,
                 particle_info partInfo, Field<Real> **fields, Site *sites,
                 int nfield, double *params, double *outputs, int noutputs)
{
    Real v[3];
    Real v2 = (*part).vel[0] * (*part).vel[0] + (*part).vel[1] * (*part).vel[1]
              + (*part).vel[2] * (*part).vel[2];
    Real e2 = v2 + params[0] * params[0];
    Real phi = 0;
    Real chi = 0;

    if (nfield >= 1)
    {
        phi = (*fields[0]) (sites[0]) * (1. - ref_dist[0]) * (1. - ref_dist[1])
              * (1. - ref_dist[2]);
        phi += (*fields[0]) (sites[0] + 0) * ref_dist[0] * (1. - ref_dist[1])
               * (1. - ref_dist[2]);
        phi += (*fields[0]) (sites[0] + 1) * (1. - ref_dist[0]) * ref_dist[1]
               * (1. - ref_dist[2]);
        phi += (*fields[0]) (sites[0] + 0 + 1) * ref_dist[0] * ref_dist[1]
               * (1. - ref_dist[2]);
        phi += (*fields[0]) (sites[0] + 2) * (1. - ref_dist[0])
               * (1. - ref_dist[1]) * ref_dist[2];
        phi += (*fields[0]) (sites[0] + 0 + 2) * ref_dist[0]
               * (1. - ref_dist[1]) * ref_dist[2];
        phi += (*fields[0]) (sites[0] + 1 + 2) * (1. - ref_dist[0])
               * ref_dist[1] * ref_dist[2];
        phi += (*fields[0]) (sites[0] + 0 + 1 + 2) * ref_dist[0] * ref_dist[1]
               * ref_dist[2];
    }

    if (nfield >= 2)
    {
        chi = (*fields[1]) (sites[1]) * (1. - ref_dist[0]) * (1. - ref_dist[1])
              * (1. - ref_dist[2]);
        chi += (*fields[1]) (sites[1] + 0) * ref_dist[0] * (1. - ref_dist[1])
               * (1. - ref_dist[2]);
        chi += (*fields[1]) (sites[1] + 1) * (1. - ref_dist[0]) * ref_dist[1]
               * (1. - ref_dist[2]);
        chi += (*fields[1]) (sites[1] + 0 + 1) * ref_dist[0] * ref_dist[1]
               * (1. - ref_dist[2]);
        chi += (*fields[1]) (sites[1] + 2) * (1. - ref_dist[0])
               * (1. - ref_dist[1]) * ref_dist[2];
        chi += (*fields[1]) (sites[1] + 0 + 2) * ref_dist[0]
               * (1. - ref_dist[1]) * ref_dist[2];
        chi += (*fields[1]) (sites[1] + 1 + 2) * (1. - ref_dist[0])
               * ref_dist[1] * ref_dist[2];
        chi += (*fields[1]) (sites[1] + 0 + 1 + 2) * ref_dist[0] * ref_dist[1]
               * ref_dist[2];
    }

    v2 = (1. + (3. - v2 / e2) * phi - chi) / sqrt (e2);

    v[0] = (*part).vel[0] * v2;
    v[1] = (*part).vel[1] * v2;
    v[2] = (*part).vel[2] * v2;

    if (nfield >= 3)
    {
        Real b[3];

        b[0] = (*fields[2]) (sites[2], 0) * (1. - ref_dist[1])
               * (1. - ref_dist[2]);
        b[1] = (*fields[2]) (sites[2], 1) * (1. - ref_dist[0])
               * (1. - ref_dist[2]);
        b[2] = (*fields[2]) (sites[2], 2) * (1. - ref_dist[0])
               * (1. - ref_dist[1]);
        b[1] += (*fields[2]) (sites[2] + 0, 1) * ref_dist[0]
                * (1. - ref_dist[2]);
        b[2] += (*fields[2]) (sites[2] + 0, 2) * ref_dist[0]
                * (1. - ref_dist[1]);
        b[0] += (*fields[2]) (sites[2] + 1, 0) * ref_dist[1]
                * (1. - ref_dist[2]);
        b[2] += (*fields[2]) (sites[2] + 1, 2) * (1. - ref_dist[0])
                * ref_dist[1];
        b[0] += (*fields[2]) (sites[2] + 2, 0) * (1. - ref_dist[1])
                * ref_dist[2];
        b[1] += (*fields[2]) (sites[2] + 2, 1) * (1. - ref_dist[0])
                * ref_dist[2];
        b[1] += (*fields[2]) (sites[2] + 2 + 0, 1) * ref_dist[0] * ref_dist[2];
        b[0] += (*fields[2]) (sites[2] + 2 + 1, 0) * ref_dist[1] * ref_dist[2];
        b[2] += (*fields[2]) (sites[2] + 1 + 0, 2) * ref_dist[0] * ref_dist[1];

        for (int l = 0; l < 3; l++)
            (*part).pos[l] += dtau * (v[l] + b[l] / params[1]);
    }
    else
    {
        for (int l = 0; l < 3; l++)
            (*part).pos[l] += dtau * v[l];
    }
}

//////////////////////////
// update_pos_Newton
//////////////////////////
// Description:
//   Update position method (Newtonian version)
//   Note that vel[3] in the particle structure is used to store q[3] in units
//   of the particle mass, such that the meaning of vel[3] is v*a.
//
// Arguments:
//   dtau       time step
//   dx         lattice unit (unused)
//   part       pointer to particle structure
//   ref_dist   distance vector to reference point (unused)
//   partInfo   global particle properties (unused)
//   fields     array of pointers to fields appearing in geodesic equation
//   (unused) sites      array of sites on the respective lattices (unused)
//   nfield     number of fields (unused)
//   params     array of additional parameters
//              params[0] = a
//   outputs    array of reduction variables (unused)
//   noutputs   number of reduction variables (unused)
//
// Returns:
//
//////////////////////////

void update_pos_Newton (double dtau, double dx, particle *part,
                        double *ref_dist, particle_info partInfo,
                        Field<Real> **fields, Site *sites, int nfield,
                        double *params, double *outputs, int noutputs)
{
    for (int l = 0; l < 3; l++)
        (*part).pos[l] += dtau * (*part).vel[l] / params[0];
}

//////////////////////////
// projectFTtheta
//////////////////////////
// Description:
//   Compute the diverge of the velocity in Fourier space
//
// Arguments:
//   thFT       reference to the Fourier image of the divergence of the
//   velocity field viFT       reference to the Fourier image of the velocity
//   field
//
// Returns:
//
//////////////////////////

void projectFTtheta (Field<Cplx> &thFT, Field<Cplx> &viFT)
{
    const int linesize = thFT.lattice ().size (1);
    int i;
    Real *gridk;
    rKSite k (thFT.lattice ());
    Cplx tmp (0., 0.);

    gridk = (Real *)malloc (linesize * sizeof (Real));

    for (i = 0; i < linesize; i++)
        gridk[i] = (Real)linesize * sin (M_PI * 2.0 * (Real)i / (Real)linesize);

    for (k.first (); k.test (); k.next ())
        thFT (k) = Cplx (0., 1.)
                   * (gridk[k.coord (0)] * viFT (k, 0)
                      + gridk[k.coord (1)] * viFT (k, 1)
                      + gridk[k.coord (2)] * viFT (k, 2));

    free (gridk);
}

//////////////////////////
// projectFTomega
//////////////////////////
// Description:
//   Compute the curl part of the velocity field in Fourier space
//
// Arguments:
//   viFT      reference to the input Fourier image of the velocity field
//             the divergence part will be projected out
//
// Returns:
//
//////////////////////////

void projectFTomega (Field<Cplx> &viFT)
{
    const int linesize = viFT.lattice ().size (1);
    int i;
    Real *gridk2;
    // Cplx *kshift;
    Real *gridk;
    rKSite k (viFT.lattice ());
    Cplx tmp (0., 0.);
    Cplx vr[3];

    gridk2 = (Real *)malloc (linesize * sizeof (Real));
    gridk = (Real *)malloc (linesize * sizeof (Real));

    for (i = 0; i < linesize; i++)
    {
        gridk[i] = (Real)linesize * sin (M_PI * 2.0 * (Real)i / (Real)linesize);
        gridk2[i] = gridk[i] * gridk[i];
    }

    k.first ();
    if (k.coord (0) == 0 && k.coord (1) == 0 && k.coord (2) == 0)
    {
        viFT (k, 0) = Cplx (0., 0.);
        viFT (k, 1) = Cplx (0., 0.);
        viFT (k, 2) = Cplx (0., 0.);
        k.next ();
    }

    for (; k.test (); k.next ())
    {
        if ((k.coord (0) == 0 || k.coord (0) == linesize / 2)
            && (k.coord (1) == 0 || k.coord (1) == linesize / 2)
            && (k.coord (2) == 0 || k.coord (2) == linesize / 2))
        {
            viFT (k, 0) = Cplx (0., 0.);
            viFT (k, 1) = Cplx (0., 0.);
            viFT (k, 2) = Cplx (0., 0.);
        }
        else
        {
            tmp = (gridk[k.coord (0)] * viFT (k, 0)
                   + gridk[k.coord (1)] * viFT (k, 1)
                   + gridk[k.coord (2)] * viFT (k, 2))
                  / (gridk2[k.coord (0)] + gridk2[k.coord (1)]
                     + gridk2[k.coord (2)]);

            vr[0] = (viFT (k, 0) - gridk[k.coord (0)] * tmp);
            vr[1] = (viFT (k, 1) - gridk[k.coord (1)] * tmp);
            vr[2] = (viFT (k, 2) - gridk[k.coord (2)] * tmp);

            viFT (k, 0)
                = Cplx (0., 1.)
                  * (gridk[k.coord (1)] * vr[2] - gridk[k.coord (2)] * vr[1]);
            viFT (k, 1)
                = Cplx (0., 1.)
                  * (gridk[k.coord (2)] * vr[0] - gridk[k.coord (0)] * vr[2]);
            viFT (k, 2)
                = Cplx (0., 1.)
                  * (gridk[k.coord (0)] * vr[1] - gridk[k.coord (1)] * vr[0]);
        }
    }

    free (gridk2);
}
}
