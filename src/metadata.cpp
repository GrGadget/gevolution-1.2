#include "gevolution/metadata.hpp"

namespace gevolution
{

/*
    This function builds a generic gadget header file, using the cosmology
    information, scale factor and the boxsize.
*/
gadget2_header construct_gadget_header(
    const cosmology cosmo, 
    const double a /* scale factor */, 
    const double boxsize /* in units of kpc/h */)
{
    gadget2_header hdr;
    hdr.num_files = 1;
    hdr.Omega0 = cosmo.Omega_m;
    hdr.OmegaLambda = cosmo.Omega_Lambda;
    hdr.HubbleParam = cosmo.h;
    hdr.BoxSize = boxsize;
    hdr.flag_sfr = 0;
    hdr.flag_cooling = 0;
    hdr.flag_feedback = 0;
    hdr.flag_age = 0;
    hdr.flag_metals = 0;
    for (int i = 0; i < 256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4
                        - 4 * 8 - 2 * 4 - 6 * 4;
         i++)
        hdr.fill[i] = 0;
    for (int i = 0; i < 6; i++)
    {
        hdr.npart[i] = 0;
        hdr.npartTotal[i] = 0;
        hdr.npartTotalHW[i] = 0;
        hdr.mass[i] = 0.;
    }

    hdr.time = a;
    hdr.redshift = (1. / a) - 1.;
    return hdr;
}
}
