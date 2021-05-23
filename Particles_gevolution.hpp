//////////////////////////
// Particles_gevolution.hpp
//////////////////////////
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen
// Mary University of London)
//
// Last modified: March 2020
//
//////////////////////////

#ifndef PARTICLES_GEVOLUTION_HEADER
#define PARTICLES_GEVOLUTION_HEADER

#ifndef PCLBUFFER
#define PCLBUFFER 1048576
#endif
#include "LATfield2.hpp"
#include "metadata.hpp"
#include <cstdlib>
#include <iostream>
#include <mpi.h>
#include <set>
#include <string>

namespace gevolution
{

using LATfield2::Field;
using LATfield2::parallel;
using LATfield2::Real;
using LATfield2::Site;

class Particles_gevolution : 
    public 
    LATfield2::Particles<
        LATfield2::part_simple, 
        LATfield2::part_simple_info, 
        LATfield2::part_simple_dataType>
{
  public:
    void saveGadget2 (std::string filename, gadget2_header &hdr,
                      const int tracer_factor = 1, double dtau_pos = 0.,
                      double dtau_vel = 0., Field<Real> *phi = NULL);
    void saveGadget2 (std::string filename, gadget2_header &hdr,
                      lightcone_geometry &lightcone, double dist, double dtau,
                      double dtau_old, double dadtau,
                      double vertex[MAX_INTERSECTS][3], const int vertexcount,
                      std::set<long> &IDbacklog, std::set<long> &IDprelog,
                      Field<Real> *phi, const int tracer_factor = 1);
    void loadGadget2 (std::string filename, gadget2_header &hdr);
};


}
#endif
