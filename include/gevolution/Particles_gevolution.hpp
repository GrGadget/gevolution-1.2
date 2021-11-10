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

#include "gevolution/config.h"
#include "LATfield2.hpp"
#include "gevolution/real_type.hpp"
#include "gevolution/metadata.hpp"
#include <cstdlib>
#include <iostream>
#include <mpi.h>
#include <set>
#include <string>


namespace gevolution
{

using LATfield2::Field;
using LATfield2::parallel;
using LATfield2::Site;

/*
    The indivual particle data structure in gevolution is 'particle',
    defined below.
*/
// TODO: template particle variables
struct particle : LATfield2::part_simple
{
    // inherited from LATfield2::part_simple we have
    // pos[]; // that's position (x)
    // vel[]; // that's velocity (dx/dt)
    // id;
    Real mass;
    std::array<Real,3> momentum{0,0,0}; // that's momentum (p)
    std::array<Real,3> force{0,0,0};    // that's force (dp/dt)
};

typedef LATfield2::part_simple_info particle_info;
typedef LATfield2::part_simple_dataType particle_dataType;

class Particles_gevolution : 
    public 
    LATfield2::Particles<
        particle, 
        particle_info, 
        particle_dataType>
{
  public:
    using value_type = particle; // helper attribute used in containers for metaprogramming

    void saveGadget2 (std::string filename, gadget2_header &hdr,
                      const int tracer_factor = 1) const ;
    void saveGadget2 (std::string filename, gadget2_header &hdr,
                      lightcone_geometry &lightcone, double dist, double dtau,
                      double dtau_old, double dadtau,
                      double vertex[MAX_INTERSECTS][3], const int vertexcount,
                      std::set<long> &IDbacklog, std::set<long> &IDprelog,
                      Field<Real> *phi, const int tracer_factor = 1);
    void loadGadget2 (std::string filename, gadget2_header &hdr);
    
    void update_mass()
    {
        this->for_each(
            [&](particle& part,const LATfield2::Site&)
            {
                part.mass = this->parts_info()->mass;
            });
    }
};

}
#endif
