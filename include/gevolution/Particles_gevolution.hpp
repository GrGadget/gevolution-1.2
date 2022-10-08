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

#include <boost/mpi/collectives.hpp>


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
    using base_particle = LATfield2::part_simple;
    using base_particle::ID;
    using base_particle::pos;
    using base_particle::vel;
    
    // inherited from LATfield2::part_simple we have
    // pos[]; // that's position (x)
    // vel[]; // that's velocity (dx/dt)
    // id;
    Real mass;
    std::array<Real,3> momentum{0,0,0}; // that's momentum (p)
    std::array<Real,3> force{0,0,0};    // that's force (dp/dt)
    
    // Metric components at the particle position
    Real Phi{0};
    std::array<Real,3> B{0,0,0};
    
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /*version*/)
    {
        // inherited from LATfield2::part_simple
        ar & ID;
        // notive that LATfield2::part_simple::pos is not an std::array
        // bad design
        ar & pos[0] ;
        ar & pos[1] ;
        ar & pos[2] ; 
        
        ar & vel[0] ;
        ar & vel[1] ;
        ar & vel[2] ;
        
        // augmented structure
        ar & mass;
        ar & momentum;
        ar & force;
        
        ar & Phi;
        ar & B;
    }
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
    
    template<typename selector_type>
    void saveGadget2 (
        const std::string filename, 
        const gadget2_header hdr,
        selector_type select_function)const
    {
        // count the particles that meet the criteria
        long long n_part = 0 ;
        this->for_each(
            [&](const particles& part, const LATfield2::Site& x)
            {
                n_part += select_function(part) ? 1 : 0;
            }
        );
         
        
        // open filename as w
        // set hdr.numpart
        // write the hdr
        
        // take care of parallel
        // write the pos
        // write the momentum
        // write the ID
    }

    
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
    
    const MPI_Comm& communicator()const
    {
        // TODO get a cartesian communicator
        return ::LATfield2::parallel.my_comm;
    }
};

}
#endif
