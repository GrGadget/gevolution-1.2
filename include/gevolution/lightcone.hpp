#pragma once

#include <gevolution/basic_types.hpp>
#include <gevolution/distributed_database.hpp>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

#include <cmath>
#include <set>
#include <map>
#include <vector>
#include <utility>

#ifdef HAVE_HEALPIX
#include "chealpix.h"
#ifndef PIXBUFFER
#define PIXBUFFER 1048576
#endif
#endif

namespace gevolution {


#ifdef HAVE_HEALPIX
struct healpix_header
{
    uint32_t Nside;
    uint32_t Npix;
    uint32_t precision;
    uint32_t Ngrid;
    double direction[3];
    double distance;
    double boxsize;
    uint32_t Nside_ring;
    char fill[256 - 5 * 4 - 5 * 8]; /* fills to 256 Bytes */
};
#endif

struct lightcone_geometry
{
    // comoving position of the observer 
    // in units of the boxsize [L]
    double vertex[3];
    
    // redshift of the observation event
    double z;
    
    // conformal time of the observation event 
    // double tau;
    
    // direction of the axis of the lightcone
    double direction[3];
    
    // cosine of the angular opening of the geometric lightcone
    double opening;
    
    // comoving radius of the lightcone shell
    // distance[0] < distance[1]
    double distance[2]; 
    
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /*version*/)
    // Boost serialization of this class.
    // needed in order to broadcast this data using MPI
    {
        ar & vertex;
        ar & z;
        ar & direction;
        ar & opening;
        ar & distance;
    }
};

class parallel_lightcone
{
    ::boost::mpi::communicator com;
    
    // geometry of the full lightcone
    // FIXME: initialize this geometry
    lightcone_geometry geometry

    
    // in this way we decide where to store the information relative to each particle in the
    // distributed database.
    class destination_func
    {
        const ::boost::mpi::communicator& c;
        
        public:
        
        int operator()(const id_type I)
        {
            return I % c.size();
        }
    };
    
    distributed_database<id_type,bool,destination_func> db;
    
    // value of the conformal time (cosmic clock) when the lightcone was called
    double last_tau; 
    
    static double sqr(double x){return x*x;}
    
    static bool inside_shell(const particle& part, const lightcone_geometry& shell)
    // returns true if the particle part lies inside the lightcone shell
    {
        // comoving position difference
        // it is better to avoid confusing unit conversion here, so we decide to
        // express the coordinates of the vertex in the same units as the
        // particle position.
        const std::array<double,3> pos_diff{
            part.pos[0]-shell.vertex[0],
            part.pos[1]-shell.vertex[1],
            part.pos[2]-shell.vertex[2]
        };
        
        // comoving distance to the observer
        const double dist = std::sqrt(
              sqr(pos_diff[0])
            + sqr(pos_diff[1])
            + sqr(pos_diff[2])
            ) ; 
        
        // modulo of the direction vector
        const double dir_mod = std::sqrt( 
            sqr(shell.direction[0]) +
            sqr(shell.direction[1]) +
            sqr(shell.direction[2]) );
        
        
        // cos(alpha) = pos_diff * direction / |pos_diff| / |direction|
        // then alpha is the angle between the line of sight and the direction
        // of the lightcone.
        const double cos_angle =  
            ( pos_diff[0]*shell.direction[0] + 
              pos_diff[1]*shell.direction[1] + 
              pos_diff[2]*shell.direction[2] )
            / dir_mod / dist;
        
        // dist must be bounded by the shell inner and outer radii, while the
        // angle alpha between the line of sight and the lightcone direction
        // must not exceed the opening angle. Here cosines are compared.
        return dist >= shell.distance[0] && dist<shell.distance[1] && cos_angle >= shell.opening;
    }
    
    public:
    parallel_lightcone(
          const ::boost::mpi::communicator& in_com, 
          lightcone_geometry lc, 
          std::string fname):
        com{in_com}, geometry{lc}, db(com,destination_func{com}),
        
        // initialized to the mininum value possible of the tau (cosmic time) for an observed event
        last_tau{std::max(0.0,geometry.tau - geometry.distance[1])}
    {}
   
    // we use the conformal time (tau) as the cosmic clock, this is because tau*c is the comoving
    // distance travelled by light and we don't need to know the cosmology to compute the shell
    // boundaries of the lightcone. Had we use the scale factor a (or the redshift z) as our cosmic
    // clock here, we would have to translate that into tau integrating the Hubble factor.
    void saveLightcone(const Particles_gevolution& pcdm, 
                       const double a /* scale factor */,
                       const cosmology cosmo,
                       const std::string fname)
    {
        const double tau = particleHorizon(a,cosmo);
        
        if(tau<last_tau)
        // events with tau smaller than last_tau are not considered, they're out of the lightcone in
        // the past
            return;
        
        // at each time step we build a lightcone shell using the distance to the observer in the
        // current and the last time step.
        lightcone_geometry thin_shell{geometry};
        
        // this operation is the intersection of the total lightcone shell and
        // the thin shell of this timestep.
        thin_shell.distance[0] = std::max(thin_shell.tau - tau, geometry.distance[0]);
        thin_shell.distance[1] = std::min(thin_shell.tau - last_tau, geometric.distance[1]);
        
        if(thin_shell.distance[0]>=thin_shell.distance[1])
        // the intersection of the shell with the full lightcone is null
            return;
       
        // save tau of this time step to build the next shell
        last_tau = tau;
        
        // local particles' ids
        std::vector<id_type> my_particles;
        
        pcdm.for_each([&](const particle& p, const LATfield2::Site&)
        {
            my_particles.push_back(p.ID);
        });
        
        // query which particles are already stored in the lightcone, this is a boolean answer.
        auto in_past_lightcone = db.query(my_particles);
        
        const auto hdr = construct_gadget_header(a);
        
        // pcdm save snapshot, condition particles are inside the thin shell and not already stored
        // in the lightcone
        pcdm.saveGadget2(fname,hdr,
            [&](const particle& part, const LATfield2::Site&)
            {
                return !in_past_lightcone.at(part.ID) && inside_shell(part,thin_shell);
            });
        
        // for each local particle that lies inside the shell, set in lightcone
        pcdm.for_each([&](const particle& part, const LATfield2::Site&)
        {
            if(inside_shell(part,thin_shell))
                in_past_lightcone[part.ID] = true;
        });
        
        // update the distributed database
        db.insert(in_past_lightcone);
    }
};


} // namespace gevolution
