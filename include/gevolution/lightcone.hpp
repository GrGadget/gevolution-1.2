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

struct lightcone_data
{
    // comoving position of the observer 
    // in units of the boxsize [L]
    double vertex[3];
    
    // redshift of the observation event
    double redshift_observer;
    
    // direction of the axis of the lightcone
    double direction[3];
    
    // cosine of the angular opening of the geometric lightcone
    double cos_opening;
    
    // redshift interval of observation
    double redshift[2]; 
};

struct lightcone_geometry
{
    // comoving position of the observer 
    // in units of the boxsize [L]
    double vertex[3];
    
    // redshift of the observation event
    // double z;
    
    // conformal time of the observation event 
    double tau;
    
    // direction of the axis of the lightcone
    double direction[3];
    
    // cosine of the angular opening of the geometric lightcone
    double cos_opening;
    
    // comoving radius of the lightcone shell
    // r_low < r_high
    // an event that occured at conformal time t_e is observed at comoving
    // distance r_e = t_o - t_e, where t_o is the conformal time of the
    // observation. Therefore we could have simply encoded the observation
    // events in conformal time as 
    // t_early = t_o - r_low 
    // t_latest = t_o - r_high
    // then only events happening in conformal time t, with 
    // t_early <= t <= t_latest
    // will be observed within this lightcone.
    double r_low, r_high; 
    
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /*version*/)
    // Boost serialization of this class.
    // needed in order to broadcast this data using MPI
    {
        ar & vertex;
        ar & tau;
        ar & direction;
        ar & cos_opening;
        ar & r_low & r_high;
    }
};

class parallel_lightcone
{
    typedef gadget_serialization::id_type id_type;
    ::boost::mpi::communicator com;
    
    // geometry of the full lightcone
    // FIXME: initialize this geometry
    lightcone_geometry geometry;

    
    // in this way we decide where to store the information relative to each particle in the
    // distributed database.
    class destination_func
    {
        const ::boost::mpi::communicator& c;
        
        public:
        
        destination_func(::boost::mpi::communicator const & com): 
            c{com} 
        {}
        
        int operator()(const id_type I) const
        {
            return I % c.size();
        }
    };
    
    distributed_database<id_type,bool,destination_func> db;
    
    // value of the conformal time (cosmic clock) when the lightcone was called
    double last_tau{0}; 
    
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
        return dist >= shell.r_low && dist<shell.r_high 
            && cos_angle >= shell.cos_opening;
    }
    
    public:
    parallel_lightcone(
          const ::boost::mpi::communicator& in_com, 
          lightcone_geometry lc):
        com{in_com}, geometry{lc}, db(com,destination_func{com}),
        
        // initialized to the mininum value possible of the tau (cosmic time) for an observed event
        last_tau{std::max(0.0,geometry.tau - geometry.r_high)}
    {}
   
    // we use the conformal time (tau) as the cosmic clock, this is because tau*c is the comoving
    // distance travelled by light and we don't need to know the cosmology to compute the shell
    // boundaries of the lightcone. Had we use the scale factor a (or the redshift z) as our cosmic
    // clock here, we would have to translate that into tau integrating the Hubble factor.
    void saveLightcone(const Particles_gevolution& pcdm, 
                       const double a /* scale factor */,
                       const cosmology cosmo,
                       const double boxsize /* boxsize in kpc/h */, 
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
        thin_shell.r_low = std::max(thin_shell.tau - tau, geometry.r_low);
        thin_shell.r_high = std::min(thin_shell.tau - last_tau, geometry.r_high);
        
        if(thin_shell.r_low>=thin_shell.r_high)
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
        std::vector< std::pair<id_type,bool> > query_result = db.query(my_particles);
       
        auto cmp_function =
            [](
                const std::pair<id_type,bool>& a, 
                const std::pair<id_type,bool>& b)
            {
                return a.first < b.first;
            };
       
        // we are going to use this vector to query particle by id, so lets
        // order it by id and the search will be O(log n) with binary search
        std::sort(query_result.begin(),query_result.end(),cmp_function);
        
        auto position_in_query = [&](id_type const id)
        {
            // search in the query array
            auto it =
            std::lower_bound(query_result.begin(),query_result.end(),std::make_pair(id,false),cmp_function);
            if (it == query_result.end())
                throw std::runtime_error("parallel_lightcone::saveLightcone id not found in query_result");
            
            return std::distance(query_result.begin(),it);
        };
        
        auto in_past_lightcone = [&](id_type const id)
        {
            return query_result.at( position_in_query(id) ).second;
        };
        
        
        const auto hdr = construct_gadget_header(cosmo,a,boxsize);
        
        // pcdm save snapshot, condition particles are inside the thin shell and not already stored
        // in the lightcone
        pcdm.saveGadget2(fname,hdr,
            [&](const particle& part, const LATfield2::Site&)
            {
                return !in_past_lightcone(part.ID) && inside_shell(part,thin_shell);
            });
        
        // for each local particle that lies inside the shell, set in lightcone
        pcdm.for_each([&](const particle& part, const LATfield2::Site&)
        {
            if(inside_shell(part,thin_shell))
                query_result.at ( position_in_query(part.ID)).second = true;
        });
        
        // update the distributed database
        db.insert(query_result);
    }
};


} // namespace gevolution
