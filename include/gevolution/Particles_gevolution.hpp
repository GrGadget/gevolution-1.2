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
#include "gevolution/basic_types.hpp"
#include "gevolution/metadata.hpp"
#include "gevolution/final_action.hpp"
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
    
    /* Like the previous saveGadget2 function. But this one serializes only the
     * subset of particles for which selector_type()==true. If the functor is
     * always true, then this function should produce the same serialization as
     * the first saveGadget2. */
    template<typename selector_type>
    void saveGadget2 (
        const std::string filename, 
        gadget2_header hdr,
        selector_type select_function) const
    // precondition: hdr must already contain cosmological data and number of particles=0
    {
    
        // count the particles that meet the criteria
        int64_t npart = 0 ;
        this->for_each(
            [&](const particle& part, const LATfield2::Site& x)
            {
                npart += (select_function(part,x)) ? 1 : 0;
            }
        );
        
        auto const & cart_com = cartesian_communicator();
        const int64_t total_npart 
            = ::boost::mpi::all_reduce(cart_com,npart,std::plus<decltype(npart)>());
        
        // the sum of all npart up to rank()-1
        const int64_t partial_count 
            = ::boost::mpi::scan(cart_com,npart,std::plus<decltype(npart)>()) - npart;
        
        if(total_npart==0)
        // if no particles are selected should we output and empy snapshot?
            return;
        
        // set the number of particles in the snapshot header
        hdr.npart[1] = hdr.npartTotal[1] = static_cast<uint32_t>(total_npart % (1LL << 32));
        hdr.npartTotalHW[1] = static_cast<uint32_t>(total_npart / (1LL<<32));
    
        MPI_File outfile;
        MPI_File_open (cart_com, filename.c_str(),
                       MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,
                       &outfile);
        
        // RAII at work here. On destruction this object will release the file
        // resource.
        auto close_mpi_file = finally([&](){
            MPI_File_close(&outfile); 
        });
        
        const MPI_Offset filesize = 
        /* 3 positions + 3 velocities + 1 id per particle*/
        total_npart*(3*sizeof(gadget_serialization::position_type)
        +3*sizeof(gadget_serialization::velocity_type)+sizeof(gadget_serialization::id_type)) 
        /* header */
         +sizeof(hdr)
        /* the file consist of 4 blocks (HDR,POS,VEL,ID) each block has a begining and end flag of 4
         * bytes */
         +4 * 2 *sizeof(int32_t);
        
        MPI_File_set_size (outfile, filesize);
        
        const MPI_Offset pos_start = // first POS guard
                            sizeof(hdr) + 2*sizeof(int32_t),
                         this_pos_start = // local particles POS data
                            pos_start + sizeof(int32_t) 
                            + 3*sizeof(gadget_serialization::position_type)*partial_count;
                            
        const MPI_Offset pos_end = // last POS guard
                            pos_start + sizeof(int32_t) 
                            + sizeof(gadget_serialization::position_type)*3*npart,
                         this_pos_end = // end of local particles POS data
                            this_pos_start + 3*sizeof(gadget_serialization::position_type)*npart;
        
        const MPI_Offset vel_start = // first VEL guard
                            pos_end + sizeof(int32_t),
                         this_vel_start = // local particles VEL data
                            vel_start + sizeof(int32_t) 
                            + 3*sizeof(gadget_serialization::velocity_type)*partial_count;
                         
        const MPI_Offset vel_end = // last VEL guard
                            vel_start + sizeof(int32_t) 
                            + sizeof(gadget_serialization::velocity_type)*3*npart,
                         this_vel_end = // end of local particles VEL data
                            this_vel_start + 3*sizeof(gadget_serialization::velocity_type)*npart;
        
        const MPI_Offset ids_start = // first IDS guard
                            vel_end + sizeof(int32_t),
                         this_ids_start = // local particles IDS data
                            ids_start + sizeof(int32_t) 
                            + sizeof(gadget_serialization::id_type)*partial_count;
                         
        const MPI_Offset ids_end = // last IDS guard
                            ids_start + sizeof(int32_t) 
                            + sizeof(gadget_serialization::id_type)*npart,
                         this_ids_end = // end of local particles IDS data
                            this_ids_start + sizeof(gadget_serialization::id_type)*npart;
        
        if(cart_com.rank()==0)
        {
            uint32_t blocksize{};
            MPI_Status status;
            
            // write header
            blocksize = sizeof (hdr);
            MPI_File_write_at (outfile, 0, &blocksize, 1, MPI_UNSIGNED,
                               &status);
            MPI_File_write_at (outfile, sizeof (uint32_t), &hdr, sizeof (hdr),
                               MPI_BYTE, &status);
            MPI_File_write_at (outfile, sizeof (hdr) + sizeof (uint32_t),
                               &blocksize, 1, MPI_UNSIGNED, &status);
            
            // TODO: note that blocksize is not correct if npart >= 2^32
            // even for smaller values of npart the blocksize fails to fit a 4 byte integer.
            // for example if npart = 711^3 then POS blocksize = 4313105172 which is larger than
            // 2^32-1
            
            // write pos guards
            blocksize = 3 * sizeof (gadget_serialization::position_type) * hdr.npart[1];
            MPI_File_write_at (outfile, pos_start,
                               &blocksize, 1, MPI_UNSIGNED, &status);
            MPI_File_write_at (outfile, pos_end,
                               &blocksize, 1, MPI_UNSIGNED, &status);
            
            // write vel guards
            blocksize = 3 * sizeof (gadget_serialization::velocity_type) * hdr.npart[1];
            MPI_File_write_at (outfile, vel_start,
                               &blocksize, 1, MPI_UNSIGNED, &status);
            MPI_File_write_at (outfile, vel_end,
                               &blocksize, 1, MPI_UNSIGNED, &status);
            
            // write ids guards
            blocksize = sizeof(gadget_serialization::id_type) * hdr.npart[1];
            MPI_File_write_at (outfile, ids_start,
                               &blocksize, 1, MPI_UNSIGNED, &status);
            MPI_File_write_at (outfile, ids_end,
                               &blocksize, 1, MPI_UNSIGNED, &status);
        }
        
        std::vector<gadget_serialization::position_type> posdata;
        std::vector<gadget_serialization::velocity_type> veldata;
        std::vector<gadget_serialization::id_type> IDs;
        
        
        auto flush_to_file = [&]()
        {
            const double inverse_a = 1.0 / hdr.time;
            
            MPI_Status status;
            
            // convert the positions from [Boxsize] to [kpc/h]
            for(auto & pos : posdata)
                pos *= hdr.BoxSize;
            
            // convert the positions from [Boxsize] to [kpc/h]
            // FIXME: vel is not actually velocity, neither momentum, it becomes
            // momentum per unit of mass divided by the scale factor:
            // vel (snapshot) = momentum / mass / a
            for(auto & vel : veldata)
            {
                vel *= cosmology::C_SPEED_OF_LIGHT;
                
                // if we want to save the momentum/mass for particles we comment
                // this line. In an homogeneous universe for non-relativistic
                // species this momentum/mass/a is equal to dx/dtau where x is
                // comoving position and tau is conformal time.
                vel *= inverse_a;
                    
            }
            MPI_File_write_at(outfile,
                              this_pos_start,
                              posdata.data(),
                              posdata.size()*sizeof(decltype(posdata)::value_type),
                              MPI_BYTE,&status);
            
            MPI_File_write_at(outfile,
                              this_vel_start,
                              veldata.data(),
                              veldata.size()*sizeof(decltype(veldata)::value_type),
                              MPI_BYTE,&status);
            
            MPI_File_write_at(outfile,
                              this_ids_start,
                              IDs.data(),
                              IDs.size()*sizeof(decltype(IDs)::value_type),
                              MPI_BYTE,&status);
            
            posdata.clear();
            veldata.clear();
            IDs.clear();
        };
        
        posdata.reserve(3*PCLBUFFER);
        veldata.reserve(3*PCLBUFFER);
        IDs.reserve(PCLBUFFER);
        
        this->for_each(
            [&](const particle& part, const LATfield2::Site& x)
            {
                if(select_function(part,x))
                {
                    // save selected particles into buffer
                    std::copy(std::begin(part.pos),
                              std::end(part.pos),
                              std::back_inserter(posdata));
                              
                    std::copy(std::begin(part.momentum),
                              std::end(part.momentum),
                              std::back_inserter(veldata));
                       
                    IDs.push_back(part.ID);
                }
                
                // buffer is full, flush to file
                if(IDs.size()>=PCLBUFFER)
                {
                    flush_to_file();
                }
            }
        );
        
        if(IDs.size()>0)
            flush_to_file();
        
    }

    
    // void saveGadget2 (std::string filename, gadget2_header &hdr,
    //                   lightcone_geometry &lightcone, double dist, double dtau,
    //                   double dtau_old, double dadtau,
    //                   double vertex[MAX_INTERSECTS][3], const int vertexcount,
    //                   std::set<long> &IDbacklog, std::set<long> &IDprelog,
    //                   Field<Real> *phi, const int tracer_factor = 1);
    
    void loadGadget2 (std::string filename, gadget2_header &hdr);
    
    void update_mass()
    {
        this->for_each(
            [&](particle& part,const LATfield2::Site&)
            {
                part.mass = this->parts_info()->mass;
            });
    }
    
    auto const & communicator()const
    {
        return ::LATfield2::parallel.cartesian_communicator();
    }
    
    auto const & cartesian_communicator()const
    {
        return ::LATfield2::parallel.cartesian_communicator();
    }
};

}
#endif
