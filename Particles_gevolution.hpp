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
#include "particle_handler.hpp"
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

/*
    The indivual particle data structure in gevolution is 'particle',
    defined below.
*/
struct particle : LATfield2::part_simple
{
    LATfield2::Real acc[3];
};

class Particles_gevolution : 
    public 
    LATfield2::Particles<
        particle, 
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

// /* Specialization of particle_handler for Gevolution's type of particle */
// class Particles_gevolution_handler : public particle_handler
// {
//   Particles_gevolution &P;
// 
//  public:
//   Particles_gevolution_handler(Particles_gevolution &ref_P) : 
//     P{ref_P} {}
// 
//   size_t size() const override { return P.NumPart; }
// 
//   std::array<long long int, 3> get_position(int idx) const override
//   {
//     int i = Sp.get_active_index(idx);
//     return {Sp.P[i].IntPos[0], Sp.P[i].IntPos[1], Sp.P[i].IntPos[2]};
//   }
//   double get_mass(int i) const override { return Sp.P[i].getMass(); }
//   void set_acceleration(int idx, std::array<double, 3> A) const override
//   {
//     int i             = Sp.get_active_index(idx);
//     Sp.P[i].GravPM[0] = A[0];
//     Sp.P[i].GravPM[1] = A[1];
//     Sp.P[i].GravPM[2] = A[2];
//   }
// };


}
#endif
