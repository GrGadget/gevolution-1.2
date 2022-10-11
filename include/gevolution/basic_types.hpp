#pragma once

#include "gevolution/config.h"
#include "LATfield2.hpp"

namespace gevolution
{
    using LATfield2::Real;
    typedef LATfield2::Imag Cplx;
    
    // type for the serialization of particles into Gadget snapshots
    // TODO: enforce this at compile time of the library of leave it as a template?
    typedef float position_type;
    typedef float velocity_type;
    typedef int64_t id_type;
    
    struct type_checks
    {
        static_assert(sizeof(LATfield2::Real)==8); // double precision on floating point
        static_assert(sizeof(id_type)==8); // serialization of IDs up to 2^64
        static_assert(sizeof(LATfield2::part_simple::ID)==8); // internal IDs up to 2^64
    };
}
