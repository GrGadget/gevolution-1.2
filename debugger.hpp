#pragma once

#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

namespace gevolution
{

struct part_data
{
    long ID;
    std::array<double,3> pos;
    std::array<double,3> acc;
    
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & ID;
        ar & pos;
        ar & acc;
    }
};

class debugger_t 
{
    std::vector< part_data > data;

    boost::mpi::communicator com;
    std::string fname;
    const int root{ 0 };

  public:
    debugger_t (boost::mpi::communicator _com, std::string _fname)
        : com{ _com }, fname{ _fname }
    {
    }

    void flush ()
    {
        #ifndef NDEBUG
        for (int p = 0; p < com.size (); ++p)
        {
            if (p == com.rank ())
            {
                std::ofstream os (fname, std::ios::binary | std::ios::app);
                boost::archive::binary_oarchive oa(os);
                oa << data;
                data.clear();
            }
            com.barrier ();
        }
        #endif
    }
    void append (long id, std::array<double, 3> Pos, std::array<double, 3> Acc)
    {
        data.push_back({id,Pos,Acc});
    }
    ~debugger_t () { flush (); }
};

class analizer_t
{
    std::vector<part_data> data;
    std::string fname;
    
    public:
    analizer_t(std::string _fname):
        fname{_fname}
    {
        std::ifstream is (fname, std::ios::binary);
        boost::archive::binary_iarchive ia(is);
        std::vector<part_data> buff;
        while(is)
        {
            ia >> buff;
            std::copy(buff.begin(),buff.end(),std::back_inserter(data));
        }
    }
};

extern debugger_t *Debugger;
} // namespace gevolution
