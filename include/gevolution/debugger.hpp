#pragma once

#include "gevolution/config.h"
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
#include <filesystem>

namespace gevolution
{

struct part_data
{
    unsigned long long ID;
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
    const double Pos_physical{1},Acc_physical{1};
    
  public:
    debugger_t (boost::mpi::communicator _com, std::string _fname,
        double Pos_fac,double Acc_fac)
        : com{ _com }, fname{ _fname },
        Pos_physical{Pos_fac}, Acc_physical{Acc_fac}
    {
        if (com.rank()==root)
            std::filesystem::remove(fname);
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
    void append (unsigned long long id, std::array<double, 3> Pos, std::array<double, 3> Acc)
    {
        for(auto &x: Pos) x *= Pos_physical;
        for(auto &a: Acc) a *= Acc_physical;
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
