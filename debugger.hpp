#pragma once

#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace gevolution
{
class debugger_t : public std::stringstream
{
    std::vector<long> ID;
    std::vector<std::array<double, 3> > pos;
    std::vector<std::array<double, 3> > acc;

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
        for (int p = 0; p < com.size (); ++p)
        {
            if (p == com.rank ())
            {
                std::ofstream os (fname, std::ios::binary | std::ios::app);
                int n = ID.size ();
                os << n;
                for (int i = 0; i < n; ++i)
                {
                    os << ID[i] << pos[i][0] << pos[i][1] << pos[i][2]
                       << acc[i][0] << acc[i][1] << acc[i][2];
                }

                ID.clear ();
                pos.clear ();
                acc.clear ();
            }
            com.barrier ();
        }
    }
    void append (long id, std::array<double, 3> Pos, std::array<double, 3> Acc)
    {
        ID.push_back (id);
        pos.emplace_back (Pos);
        acc.emplace_back (Acc);
    }
    ~debugger_t () { flush (); }
};
extern debugger_t *Debugger;
} // namespace gevolution
