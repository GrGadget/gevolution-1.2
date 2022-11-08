#pragma once

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

#include <cmath>
#include <set>
#include <map>
#include <vector>
#include <utility>

namespace gevolution {


template<typename key_type, typename value_type, typename destination_func_type>
class distributed_database
{
    typedef std::pair<key_type,value_type> data_type;
    
    ::boost::mpi::communicator com;
    destination_func_type destination;
    
    std::map<key_type,value_type> my_values;
    
    public:
    distributed_database(const mpi::communicator& in_com,const destination_func_type dest):
        com{in_com}, destination{dest}
    {}
    
    auto local_query(const key_type x)const
    {
        assert(destination(x)==com.rank());
        return my_values.at(x);
    }
    void local_insert(const key_type x, const value_type v)
    {
        assert(destination(x)==com.rank());
        my_values[x]=v;
    }
    
    void insert( const std::vector<data_type>& batch )
    {
        std::vector< std::vector<data_type> > send(com.size()),recv;
        
        for(auto data: batch)
        {
            const auto proc = destination(data.first);
            send[proc].push_back(data);
        }
        
        mpi::all_to_all(com,send,recv);
        
        for(const auto & v: recv)
            for(const auto &data : v)
                local_insert(data.first,data.second);
    }
    auto query(const std::vector<key_type>& batch)const
    {
        std::vector< std::vector<data_type> > sendrecv_result(com.size());
        
        {
            std::vector< std::vector<key_type> > sendrecv_query(com.size());
            for(auto data: batch)
            {
                const auto proc = destination(data);
                sendrecv_query[proc].push_back(data);
            }
            
            mpi::all_to_all(com,sendrecv_query,sendrecv_query);
            
           for(int proc=0;proc<com.size();++proc)
               for(const auto &key : sendrecv_query[proc])
               {
                   sendrecv_result[proc].push_back( {key,local_query(key)} );
               }
           
           mpi::all_to_all(com,sendrecv_result,sendrecv_result);
            
        } // clean sendrecv_query
        
        std::vector< data_type > result;
        
        for(const auto &v : sendrecv_result)
            for(const auto &d : v)
                result.push_back(d);
        
        return result;
    }
};

} // namespace gevolution

