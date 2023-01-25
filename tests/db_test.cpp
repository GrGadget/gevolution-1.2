#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
namespace mpi = boost::mpi;


#include <gevolution/distributed_database.hpp>
#include <map>
#include <random>
#include <iostream>
#include <vector>

bool is_prime(const int x)
{
    if(x==0 || x==1)return false;
    if(x==2 || x==3)return true;
    
    for(int i=2;i*i<=x;++i)
        if(x % i == 0)
            return false;
    
    return true;
}

bool oracle(const int x)
{
    return is_prime(x);
}

int main()
{
    
    mpi::environment env;
    mpi::communicator world;
    
    std::default_random_engine gen(world.rank()+10);
    gen();
    std::uniform_int_distribution<int> random_n(0,1000);
    
    
    auto dest_F = [world](const int x){
        return x % world.size();
    };
    
    gevolution::distributed_database<int,bool,decltype(dest_F)> db(world,dest_F);
    
    std::vector< int > my_numbers;
    
    for(int i=0;i<10;++i)
    {
        my_numbers.push_back( random_n(gen) );
    }
    
    std::vector< std::pair<int,bool> > my_input;
    for(auto x: my_numbers)
    {
        my_input.push_back({x,oracle(x)});
    }
    
    db.insert(my_input);
    
    auto my_query = db.query(my_numbers);
    
    for(auto [k,v] : my_query )
    {
        assert(v==oracle(k));
    }
    
    for(int i=0;i<world.size();++i)
    {
        if(i==world.rank())
        {
            std::cout << "proc #"<<i<<"\n";
            std::cout << " have numbers: ";
            for(auto x: my_numbers)
                std::cout << x << " ";
            std::cout<<std::endl;
            
            std::cout << " they correspond to proc: ";
            for(auto x: my_numbers)
                std::cout << dest_F(x) << " ";
            std::cout<<std::endl;
            
            std::cout << " oracle: ";
            for(auto [x,v]: my_query)
                std::cout <<"(" << x << ", " << oracle(x) << ") ";
            std::cout<<std::endl;
            
            
            std::cout << " query: ";
            for(auto [x,v]: my_query)
                std::cout <<"(" << x << ", " << v << ") ";
            std::cout<<std::endl;
            
            world.barrier();
        }
    }
    
    return 0;
}


