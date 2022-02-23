#pragma once

#include <complex> // norm
#include <vector>  // vector
#include <utility> // pair
#include <exception> // runtime_error
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/inplace.hpp>
#include <boost/serialization/utility.hpp>
#include "LATfield2.hpp"

namespace gevolution
{
    namespace detail
    {
        // define this because boost::mpi does not support lambdas yet
        template<typename T1, typename T2>
        struct my_pair_sum
        {
            using value_type = typename std::pair<T1,T2>;
            
            value_type operator()(const value_type& x, const value_type& y)
            {
                return value_type{x.first+y.first,x.second+y.second};
            }
        };
    } // namespace detail
    // power spectrum of a scalar field
    template<class complex_field_type>
    auto power_spectrum(const complex_field_type& F)
    {
        using complex_type = typename complex_field_type::value_type;
        using real_type = typename complex_type::value_type;
        using std::norm;
        const int N_global = F.lattice().size(1);
        const real_type N6_inv = std::pow( 1.0 / N_global , 6 );
        const int k_nyquist = (N_global - 1)/2;
        std::vector< std::pair<int,real_type> > pw(k_nyquist+1,{0,0});
        
        auto signed_mode = [k_nyquist,N_global](int n)
        {
            return n <= k_nyquist ? n : n-N_global;
        };
        
        // avoid doubly counting modes, this happens because we have a source
        // field which is real and the complex FFT has been done with FFTW
        // scheme of R2C, hence some modes are omited because of the halcomplex
        // symmetry, while other modes with complex conjugate counterpart are
        // still present.
        auto unique_mode = [](std::array<int,3> k_mode)
        {
            bool answer = true;
            for(auto k: k_mode)
            {
                if(k<0)
                    answer = false;
                else if(k>0)
                    break;
            }
            
            return answer;
        };
        
        if(F.components()>1)
            throw std::runtime_error(
                "gevolution::power_spectrum only works for scalar fields");
        
        // accumulate modes on local grid
        F.for_each( 
            [&pw,&signed_mode,&unique_mode](
                const complex_type& value, const LATfield2::Site & x)
            {
                std::array<int,3> k_modes;
                double global_mode=0;
                for(int i=0;i<3;++i)
                {
                    int k_i = signed_mode( x.coord(i) );
                    k_modes[i] = k_i;
                    global_mode += k_i*k_i;
                }
                size_t index = std::floor( std::sqrt(global_mode) + 0.5 );
                
                if(unique_mode(k_modes) and index<pw.size())
                {
                    auto& [count,sum] = pw[index];
                    ++count;
                    using std::norm;
                    auto z2 = norm(value);
                    sum += z2;
                    
                    // if(index==1)
                    // {
                    //     std::cerr << "gevolution::power_spectrum mode: " 
                    //         << k_modes[0] << " " << k_modes[1] << " " << k_modes[2]
                    //         << " value = " << z2
                    //         << '\n';
                    // }
                }
            });
        
        const boost::mpi::communicator& com = LATfield2::parallel.my_comm;
        ::boost::mpi::all_reduce(com,
            ::boost::mpi::inplace_t< std::pair<int,real_type>* >(pw.data()),pw.size(),
            detail::my_pair_sum<int,real_type>());
        std::vector<real_type> average(pw.size());
        for(auto i=0U;i<pw.size();++i)
        {
            if(pw[i].first>0)
                pw[i].second /= pw[i].first;
            pw[i].second *= N6_inv;
        }
        return pw;
    }
} // namespace gevolution
