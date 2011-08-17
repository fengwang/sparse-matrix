#ifndef _BICONJUGATE_GRADIENT_HPP_INCLUDED_SFODIUJ34OIUSFUY84987SFKIU4IOUSFDKLJHERLKJSFLKJSDFLKJSLKFJSKJDFLKJSDKJFDSDFDS
#define _BICONJUGATE_GRADIENT_HPP_INCLUDED_SFODIUJ34OIUSFUY84987SFKIU4IOUSFDKLJHERLKJSFLKJSDFLKJSLKFJSKJDFLKJSDKJFDSDFDS

#include <parallel_organizer.hpp>
#include <sparse_matrix.hpp>

#include <cmath>
#include <cstddef>
#include <cassert>

#include <vector>
#include <valarray>
#include <numeric>
#include <iterator>
#include <algorithm>

namespace sm
{

    template< typename T, typename Array=std::valarray<T>, typename Matrix=sm::sparse_matrix<T> >
    struct biconjugate_method
    {
        typedef biconjugate_method                  self_type;
        typedef Matrix                              matrix_type;
        typedef Array                               array_type;
        typedef typename matrix_type::value_type    value_type;
        typedef std::size_t                         size_type;

        //Ax = b -- A[n][n] x[n] = b[n]
        const matrix_type&                          A;
        array_type&                                 x;
        const array_type&                           b;

        value_type                                  eps;

        biconjugate_method( const matrix_type&      A_, 
                            array_type&             x_, 
                            const array_type&       b_, 
                            const value_type        eps_ = value_type(1.0E-80L) ) : 
                            A(A_), x(x_), b(b_), eps(eps_) 
        {
            assert( A.row() == A.col() );
            assert( A.row() == b.size() );
            assert( b.size() == x.size() );
            assert( eps >= value_type() );
        }

        void
        operator() () const 
        {
            cg_impl();
        }

        private:
        void cg_impl() const 
        {
            const size_type n = b.size();
            array_type r = b - A * x;
            array_type p = r;
            value_type alpha;
            value_type beta;

            value_type tmp1;
            value_type tmp2;

            std::vector<value_type> cache;

            const value_type rem = eps * std::inner_product( &b[0], &b[n], &b[0], value_type() ) / n;

            for(size_type i = 0; i < n; ++i )
            {
                // paralleled <- HOT SPOT
                const array_type&& w     = A * p;

                ////-----------------------------------------------------------------------------------------------
                ////                              << None Parallel Impl >>
                tmp1  = std::inner_product( &r[0], &r[n], &r[0], value_type() );
                ////***********************************************************************************************
                ////                              << Parallel Impl >>
                //cache.clear();
                //parallel_organizer()(  []( value_type* first1, value_type* last1, value_type* first2 ){ return std::inner_product(first1, last1, first2, value_type()); },
                //                       &r[0], &r[n], &r[0], std::back_inserter( cache ) );
                //tmp1 = std::accumulate( cache.begin(), cache.end(), value_type() );
                ////-----------------------------------------------------------------------------------------------

                ////-----------------------------------------------------------------------------------------------
                ////                              << None Parallel Impl >>
                tmp2  = std::inner_product( &p[0], &p[n], &w[0], value_type() );
                ////***********************************************************************************************
                ////                              << Parallel Impl >>
                //cache.clear();
                //parallel_organizer()( do_inner_production_const, &p[0], &p[n], &w[0], std::back_inserter( cache ) );
                //tmp2 = std::accumulate( cache.begin(), cache.end(), value_type() );
                ////-----------------------------------------------------------------------------------------------

                alpha = tmp1 / tmp2;

                if ( isnan(alpha) || isinf(alpha) ) break;

                // needs paralleling
                x    += alpha * p;
                // needs paralleling
                r    -= alpha * w;

                // needs paralleling
                tmp2  = std::inner_product( &r[0], &r[n], &r[0], value_type() );

                if ( tmp2 < rem ) break;

                beta  = tmp2 / tmp1;

                if ( isnan(beta) || isinf(beta) ) break;

                // needs paralleling
                p    *= beta;
                // needs paralleling
                p    += r;
            }
        }


    };

}//namespace feng

#endif//_BICONJUGATE_GRADIENT_HPP_INCLUDED_SFODIUJ34OIUSFUY84987SFKIU4IOUSFDKLJHERLKJSFLKJSDFLKJSLKFJSKJDFLKJSDKJFDSDFDS

