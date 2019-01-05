/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


/**
 * @file ExponentialRelations_impl.hpp
 */

#ifndef GEOSX_EXPONENTIALRELATION_IMPL_HPP_
#define GEOSX_EXPONENTIALRELATION_IMPL_HPP_

#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#include <cmath>

namespace geosx
{

namespace constitutive
{

template< typename LAMBDA >
static auto ExponentApproximationTypeSwitchBlock( ExponentApproximationType const eat,
                                                  LAMBDA&& lambda )
{
  switch( eat )
  {
    case ExponentApproximationType::Full:
    {
      return lambda( ExponentApproximationTypeWrapper<ExponentApproximationType::Full>() );
    }
    case ExponentApproximationType::Quadratic:
    {
      return lambda( ExponentApproximationTypeWrapper<ExponentApproximationType::Quadratic>() );
    }

    case ExponentApproximationType::Linear:
    {
      return lambda( ExponentApproximationTypeWrapper<ExponentApproximationType::Linear>() );
    }
    default:
    {
      GEOS_ERROR( "ExponentApproximationTypeSwitchBlock() ExponentApproximationType is invalid!" );
    }
  }
}


template<typename T>
ExponentialRelation<T>::ExponentialRelation()
  : ExponentialRelation( ExponentApproximationType::Full )
{}

template<typename T>
ExponentialRelation<T>::ExponentialRelation( ExponentApproximationType type )
  : ExponentialRelation( type, T( 0 ), T( 1 ), T( 1 ))
{}

template<typename T>
ExponentialRelation<T>::ExponentialRelation( ExponentApproximationType type, T x0, T y0, T alpha )
{
  SetParameters( type, x0, y0, alpha );
}

template<typename T>
void ExponentialRelation<T>::SetApproximationType( ExponentApproximationType type )
{
  m_approximationType = type;
}

template<typename T>
void ExponentialRelation<T>::SetCoefficients( T x0, T y0, T alpha )
{
  m_x0 = x0;
  m_y0 = y0;
  m_alpha = alpha;
}

template<typename T>
void ExponentialRelation<T>::SetParameters( ExponentApproximationType type, T x0, T y0, T alpha )
{
  SetApproximationType( type );
  SetCoefficients( x0, y0, alpha );
}

namespace detail
{

using std::exp;
using std::log;
using std::sqrt;

template<typename T, ExponentApproximationType ApproxType>
struct ExponentialCompute
{};

template<typename T>
struct ExponentialCompute<T, ExponentApproximationType::Full>
{
  inline static void Compute( const T & x0, const T & y0, const T & alpha, const T & x, T & y )
  {
    y = y0 * exp( alpha * (x - x0));
  }

  inline static void Compute( const T & x0, const T & y0, const T & alpha, const T & x, T & y, T & dy_dx )
  {
    y = y0 * exp( alpha * (x - x0));
    dy_dx = alpha * y;
  }

  inline static void Inverse( const T & x0, const T & y0, const T & alpha, const T & y, T & x )
  {
    x = x0 + log( y / y0 ) / alpha;
  }

  inline static void Inverse( const T & x0, const T & y0, const T & alpha, const T & y, T & x, T & dx_dy )
  {
    const T alpha_inv = T( 1.0 ) / alpha;
    x = x0 + alpha_inv * log( y / y0 );
    dx_dy = alpha_inv / y;
  }
};

template<typename T>
struct ExponentialCompute<T, ExponentApproximationType::Quadratic>
{
  inline static void Compute( const T & x0, const T & y0, const T & alpha, const T & x, T & y )
  {
    const T z = T( 1.0 ) + alpha * (x - x0);
    y = y0 / T( 2.0 ) * (T( 1.0 ) + z * z);
  }

  inline static void Compute( const T & x0, const T & y0, const T & alpha, const T & x, T & y, T & dy_dx )
  {
    const T z = T( 1.0 ) + alpha * (x - x0);
    y = y0 / T( 2.0 ) * (T( 1.0 ) + z * z);
    dy_dx = alpha * y0 * z;
  }

  inline static void Inverse( const T & x0, const T & y0, const T & alpha, const T & y, T & x )
  {
    const T z = sqrt( T( 2.0 ) * y / y0 - T( 1.0 ));
    x = x0 + (z - 1) / alpha;
  }

  inline static void Inverse( const T & x0, const T & y0, const T & alpha, const T & y, T & x, T & dx_dy )
  {
    const T alpha_inv = T( 1.0 ) / alpha;
    const T z = sqrt( T( 2.0 ) * y / y0 - T( 1.0 ));
    x = x0 + alpha_inv * (z - 1);
    dx_dy = alpha_inv / (y0 * z);
  }
};

template<typename T>
struct ExponentialCompute<T, ExponentApproximationType::Linear>
{
  inline static void Compute( const T & x0, const T & y0, const T & alpha, const T & x, T & y )
  {
    y = y0 * (T( 1.0 ) + alpha * (x - x0));
  }

  inline static void Compute( const T & x0, const T & y0, const T & alpha, const T & x, T & y, T & dy_dx )
  {
    y = y0 * (T( 1.0 ) + alpha * (x - x0));
    dy_dx = alpha * y0;
  }

  inline static void Inverse( const T & x0, const T & y0, const T & alpha, const T & y, T & x )
  {
    x = x0 + (y / y0 - 1) / alpha;
  }

  inline static void Inverse( const T & x0, const T & y0, const T & alpha, const T & y, T & x, T & dx_dy )
  {
    const T alpha_inv = T( 1.0 ) / alpha;
    x = x0 + alpha_inv * (y / y0 - 1);
    dx_dy = alpha_inv / y0;
  }
};

}

template<typename T>
void ExponentialRelation<T>::Compute( const T & x, T & y ) const
{
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&] ( auto const eat )
  {
    detail::ExponentialCompute<T, decltype(eat)::value>::Compute( m_x0, m_y0, m_alpha, x, y );
  } );
}

template<typename T>
void ExponentialRelation<T>::Inverse( const T & y, T & x ) const
{
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&] ( auto const eat )
  {
    detail::ExponentialCompute<T, decltype(eat)::value>::Inverse( m_x0, m_y0, m_alpha, y, x );
  } );
}

template<typename T>
template<typename Policy>
void ExponentialRelation<T>::Compute( arrayView1d<T const> const & x, arrayView1d<T> & y ) const
{
  GEOS_ERROR_IF( y.size(0) != x.size(0), "Array size mismatch" );
  T x0 = m_x0;
  T y0 = m_y0;
  T alpha = m_alpha;
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&] ( auto const eat )
  {
    forall_in_range<Policy>( 0, x.size(0), GEOSX_LAMBDA ( auto i )
    {
      detail::ExponentialCompute<T, decltype(eat)::value>::Compute( x0, y0, alpha, x[i], y[i] );
    } );
  } );
}

template<typename T>
template<typename Policy>
void ExponentialRelation<T>::Inverse( arrayView1d<T const> const & y, arrayView1d<T> & x ) const
{
  GEOS_ERROR_IF( x.size(0) != y.size(0), "Array size mismatch" );
  T x0 = m_x0;
  T y0 = m_y0;
  T alpha = m_alpha;
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&] ( auto const eat )
  {
    forall_in_range<Policy>( 0, y.size(0), GEOSX_LAMBDA ( auto i )
    {
      detail::ExponentialCompute<T, decltype(eat)::value>::Inverse( x0, y0, alpha, y[i], x[i] );
    } );
  } );
}

template<typename T>
template<typename Policy>
void ExponentialRelation<T>::Compute( arrayView1d<T const> const & x, arrayView2d<T> & y ) const
{
  GEOS_ERROR_IF( y.size(0) != x.size(0), "Array size mismatch" );
  T x0 = m_x0;
  T y0 = m_y0;
  T alpha = m_alpha;
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&] ( auto const eat )
  {
    auto const numJ = y.size(1);
    forall_in_range<Policy>( 0, x.size(0), GEOSX_LAMBDA ( auto i )
    {
      for (localIndex j = 0; j < numJ; ++j)
      {
        detail::ExponentialCompute<T, decltype(eat)::value>::Compute( x0, y0, alpha, x[i], y[i][j] );
      }
    } );
  } );
}

template<typename T>
template<typename Policy>
void ExponentialRelation<T>::Inverse( arrayView1d<T const> const & y, arrayView2d<T> & x ) const
{
  GEOS_ERROR_IF( x.size(0) != y.size(0), "Array size mismatch" );
  T x0 = m_x0;
  T y0 = m_y0;
  T alpha = m_alpha;
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&] ( auto const eat )
  {
    auto const numJ = x.size(1);
    forall_in_range<Policy>( 0, y.size(0), GEOSX_LAMBDA ( auto i )
    {
      for (localIndex j = 0; j < numJ; ++j)
      {
        detail::ExponentialCompute<T, decltype(eat)::value>::Inverse( x0, y0, alpha, y[i], x[i][j] );
      }
    } );
  } );
}

template<typename T>
template<typename Policy>
void ExponentialRelation<T>::Compute( arrayView2d<T const> const & x, arrayView2d<T> & y ) const
{
  GEOS_ERROR_IF( y.size(0) != x.size(0), "Array size mismatch" );
  GEOS_ERROR_IF( y.size(1) != x.size(1), "Array size mismatch" );
  T x0 = m_x0;
  T y0 = m_y0;
  T alpha = m_alpha;
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&] ( auto const eat )
  {
    auto const numJ = x.size(1);
    forall_in_range<Policy>( 0, x.size(0), GEOSX_LAMBDA ( auto i )
    {
      for (localIndex j = 0; j < numJ; ++j)
      {
        detail::ExponentialCompute<T, decltype(eat)::value>::Compute( x0, y0, alpha, x[i][j], y[i][j] );
      }
    } );
  } );
}

template<typename T>
template<typename Policy>
void ExponentialRelation<T>::Inverse( arrayView2d<T const> const & y, arrayView2d<T> & x ) const
{
  GEOS_ERROR_IF( x.size(0) != y.size(0), "Array size mismatch" );
  GEOS_ERROR_IF( x.size(1) != y.size(1), "Array size mismatch" );
  T x0 = m_x0;
  T y0 = m_y0;
  T alpha = m_alpha;
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&] ( auto const eat )
  {
    auto const numJ = y.size(1);
    forall_in_range<Policy>( 0, y.size(0), GEOSX_LAMBDA ( auto i )
    {
      for (localIndex j = 0; j < numJ; ++j)
      {
        detail::ExponentialCompute<T, decltype(eat)::value>::Inverse( x0, y0, alpha, y[i][j], x[i][j] );
      }
    } );
  } );
}

template<typename T>
void ExponentialRelation<T>::Compute( const T & x, T & y, T & dy_dx ) const
{
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&] ( auto const eat )
  {
    detail::ExponentialCompute<T, decltype(eat)::value>::Compute( m_x0, m_y0, m_alpha, x, y, dy_dx );
  } );
}

template<typename T>
void ExponentialRelation<T>::Inverse( const T & y, T & x, T & dx_dy ) const
{
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&] ( auto const eat )
  {
    detail::ExponentialCompute<T, decltype(eat)::value>::Inverse( m_x0, m_y0, m_alpha, y, x, dx_dy );
  } );
}

template<typename T>
template<typename Policy>
void ExponentialRelation<T>::Compute( arrayView1d<T const> const & x, arrayView1d<T> & y, arrayView1d<T> & dy_dx ) const
{
  GEOS_ERROR_IF( y.size(0) != x.size(0) || dy_dx.size(0) != x.size(0), "Array size mismatch" );
  T x0 = m_x0;
  T y0 = m_y0;
  T alpha = m_alpha;
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&] ( auto const eat )
  {
    forall_in_range<Policy>( 0, x.size(0), GEOSX_LAMBDA ( auto i )
    {
      detail::ExponentialCompute<T, decltype(eat)::value>::Compute( x0, y0, alpha, x[i], y[i], dy_dx[i] );
    } );
  } );
}

template<typename T>
template<typename Policy>
void ExponentialRelation<T>::Inverse( arrayView1d<T const> const & y, arrayView1d<T> & x, arrayView1d<T> & dx_dy ) const
{
  GEOS_ERROR_IF( x.size(0) != y.size(0) || dx_dy.size(0) != y.size(0), "Array size mismatch" );
  T x0 = m_x0;
  T y0 = m_y0;
  T alpha = m_alpha;
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&] ( auto const eat )
  {
    forall_in_range<Policy>( 0, y.size(0), GEOSX_LAMBDA ( auto i )
    {
      detail::ExponentialCompute<T, decltype(eat)::value>::Inverse( x0, y0, alpha, y[i], x[i], dx_dy[i] );
    } );
  } );
}

template<typename T>
template<typename Policy>
void ExponentialRelation<T>::Compute( arrayView1d<T const> const & x, arrayView2d<T> & y, arrayView2d<T> & dy_dx ) const
{
  GEOS_ERROR_IF( y.size(0) != x.size(0) || dy_dx.size(0) != x.size(0), "Array size mismatch" );
  GEOS_ERROR_IF( dy_dx.size(1) != y.size(1), "Array size mismatch" );
  T x0 = m_x0;
  T y0 = m_y0;
  T alpha = m_alpha;
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&] ( auto const eat )
  {
    auto const numJ = y.size(1);
    forall_in_range<Policy>( 0, x.size(0), GEOSX_LAMBDA ( auto i )
    {
      for (localIndex j = 0; j < numJ; ++j)
      {
        detail::ExponentialCompute<T, decltype(eat)::value>::Compute( x0, y0, alpha, x[i], y[i][j], dy_dx[i][j] );
      }
    } );
  } );
}

template<typename T>
template<typename Policy>
void ExponentialRelation<T>::Inverse( arrayView1d<T const> const & y, arrayView2d<T> & x, arrayView2d<T> & dx_dy ) const
{
  GEOS_ERROR_IF( x.size(0) != y.size(0) || dx_dy.size(0) != y.size(0), "Array size mismatch" );
  GEOS_ERROR_IF( dx_dy.size(1) != x.size(1), "Array size mismatch" );
  T x0 = m_x0;
  T y0 = m_y0;
  T alpha = m_alpha;
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&] ( auto const eat )
  {
    auto const numJ = x.size(1);
    forall_in_range<Policy>( 0, y.size(0), GEOSX_LAMBDA ( auto i )
    {
      for (localIndex j = 0; j < numJ; ++j)
      {
        detail::ExponentialCompute<T, decltype(eat)::value>::Inverse( x0, y0, alpha, y[i], x[i][j], dx_dy[i][j] );
      }
    } );
  } );
}

template<typename T>
template<typename Policy>
void ExponentialRelation<T>::Compute( arrayView2d<T const> const & x, arrayView2d<T> & y, arrayView2d<T> & dy_dx ) const
{
  GEOS_ERROR_IF( y.size(0) != x.size(0) || dy_dx.size(0) != x.size(0), "Array size mismatch" );
  GEOS_ERROR_IF( y.size(1) != x.size(1) || dy_dx.size(1) != x.size(1), "Array size mismatch" );
  T x0 = m_x0;
  T y0 = m_y0;
  T alpha = m_alpha;
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&] ( auto const eat )
  {
    auto const numJ = x.size(1);
    forall_in_range<Policy>( 0, x.size(0), GEOSX_LAMBDA ( auto i )
    {
      for (localIndex j = 0; j < numJ; ++j)
      {
        detail::ExponentialCompute<T, decltype(eat)::value>::Compute( x0, y0, alpha, x[i][j], y[i][j], dy_dx[i][j] );
      }
    } );
  } );
}

template<typename T>
template<typename Policy>
void ExponentialRelation<T>::Inverse( arrayView2d<T const> const & y, arrayView2d<T> & x, arrayView2d<T> & dx_dy ) const
{
  GEOS_ERROR_IF( x.size(0) != y.size(0) || dx_dy.size(0) != y.size(0), "Array size mismatch" );
  GEOS_ERROR_IF( x.size(1) != y.size(1) || dx_dy.size(1) != y.size(1), "Array size mismatch" );
  T x0 = m_x0;
  T y0 = m_y0;
  T alpha = m_alpha;
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&] ( auto const eat )
  {
    auto const numJ = y.size(1);
    forall_in_range<Policy>( 0, y.size(0), GEOSX_LAMBDA ( auto i )
    {
      for (localIndex j = 0; j < numJ; ++j)
      {
        detail::ExponentialCompute<T, decltype(eat)::value>::Inverse( x0, y0, alpha, y[i][j], x[i][j], dx_dy[i][j] );
      }
    } );
  } );
}

} // namespace constitutive

} // namespace geosx


#endif // GEOSX_EXPONENTIALRELATION_IMPL_HPP_
