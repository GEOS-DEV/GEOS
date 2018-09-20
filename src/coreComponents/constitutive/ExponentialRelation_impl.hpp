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

#include "RAJA/RAJA.hpp"
#include "RAJA/index/RangeSegment.hpp"

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


template<typename IndexType, typename RealType>
ExponentialRelation<IndexType, RealType>::ExponentialRelation()
  : ExponentialRelation( ExponentApproximationType::Full )
{}

template<typename IndexType, typename RealType>
ExponentialRelation<IndexType, RealType>::ExponentialRelation( ExponentApproximationType type )
  : ExponentialRelation( type, RealType( 0 ), RealType( 1 ), RealType( 1 ))
{}

template<typename IndexType, typename RealType>
ExponentialRelation<IndexType, RealType>::ExponentialRelation( ExponentApproximationType type,
                                                               RealType x0, RealType y0, RealType alpha )
{
  SetParameters( type, x0, y0, alpha );
}

template<typename IndexType, typename RealType>
void ExponentialRelation<IndexType, RealType>::SetApproximationType( ExponentApproximationType type )
{
  m_approximationType = type;
}

template<typename IndexType, typename RealType>
void ExponentialRelation<IndexType, RealType>::SetCoefficients( RealType x0, RealType y0, RealType alpha )
{
  m_x0 = x0;
  m_y0 = y0;
  m_alpha = alpha;
}

template<typename IndexType, typename RealType>
void ExponentialRelation<IndexType, RealType>::SetParameters( ExponentApproximationType type,
                                                              RealType x0, RealType y0, RealType alpha )
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

  inline static void Compute2( const T & x0, const T & y0, const T & alpha, const T & x, T & y, T & dy_dx )
  {
    y = y0 * exp( alpha * (x - x0));
    dy_dx = alpha * y;
  }

  inline static void Inverse( const T & x0, const T & y0, const T & alpha, const T & y, T & x )
  {
    x = x0 + log( y / y0 ) / alpha;
  }

  inline static void Inverse2( const T & x0, const T & y0, const T & alpha, const T & y, T & x, T & dx_dy )
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

  inline static void Compute2( const T & x0, const T & y0, const T & alpha, const T & x, T & y, T & dy_dx )
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

  inline static void Inverse2( const T & x0, const T & y0, const T & alpha, const T & y, T & x, T & dx_dy )
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

  inline static void Compute2( const T & x0, const T & y0, const T & alpha, const T & x, T & y, T & dy_dx )
  {
    y = y0 * (T( 1.0 ) + alpha * (x - x0));
    dy_dx = alpha * y0;
  }

  inline static void Inverse( const T & x0, const T & y0, const T & alpha, const T & y, T & x )
  {
    x = x0 + (y / y0 - 1) / alpha;
  }

  inline static void Inverse2( const T & x0, const T & y0, const T & alpha, const T & y, T & x, T & dx_dy )
  {
    const T alpha_inv = T( 1.0 ) / alpha;
    x = x0 + alpha_inv * (y / y0 - 1);
    dx_dy = alpha_inv / y0;
  }
};

template<typename Policy, typename I, typename Lambda>
void forall_compute( I first, I last, Lambda lambda )
{
  RAJA::RangeSegment seg( first, last );
  RAJA::forall<Policy>( seg, [=] ( I index ) mutable -> void
        {
          lambda( index );
        } );
}

}

template<typename IndexType, typename RealType>
void ExponentialRelation<IndexType, RealType>::Compute( const RealType & x, RealType & y ) const
{
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&]( auto const eat ) -> void
      {
        detail::ExponentialCompute<RealType, eat>::Compute( m_x0, m_y0, m_alpha, x, y );
      } );
}

template<typename IndexType, typename RealType>
void ExponentialRelation<IndexType, RealType>::Inverse( const RealType & y, RealType & x ) const
{
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&]( auto const eat ) -> void
      {
        detail::ExponentialCompute<RealType, eat>::Inverse( m_x0, m_y0, m_alpha, y, x );
      } );
}

template<typename IndexType, typename RealType>
template<typename Policy>
void ExponentialRelation<IndexType, RealType>::BatchCompute( IndexType size, ConstPtrType x_ptr, PtrType y_ptr ) const
{
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&]( auto const eat ) -> void
      {
        detail::forall_compute<Policy>( 0, size, [&] ( IndexType i ) -> void
        {
          detail::ExponentialCompute<RealType, decltype(eat)::value >::
          Compute( m_x0, m_y0, m_alpha, x_ptr[i], y_ptr[i] );
        } );
      } );
}

template<typename IndexType, typename RealType>
template<typename Policy>
void ExponentialRelation<IndexType, RealType>::BatchInverse( IndexType size, ConstPtrType y_ptr, PtrType x_ptr ) const
{
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&]( auto const eat ) -> void
      {
        detail::forall_compute<Policy>( 0, size, [&] ( IndexType i ) -> void
        {
          detail::ExponentialCompute<RealType, decltype(eat)::value >::
          Inverse( m_x0, m_y0, m_alpha, y_ptr[i], x_ptr[i] );
        } );
      } );
}

template<typename IndexType, typename RealType>
void ExponentialRelation<IndexType, RealType>::Compute( const RealType & x, RealType & y, RealType & dy_dx ) const
{
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&]( auto const eat ) -> void
      {
        detail::ExponentialCompute<RealType, decltype(eat)::value >::
        Compute2( m_x0, m_y0, m_alpha, x, y, dy_dx );
      } );
}

template<typename IndexType, typename RealType>
void ExponentialRelation<IndexType, RealType>::Inverse( const RealType & y, RealType & x, RealType & dx_dy ) const
{
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&]( auto const eat ) -> void
      {
        detail::ExponentialCompute<RealType, decltype(eat)::value >::
        Inverse2( m_x0, m_y0, m_alpha, y, x, dx_dy );
      } );
}

template<typename IndexType, typename RealType>
template<typename Policy>
void ExponentialRelation<IndexType, RealType>::BatchCompute( IndexType size, ConstPtrType x_ptr,
                                                             PtrType y_ptr, PtrType dy_dx_ptr ) const
{
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&]( auto const eat ) -> void
      {
        detail::forall_compute<Policy>( 0, size, [&] ( IndexType i ) -> void
        {
          detail::ExponentialCompute<RealType, decltype(eat)::value >::
          Compute2( m_x0, m_y0, m_alpha, x_ptr[i], y_ptr[i], dy_dx_ptr[i] );
        } );
      } );
}

template<typename IndexType, typename RealType>
template<typename Policy>
void ExponentialRelation<IndexType, RealType>::BatchInverse( IndexType size, ConstPtrType y_ptr,
                                                             PtrType x_ptr, PtrType dx_dy_ptr ) const
{
  ExponentApproximationTypeSwitchBlock( m_approximationType, [&]( auto const eat ) -> void
      {
        detail::forall_compute<Policy>( 0, size, [&] ( IndexType i ) -> void
        {
          detail::ExponentialCompute<RealType, decltype(eat)::value >::
          Inverse2( m_x0, m_y0, m_alpha, y_ptr[i], x_ptr[i], dx_dy_ptr[i] );
        } );
      } );
}

} // namespace constitutive

} // namespace geosx


#endif // GEOSX_EXPONENTIALRELATION_IMPL_HPP_
