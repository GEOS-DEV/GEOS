/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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
 * @file ExponentialRelations.hpp
 */

#ifndef GEOSX_EXPONENTIALRELATION_HPP_
#define GEOSX_EXPONENTIALRELATION_HPP_

#include <common/DataTypes.hpp>

#include <cmath>

namespace geosx
{

namespace constitutive
{

/**
 * @enum Enumeration describing available approximations of the exponent
 */
enum class ExponentApproximationType
{
  Full,
  Linear,
  Quadratic
};

/**
 * @class ExponentialRelation
 *
 * Describes a simple exponential relationship with three parameters of the form: \f$ y = y0 * exp(a * (x - x0)) \f$.
 * Optional linear and quadratic Taylor series approximation of the exponent.
 * Implements point and batch direct (y(x)) and inverse (x(y)) compute functions and provides derivatives.
 *
 * @tparam T scalar real-valued type used in computation
 * @tparam EAT the type/order of exponent approximation (linear, quadratic, full)
 */
template<typename T, ExponentApproximationType EAT>
class ExponentialRelation
{
public:

  static constexpr ExponentApproximationType type = EAT;

  // *** constructors ***

  /**
   * @brief Default constructor. Sets \f$ x0 = 0, y0 = 1, alpha = 1 \f$
   */
  explicit ExponentialRelation();

  /**
   * @brief Constructor with arguments
   * @param x0 variable shift
   * @param y0 scaling coefficient
   * @param alpha exponential coefficient
   */
  explicit ExponentialRelation( T x0, T y0, T alpha );

  // *** setters ***

  /**
   * @brief Sets the relation coefficients
   * @param x0 variable shift
   * @param y0 scaling coefficient
   * @param alpha exponential coefficient
   */
  void SetCoefficients( T x0, T y0, T alpha );

  // *** no-derivative computes ***

  /**
   * @brief Compute value
   * @param x
   * @param y
   */
  void Compute( const T & x, T & y ) const;

  /**
   * @brief Compute inverse value
   * @param y
   * @param x
   */
  void Inverse( const T & y, T & x ) const;

  // *** derivative computes ***

  /**
   * @brief Compute value and derivative
   * @param x
   * @param y
   * @param dy_dx
   */
  void Compute( const T & x, T & y, T & dy_dx ) const;

  /**
   * @brief Compute inverse value and derivative
   * @param y
   * @param x
   * @param dx_dy
   */
  void Inverse( const T & y, T & x, T & dx_dy ) const;

private:

  /// variable shift
  T m_x0;

  /// scaling coefficient
  T m_y0;

  /// exponential coefficient
  T m_alpha;

};

template< ExponentApproximationType EAT >
struct ExponentApproximationTypeWrapper
{
  static constexpr ExponentApproximationType value = EAT;
};

/**
 * @brief Function that implements the switchyard for ExponentApproximationType
 * @tparam LAMBDA the type of the user-provided generic lambda
 * @param type type/order of the approximation (linear, quadratic, full)
 * @param lambda user-provided generic lambda to be called with an instance of EAT wrapper type
 */
template<typename T, typename LAMBDA>
void ExponentApproximationTypeSwitchBlock( ExponentApproximationType const type, T const & x0, T const & y0, T const & alpha, LAMBDA && lambda )
{
  switch( type )
  {
    case ExponentApproximationType::Full:
    {
      return lambda( ExponentialRelation<T, ExponentApproximationType::Full>(x0, y0, alpha) );
    }
    case ExponentApproximationType::Quadratic:
    {
      return lambda( ExponentialRelation<T, ExponentApproximationType::Quadratic>(x0, y0, alpha) );
    }
    case ExponentApproximationType::Linear:
    {
      return lambda( ExponentialRelation<T, ExponentApproximationType::Linear>(x0, y0, alpha) );
    }
    default:
    {
      GEOS_ERROR( "ExponentApproximationTypeSwitchBlock() ExponentApproximationType is invalid!" );
    }
  }
}

/**
 * @brief Same as above, but with non-default constructor arguments (relation coefficients)
 * @tparam T scalar real-valued type
 * @tparam LAMBDA type of the user-provided generic lambda
 * @param type runtime choice of ExponentApproximationType
 * @param x0 variable shift
 * @param y0 scaling coefficient
 * @param alpha exponential coefficient
 * @param lambda user-provided generic lambda that will be passed an instance of ExponentialRelation<T,?>
 */
template<typename T, typename LAMBDA>
void makeExponentialRelation( ExponentApproximationType type,
                              T const & x0, T const & y0, T const & alpha,
                              LAMBDA && lambda )
{
  ExponentApproximationTypeSwitchBlock( type, x0, y0, alpha, std::move(lambda));
}

namespace detail
{

// bring into scope for ADL
using std::exp;
using std::log;
using std::sqrt;

/**
 * @brief Specializations of this struct provide compute methods for different exponent approximations
 * @tparam T the data type
 * @tparam EAT the type/order of exponent approximation (linear, quadratic, full)
 * @tparam INVERSE whether to compute the inverse of the exponent
 */
template<typename T, ExponentApproximationType EAT, bool INVERSE = false>
struct ExponentialCompute
{};

template<typename T>
struct ExponentialCompute<T, ExponentApproximationType::Full, false>
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
};

template<typename T>
struct ExponentialCompute<T, ExponentApproximationType::Full, true>
{
  inline static void Compute( const T & x0, const T & y0, const T & alpha, const T & y, T & x )
  {
    x = x0 + log( y / y0 ) / alpha;
  }

  inline static void Compute( const T & x0, const T & y0, const T & alpha, const T & y, T & x, T & dx_dy )
  {
    const T alpha_inv = T( 1.0 ) / alpha;
    x = x0 + alpha_inv * log( y / y0 );
    dx_dy = alpha_inv / y;
  }
};

template<typename T>
struct ExponentialCompute<T, ExponentApproximationType::Quadratic, false>
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
};

template<typename T>
struct ExponentialCompute<T, ExponentApproximationType::Quadratic, true>
{
  inline static void Compute( const T & x0, const T & y0, const T & alpha, const T & y, T & x )
  {
    const T z = sqrt( T( 2.0 ) * y / y0 - T( 1.0 ));
    x = x0 + (z - 1) / alpha;
  }

  inline static void Compute( const T & x0, const T & y0, const T & alpha, const T & y, T & x, T & dx_dy )
  {
    const T alpha_inv = T( 1.0 ) / alpha;
    const T z = sqrt( T( 2.0 ) * y / y0 - T( 1.0 ));
    x = x0 + alpha_inv * (z - 1);
    dx_dy = alpha_inv / (y0 * z);
  }
};

template<typename T>
struct ExponentialCompute<T, ExponentApproximationType::Linear, false>
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
};

template<typename T>
struct ExponentialCompute<T, ExponentApproximationType::Linear, true>
{
  inline static void Compute( const T & x0, const T & y0, const T & alpha, const T & y, T & x )
  {
    x = x0 + (y / y0 - 1) / alpha;
  }

  inline static void Compute( const T & x0, const T & y0, const T & alpha, const T & y, T & x, T & dx_dy )
  {
    const T alpha_inv = T( 1.0 ) / alpha;
    x = x0 + alpha_inv * (y / y0 - 1);
    dx_dy = alpha_inv / y0;
  }
};

}

template<typename T, ExponentApproximationType EAT>
ExponentialRelation<T, EAT>::ExponentialRelation()
  : ExponentialRelation(T(0), T(1), T(1))
{

}

template<typename T, ExponentApproximationType EAT>
ExponentialRelation<T, EAT>::ExponentialRelation( T x0, T y0, T alpha )
{
  SetCoefficients( x0, y0, alpha );
}

template<typename T, ExponentApproximationType EAT>
void ExponentialRelation<T, EAT>::SetCoefficients( T x0, T y0, T alpha )
{
  m_x0 = x0;
  m_y0 = y0;
  m_alpha = alpha;
}

template<typename T, ExponentApproximationType EAT>
void ExponentialRelation<T, EAT>::Compute( const T & x, T & y ) const
{
  detail::ExponentialCompute<T, EAT>::Compute( m_x0, m_y0, m_alpha, x, y );
}

template<typename T, ExponentApproximationType EAT>
void ExponentialRelation<T, EAT>::Inverse( const T & y, T & x ) const
{
  detail::ExponentialCompute<T, EAT, true>::Compute( m_x0, m_y0, m_alpha, y, x );
}

template<typename T, ExponentApproximationType EAT>
void ExponentialRelation<T, EAT>::Compute( const T & x, T & y, T & dy_dx ) const
{
  detail::ExponentialCompute<T, EAT>::Compute( m_x0, m_y0, m_alpha, x, y, dy_dx );
}

template<typename T, ExponentApproximationType EAT>
void ExponentialRelation<T, EAT>::Inverse( const T & y, T & x, T & dx_dy ) const
{
  detail::ExponentialCompute<T, EAT, true>::Compute( m_x0, m_y0, m_alpha, y, x, dx_dy );
}

} // namespace constitutive

} // namespace geosx

#endif // GEOSX_EXPONENTIALRELATION_HPP_
