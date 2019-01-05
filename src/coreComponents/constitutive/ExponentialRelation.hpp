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
 * @file ExponentialRelations.hpp
 */

#ifndef GEOSX_EXPONENTIALRELATION_HPP_
#define GEOSX_EXPONENTIALRELATION_HPP_

#include <common/DataTypes.hpp>

namespace geosx
{

namespace constitutive
{

/**
 * @enum enumeration describing available approximations of the exponent
 */
enum class ExponentApproximationType
{
  Full,
  Linear,
  Quadratic
};

template< ExponentApproximationType EAT >
struct ExponentApproximationTypeWrapper
{
  static constexpr ExponentApproximationType value = EAT;
};

template< typename LAMBDA >
static auto ExponentApproximationTypeSwitchBlock( ExponentApproximationType const eat,
                                                  LAMBDA&& lambda );

/**
 * @class ExponentialRelation
 *
 * Describes a simple exponential relationship with three parameters of the form: \f$ y = y0 * exp(a * (x - x0)) \f$.
 * Optional linear and quadratic Taylor series approximation of the exponent.
 * Implements point and batch direct (y(x)) and inverse (x(y)) compute functions and provides derivatives.
 *
 * @tparam T scalar real type used in computation
 */
template<typename T>
class ExponentialRelation
{
public:

  // *** constructors ***

  explicit ExponentialRelation();

  explicit ExponentialRelation( ExponentApproximationType type );

  explicit ExponentialRelation( ExponentApproximationType type, T x0, T y0, T alpha );

  // *** setters ***

  void SetApproximationType( ExponentApproximationType type );

  void SetCoefficients( T x0, T y0, T alpha );

  void SetParameters( ExponentApproximationType type, T x0, T y0, T alpha );

  // *** no-derivative computes ***

  void Compute( const T & x, T & y ) const;

  void Inverse( const T & y, T & x ) const;

  template<typename Policy>
  void Compute( arrayView1d<T const> const & x, arrayView1d<T> & y ) const;

  template<typename Policy>
  void Inverse( arrayView1d<T const> const & y, arrayView1d<T> & x ) const;

  template<typename Policy>
  void Compute( arrayView1d<T const> const & x, arrayView2d<T> & y ) const;

  template<typename Policy>
  void Inverse( arrayView1d<T const> const & y, arrayView2d<T> & x ) const;

  template<typename Policy>
  void Compute( arrayView2d<T const> const & x, arrayView2d<T> & y ) const;

  template<typename Policy>
  void Inverse( arrayView2d<T const> const & y, arrayView2d<T> & x ) const;

  // *** derivative computes ***

  void Compute( const T & x, T & y, T & dy_dx ) const;

  void Inverse( const T & y, T & x, T & dx_dy ) const;

  template<typename Policy>
  void Compute( arrayView1d<T const> const & x, arrayView1d<T> & y, arrayView1d<T> & dy_dx ) const;

  template<typename Policy>
  void Inverse( arrayView1d<T const> const & y, arrayView1d<T> & x, arrayView1d<T> & dx_dy ) const;

  template<typename Policy>
  void Compute( arrayView1d<T const> const & x, arrayView2d<T> & y, arrayView2d<T> & dy_dx ) const;

  template<typename Policy>
  void Inverse( arrayView1d<T const> const & y, arrayView2d<T> & x, arrayView2d<T> & dx_dy ) const;

  template<typename Policy>
  void Compute( arrayView2d<T const> const & x, arrayView2d<T> & y, arrayView2d<T> & dy_dx ) const;

  template<typename Policy>
  void Inverse( arrayView2d<T const> const & y, arrayView2d<T> & x, arrayView2d<T> & dx_dy ) const;

private:

  /// indicator of exponent approximation to use
  ExponentApproximationType m_approximationType;

  /// coefficients of exponential relations
  T m_x0;
  T m_y0;
  T m_alpha;

};

} // namespace constitutive

} // namespace geosx

#include "constitutive/ExponentialRelation_impl.hpp"

#endif // GEOSX_EXPONENTIALRELATION_HPP_
