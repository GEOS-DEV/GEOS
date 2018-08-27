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
 * @tparam IndexType index type used in batch updates
 * @tparam RealType scalar real type used in computation
 */
template<typename IndexType, typename RealType>
class ExponentialRelation
{
public:

  using PtrType = RealType *;
  using ConstPtrType = RealType const *;

  // *** constructors ***

  explicit ExponentialRelation();

  explicit ExponentialRelation( ExponentApproximationType type );

  explicit ExponentialRelation( ExponentApproximationType type, RealType x0, RealType y0, RealType alpha );

  // *** setters ***

  void SetApproximationType( ExponentApproximationType type );

  void SetCoefficients( RealType x0, RealType y0, RealType alpha );

  void SetParameters( ExponentApproximationType type, RealType x0, RealType y0, RealType alpha );

  // *** no-derivative computes ***

  void Compute( const RealType & x, RealType & y );

  void Inverse( const RealType & y, RealType & x );

  template<typename Policy>
  void BatchCompute( IndexType size, ConstPtrType x_ptr, PtrType y_ptr );

  template<typename Policy>
  void BatchInverse( IndexType size, ConstPtrType y_ptr, PtrType x_ptr );

  // *** derivative computes ***

  void Compute( const RealType & x, RealType & y, RealType & dy_dx );

  void Inverse( const RealType & y, RealType & x, RealType & dx_dy );

  template<typename Policy>
  void BatchCompute( IndexType size, ConstPtrType x_ptr, PtrType y_ptr, PtrType dy_dx_ptr );

  template<typename Policy>
  void BatchInverse( IndexType size, ConstPtrType y_ptr, PtrType x_ptr, PtrType dx_dy_ptr );

private:

  /// indicator of exponent approximation to use
  ExponentApproximationType m_approximationType;

  /// coefficients of exponential relations
  RealType m_x0;
  RealType m_y0;
  RealType m_alpha;

};

} // namespace constitutive

} // namespace geosx

#include "constitutive/ExponentialRelation_impl.hpp"

#endif // GEOSX_EXPONENTIALRELATION_HPP_
