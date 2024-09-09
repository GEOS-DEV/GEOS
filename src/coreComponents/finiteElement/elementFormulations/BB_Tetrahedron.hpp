/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file BB_Tetrahedron.hpp
 */

#ifndef GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_BBTETRAHEDRON_HPP_
#define GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_BBTETRAHEDRON_HPP_

#include "FiniteElementBase.hpp"
#include <utility>



namespace geos
{
namespace finiteElement
{

/**
 * This class contains the kernel accessible functions specific to the
 * Bernstein-Bézier (BB) modal any-order tetrahedron finite element with 
 * Gaussian quadrature rules. Available functions are tailored for
 * Discontinuous Galerkin (DG) applications.
 */
template< int k >
class BB_Tetrahedron final : public FiniteElementBase
{
public:

  /// The number of shape functions per element.
  constexpr static localIndex numNodes = ( k + 1 ) * ( k + 2 ) * ( k + 3 ) / 6;

  /// The number of shape functions per face
  constexpr static localIndex numNodesPerFace = ( k + 1 ) * ( k + 2 ) / 2;

  /// The maximum number of support points per element.
  constexpr static localIndex maxSupportPoints = numNodes;

  /// The number of quadrature points per element.
  constexpr static localIndex numQuadraturePoints = numNodes;

  /** @cond Doxygen_Suppress */
  USING_FINITEELEMENTBASE
  /** @endcond Doxygen_Suppress */

  virtual ~BB_Tetrahedron() override
  {}

  GEOS_HOST_DEVICE
  virtual localIndex getNumQuadraturePoints() const override
  {
    return numQuadraturePoints;
  }

  /**
   * @brief Get the number of quadrature points.
   * @param stack Stack variables as filled by @ref setupStack.
   * @return The number of quadrature points.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static localIndex getNumQuadraturePoints( StackVariables const & stack )
  {
    GEOS_UNUSED_VAR( stack );
    return numQuadraturePoints;
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual localIndex getNumSupportPoints() const override
  {
    return numNodes;
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual localIndex getMaxSupportPoints() const override
  {
    return maxSupportPoints;
  }

  /**
   * @brief Get the number of support points.
   * @param stack Object that holds stack variables.
   * @return The number of support points.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static localIndex getNumSupportPoints( StackVariables const & stack )
  {
    GEOS_UNUSED_VAR( stack );
    return numNodes;
  }


  /**
   * @brief Calculate shape functions values at a single point in the reference element using De Casteljau's algorithm.
   * @param[in] coords The parent coordinates at which to evaluate the shape function value, in the reference element
   * @param[out] N The shape function values.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void calcN( real64 const (&coords)[3],
                     real64 (& N)[numNodes] )
  {
    return calcN( { 1.0 - coords[ 0 ] - coords[ 1 ] - coords[ 2 ],
                    coords[ 0 ], coords[ 1 ], coords[ 2 ] }, N );
  }

  /**
   * @brief Calculate shape functions values at a single point, given the coorginates of the tetrahedron vertices, using De Casteljau's algorithm.
   * @param[in] coords The parent coordinates at which to evaluate the shape function value, in the reference element
   * @param[out] N The shape function values.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void calcN( real64 const (&X)[4][3],
                     real64 const (&coords)[3],
                     real64 (& N)[numNodes] )
  {
    real64 lambda[4] = {};
    real64 m[3][3] = {};
    for( int i = 0; i < 3; i++ )
    {
      for( int j = 0; j < 3; j++ )
      {
        m[ i ][ j ] = X[ i + 1 ][ j ] - X[ 0 ][ j ];
      }
    }
    real64 den = LvARray::math::abs( LvArray::tensorOps::determinant< 3 >( m ) );
    for( int i = 0; i < 3; i++ )
    {
      for( int j = 0; j < 3; j++ )
      {
        m[ i ][ j ] = coords[ j ] - X[ 0 ][ j ];
      }
      lambda[ i + 1 ] = LvArray::math::abs( LvArray::tensorOps::determinant< 3 >( m ) ) / den;
      for( int j = 0; j < 3; j++ )
      {
        m[ i ][ j ] = X[ i + 1 ][ j ] - X[ 0 ][ j ];
      }
    }
    lambda[ 0 ] = 1.0 - lambda[ 1 ] - lambda[ 2 ] - lambda[ 3 ];
    return calcN( lambda, N );
  }

  /**
   * @brief Calculate shape functions values at a single point using De Casteljau's algorithm.
   * @param[in] lambda barycentric coordinates of the point in thetetrahedron
   * @param[out] N The shape function values.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void calcN( real64 const ( & lambda)[4],
                     real64 (& N)[numNodes] )
  {
    N[ 0 ] = 6.0;
    int prev;
    int c;
    int limits[ 4 ] = { 1, 1, 1, 1 };
    for( int kp = 1; kp <= k; kp++)
    {
      prev = kp * ( kp + 1 ) * ( kp + 2 ) / 6 - 1; 
      c = ( kp + 1 ) * ( kp + 2 ) * ( kp + 3 ) / 6 - 1; 
      for( int i = 0; i < 4; i++ )
      {
        int denominator = i == 0 ? kp : 1;
        int offset = 0;
        int c1 = kp - 1;
        int c2 = i + kp - 2;
        int repetitionCount = i == 0 ? 1 : limits[ i - 1 ];
        for( int j = 0; j < limits[ i ] ; j++ )
        {
          if( j - offset >=  repetitionCount )
          {
            denominator++;
            offset += repetitionCount;
            repetitionCount = repetitionCount * c1 / c2;
            c1--;
            c2--;
          }
          N[ c-- ] = N[ prev - j ] * lambda[3 - i] * ( kp + 3 ) / denominator;
        }
      }
      for( int i = 1; i < 4; i++ )
      {
        limits[ i ] += limits[ i - 1 ];
      }
    } 
  }
 
  /**
   * @brief Calculate shape functions values at a single point in the reference element using De Casteljau's algorithm.
   * @param[in] coords The parent coordinates at which to evaluate the shape function value, in the reference element
   * @param[out] N The shape function values.
   * @param[out] gradN The derivatives of the shape function values.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void calcNandGradN( real64 const (&coords)[3],
                             real64 (& N)[numNodes] )
                             real64 (& gradN)[numNodes][3] )
  {
    real64 dNdLambda[numNodes][4];
    calcNandGradN( { 1.0 - coords[ 0 ] - coords[ 1 ] - coords[ 2 ],
                     coords[ 0 ], coords[ 1 ], coords[ 2 ] }, N, dNdLambda );
    for( int i = 0; i < numNodes; i++ )
    {
      for( int j = 0; j < 3; j++ )
      {
        gradN[ i ][ j ] = dNdLambda[ i ][ j + 1 ] - dNdLambda[ i ][ 0 ]
      }
    }
  }

  /**
   * @brief Calculate shape functions values at a single point, given the coorginates of the tetrahedron vertices, using De Casteljau's algorithm.
   * @param[in] coords The parent coordinates at which to evaluate the shape function value, in the reference element
   * @param[out] N The shape function values.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void calcN( real64 const (&X)[4][3],
                     real64 const (&coords)[3],
                     real64 (& N)[numNodes] )
  {
    real64 lambda[4] = {};
    real64 m[3][3] = {};
    real64 dNdLambda[numNodes][4];
    real64 dLambdadX[4][4];
    for( int i = 0; i < 3; i++ )
    {
      for( int j = 0; j < 3; j++ )
      {
        m[ i ][ j ] = X[ i + 1 ][ j ] - X[ 0 ][ j ];
      }
    }
    real64 den = LvARray::math::abs( LvArray::tensorOps::determinant< 3 >( m ) );
    for( int i = 0; i < 3; i++ )
    {
      for( int j = 0; j < 3; j++ )
      {
        m[ i ][ j ] = coords[ j ] - X[ 0 ][ j ];
      }
      lambda[ i + 1 ] = LvArray::math::abs( LvArray::tensorOps::determinant< 3 >( m ) ) / den;
      for( int j = 0; j < 3; j++ )
      {
        m[ i ][ j ] = X[ i + 1 ][ j ] - X[ 0 ][ j ];
      }
    }
    lambda[ 0 ] = 1.0 - lambda[ 1 ] - lambda[ 2 ] - lambda[ 3 ];
    calcNandGradN( lambda, N, dNdLambda );
    for( int i = 0; i < numNodes; i++ )
    {
      for( int j = 0; j < 3; j++ )
      {
        gradN[ i ][ j ] = ( ( ( X[ 2 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ]) * ( X[ 3 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) - ( X[ 3 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ]) * ( X[ 2 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) )* ( dNdLambda[ i ][ 1 ] - dNdLambda[ i ][ 0 ] ) +
                            ( ( X[ 3 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ]) * ( X[ 1 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) - ( X[ 1 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ]) * ( X[ 3 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) )* ( dNdLambda[ i ][ 2 ] - dNdLambda[ i ][ 0 ] ) +
                            ( ( X[ 1 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ]) * ( X[ 2 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) - ( X[ 2 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ]) * ( X[ 1 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) )* ( dNdLambda[ i ][ 3 ] - dNdLambda[ i ][ 0 ] )) / den 
      }
    }
  }


  /**
   * @brief Calculate the derivatives of shape functions with respect to barycentric coordinates at a single point using De Casteljau's algorithm.
   * @param[in] lambda barycentric coordinates of the point in thetetrahedron
   * @param[out] N The shape function values.
   * @param[out] gradN The derivatives of the shape functions with respect to the lambdas
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void calcNandGradN( real64 const ( & lambda)[4],
                             real64 const ( & N)[numNodes],
                             real64 (& gradN)[numNodes][ 4 ] )
  {                  
    gradN[ 0 ][ 0 ] = 0.0;
    gradN[ 0 ][ 1 ] = 0.0;
    gradN[ 0 ][ 2 ] = 0.0;
    gradN[ 0 ][ 3 ] = 0.0;
    N[ 0 ] = 6.0;
    int prev;
    int c;
    int limits[ 4 ] = { 1, 1, 1, 1 };
    for( int kp = 1; kp <= k; kp++)
    {
      prev = kp * ( kp + 1 ) * ( kp + 2 ) / 6 - 1; 
      c = ( kp + 1 ) * ( kp + 2 ) * ( kp + 3 ) / 6 - 1; 
      for( int i = 0; i < 4; i++ )
      {
        int denominator = i == 0 ? kp : 1;
        int offset = 0;
        int c1 = kp - 1;
        int c2 = i + kp - 2;
        int repetitionCount = i == 0 ? 1 : limits[ i - 1 ];
        for( int j = 0; j < limits[ i ] ; j++ )
        {
          if( j - offset >=  repetitionCount )
          {
            denominator++;
            offset += repetitionCount;
            repetitionCount = repetitionCount * c1 / c2;
            c1--;
            c2--;
          }
          gradN[ c ][ 0 ] = gradN[ prev - j ][ 0 ] * lambda[ 3 - i ] * ( kp + 3 ) / denominator;
          gradN[ c ][ 1 ] = gradN[ prev - j ][ 1 ] * lambda[ 3 - i ] * ( kp + 3 ) / denominator;
          gradN[ c ][ 2 ] = gradN[ prev - j ][ 2 ] * lambda[ 3 - i ] * ( kp + 3 ) / denominator;
          gradN[ c ][ 3 ] = gradN[ prev - j ][ 3 ] * lambda[ 3 - i ] * ( kp + 3 ) / denominator;
          gradN[ c ][ 3 - i ] += N[ prev - j ] * ( kp + 3 ) / denominator;
          N[ c-- ] = N[ prev - j ] * lambda[ 3 - i ] * ( kp + 3 ) / denominator;
        }
      }
      for( int i = 1; i < 4; i++ )
      {
        limits[ i ] += limits[ i - 1 ];
      }
    } 
  }

  /**
   * @brief Calculate shape functions values at a single point in the reference element using De Casteljau's algorithm.
   * @param[in] face The index of the vertex opposite to the dersired face
   * @param[in] coords The parent coordinates at which to evaluate the shape function value, in the reference element
   * @param[out] N The shape function values.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void calcN( real64 const (&coords)[3],
                     real64 (& N)[numNodes] )
  {
    return calcN( { 1.0 - coords[ 0 ] - coords[ 1 ] - coords[ 2 ],
                    coords[ 0 ], coords[ 1 ], coords[ 2 ] }, N );
  }

  /**
   * @brief Calculate shape functions values at a single point, given the coorginates of the tetrahedron vertices, using De Casteljau's algorithm.
   * @param[in] coords The parent coordinates at which to evaluate the shape function value, in the reference element
   * @param[out] N The shape function values.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void calcN( real64 const (&X)[4][3],
                     real64 const (&coords)[3],
                     real64 (& N)[numNodes] )
  {
    real64 lambda[4] = {};
    real64 m[3][3] = {};
    for( int i = 0; i < 3; i++ )
    {
      for( int j = 0; j < 3; j++ )
      {
        m[ i ][ j ] = X[ i + 1 ][ j ] - X[ 0 ][ j ];
      }
    }
    real64 den = LvARray::math::abs( LvArray::tensorOps::determinant< 3 >( m ) );
    for( int i = 0; i < 3; i++ )
    {
      for( int j = 0; j < 3; j++ )
      {
        m[ i ][ j ] = coords[ j ] - X[ 0 ][ j ];
      }
      lambda[ i + 1 ] = LvArray::math::abs( LvArray::tensorOps::determinant< 3 >( m ) ) / den;
      for( int j = 0; j < 3; j++ )
      {
        m[ i ][ j ] = X[ i + 1 ][ j ] - X[ 0 ][ j ];
      }
    }
    lambda[ 0 ] = 1.0 - lambda[ 1 ] - lambda[ 2 ] - lambda[ 3 ];
    return calcN( lambda, N );
  }
  /**
   * @brief Calculate shape functions values at a single point on a face of the reference element using De Casteljau's algorithm.
   * @param[in] face The index of the vertex opposite to the desired face
   * @param[in] coords The parent coordinates at which to evaluate the shape function value, in the reference element
   * @param[out] N The shape function values.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void calcN( localIndex face,
                     real64 const (&coords)[3],
                     real64 (& N)[numNodesPerFace] )
  {
    switch( face )
    {
      case 0: 
        return calcN( { coords[ 0 ], coords[ 1 ], coords[ 2 ] }, N );
      case 1: 
        return calcN( { 1.0 - coords[ 0 ] - coords[ 1 ] - coords[ 2 ], coords[ 1 ], coords[ 2 ] }, N );
      case 2: 
        return calcN( { 1.0 - coords[ 0 ] - coords[ 1 ] - coords[ 2 ], coords[ 0 ], coords[ 2 ] }, N );
      case 3: 
        return calcN( { 1.0 - coords[ 0 ] - coords[ 1 ] - coords[ 2 ], coords[ 0 ], coords[ 1 ] }, N );
      default:
        return 0;
    }
  }

  /**
   * @brief Calculate shape functions values at a single point on a face using De Casteljau's algorithm.
   * @param[in] lambda barycentric coordinates of the point in triangle
   * @param[out] N The shape function values.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void calcN( real64 const ( & lambda)[3],
                     real64 (& N)[numNodesPerFace] )
  {
    N[ 0 ] = 2.0;
    int prev;
    int c;
    int limits[ 3 ] = { 1, 1, 1 };
    for( int kp = 1; kp <= k; kp++)
    {
      prev = kp * ( kp + 1 ) / 2 - 1; 
      c = ( kp + 1 ) * ( kp + 2 ) / 2 - 1; 
      for( int i = 0; i < 3; i++ )
      {
        int denominator = i == 0 ? kp : 1;
        int offset = 0;
        int c1 = kp - 1;
        int c2 = i + kp - 2;
        int repetitionCount = i == 0 ? 1 : limits[ i - 1 ];
        for( int j = 0; j < limits[ i ] ; j++ )
        {
          if( j - offset >=  repetitionCount )
          {
            denominator++;
            offset += repetitionCount;
            repetitionCount = repetitionCount * c1 / c2;
            c1--;
            c2--;
          }
          N[ c-- ] = N[ prev - j ] * lambda[2 - i] * ( kp + 2 ) / denominator;
        }
      }
      for( int i = 1; i < 3; i++ )
      {
        limits[ i ] += limits[ i - 1 ];
      }
    } 
  }

  /**
   * @brief Calculate the derivatives of shape functions with respect to barycentric coordinates at a single point on a face using De Casteljau's algorithm.
   * @param[in] lambda barycentric coordinates of the point in the triangle
   * @param[out] N The shape function values.
   * @param[out] gradN The derivatives of the shape functions with respect to the lambdas
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void calcNandGradN( real64 const ( & lambda)[ 3 ],
                             real64 const ( & N)[ numNodes ],
                             real64 (& gradN)[ numNodes ][ 3 ] )
  {                  
    gradN[ 0 ][ 0 ] = 0.0;
    gradN[ 0 ][ 1 ] = 0.0;
    gradN[ 0 ][ 2 ] = 0.0;
    N[ 0 ] = 2.0;
    int prev;
    int c;
    int limits[ 3 ] = { 1, 1, 1 };
    for( int kp = 1; kp <= k; kp++)
    {
      prev = kp * ( kp + 1 ) / 2 - 1; 
      c = ( kp + 1 ) * ( kp + 2 ) / 2 - 1; 
      for( int i = 0; i < 3; i++ )
      {
        int denominator = i == 0 ? kp : 1;
        int offset = 0;
        int c1 = kp - 1;
        int c2 = i + kp - 2;
        int repetitionCount = i == 0 ? 1 : limits[ i - 1 ];
        for( int j = 0; j < limits[ i ] ; j++ )
        {
          if( j - offset >=  repetitionCount )
          {
            denominator++;
            offset += repetitionCount;
            repetitionCount = repetitionCount * c1 / c2;
            c1--;
            c2--;
          }
          gradN[ c ][ 0 ] = gradN[ prev - j ][ 0 ] * lambda[ 2 - i ] * ( kp + 2 )/ denominator;
          gradN[ c ][ 1 ] = gradN[ prev - j ][ 1 ] * lambda[ 2 - i ] * ( kp + 2 )/ denominator;
          gradN[ c ][ 2 ] = gradN[ prev - j ][ 2 ] * lambda[ 2 - i ] * ( kp + 2 )/ denominator;
          gradN[ c ][ 2 - i ] += N[ prev - j ] * ( kp + 2 ) / denominator;
          N[ c-- ] = N[ prev - j ] * lambda[ 2 - i ] * ( kp + 2 ) / denominator;
        }
      }
      for( int i = 1; i < 3; i++ )
      {
        limits[ i ] += limits[ i - 1 ];
      }
    } 
  }

  /**
   * @brief Calculate the shape functions derivatives wrt the physical
   *   coordinates at a single point.
   * @param[in] coords The parent coordinates at which to evaluate the shape function value
   * @param[in] X Array containing the coordinates of the support points.
   * @param[out] gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   * @return The determinant of the parent/physical transformation matrix.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 calcGradN( real64 const (&coords)[3],
                           real64 const (&X)[numNodes][3],
                           real64 ( &gradN )[numNodes][3] );

  /**
   * @brief Calculate the shape functions derivatives wrt the physical
   *   coordinates.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the mesh support points.
   * @param stack Variables allocated on the stack as filled by @ref setupStack.
   * @param gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   * @return The determinant of the parent/physical transformation matrix.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 calcGradN( localIndex const q,
                           real64 const (&X)[numNodes][3],
                           StackVariables const & stack,
                           real64 ( &gradN )[numNodes][3] );

  /**
   * @brief Calculate the integration weights for a quadrature point.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the support points.
   * @return The product of the quadrature rule weight and the determinate of
   *   the parent/physical transformation matrix.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 transformedQuadratureWeight( localIndex const q,
                                             real64 const (&X)[numNodes][3] );

  /**
   * @brief Calculates the isoparametric "Jacobian" transformation
   *   matrix/mapping from the parent space to the physical space on a 2D domain (face).
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param X Array containing the coordinates of the mesh support points.
   * @param J Array to store the Jacobian transformation.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void jacobianTransformation2d( int const qa,
                                        int const qb,
                                        real64 const (&X)[4][3],
                                        real64 ( &J )[3][2] );


  /**
   * @brief Calculates the isoparametric "Jacobian" transformation
   *   matrix/mapping from the parent space to the physical space.
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param qc The 1d quadrature point index in xi2 direction (0,1)
   * @param X Array containing the coordinates of the mesh support points.
   * @param J Array to store the Jacobian transformation.
   * @return The determinant of the Jacobian transformation matrix.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 invJacobianTransformation( int const qa,
                                           int const qb,
                                           int const qc,
                                           real64 const (&X)[8][3],
                                           real64 ( & J )[3][3] )
  {
    jacobianTransformation( qa, qb, qc, X, J );
    return LvArray::tensorOps::invert< 3 >( J );
  }

  /**
   * @brief Calculates the isoparametric "Jacobian" transformation
   *   matrix/mapping from the parent space to the physical space.
   * @param q The quadrature point index
   * @param X Array containing the coordinates of the mesh support points.
   * @param J Array to store the Jacobian transformation.
   * @return The determinant of the Jacobian transformation matrix.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 invJacobianTransformation( int const q,
                                           real64 const (&X)[8][3],
                                           real64 ( & J )[3][3] )
  {
    int qa, qb, qc;
    GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );
    return invJacobianTransformation( qa, qb, qc, X, J );
  }


  /**
   * @brief Calculate the symmetric gradient of a vector valued support field
   *   at a quadrature point using the stored inverse of the Jacobian
   *   transformation matrix.
   * @param q The quadrature point index
   * @param invJ The inverse of the Jacobian transformation matrix.
   * @param var The vector valued support field to apply the gradient
   *   operator on.
   * @param grad The symmetric gradient in Voigt notation.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void symmetricGradient( int const q,
                                 real64 const (&invJ)[3][3],
                                 real64 const (&var)[numNodes][3],
                                 real64 ( &grad )[6] )
{
  supportLoop( q, [] GEOS_HOST_DEVICE ( real64 const (&dNdXi)[3],
                                        int const nodeIndex,
                                        real64 const (&invJ)[3][3],
                                        real64 const (&var)[numNodes][3],
                                        real64 (& grad)[6] )
  {

    real64 gradN[3] = {0, 0, 0};
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        gradN[i] = gradN[i] + dNdXi[ j ] * invJ[j][i];
      }
    }

    grad[0] = grad[0] + gradN[0] * var[ nodeIndex ][0];
    grad[1] = grad[1] + gradN[1] * var[ nodeIndex ][1];
    grad[2] = grad[2] + gradN[2] * var[ nodeIndex ][2];
    grad[3] = grad[3] + gradN[2] * var[ nodeIndex ][1] + gradN[1] * var[ nodeIndex ][2];
    grad[4] = grad[4] + gradN[2] * var[ nodeIndex ][0] + gradN[0] * var[ nodeIndex ][2];
    grad[5] = grad[5] + gradN[1] * var[ nodeIndex ][0] + gradN[0] * var[ nodeIndex ][1];
  }, invJ, var, grad );
}


  /**
   * @brief Calculate the gradient of a vector valued support field at a point
   *   using the stored basis function gradients for all support points.
   * @param q The quadrature point index
   * @param invJ The inverse of the Jacobian transformation matrix.
   * @param var The vector valued support field to apply the gradient
   *   operator on.
   * @param grad The gradient.
   *
   * More precisely, the operator is defined as:
   * \f[
   * grad_{ij}  = \sum_a^{nSupport} \left ( \frac{\partial N_a}{\partial X_j} var_{ai}\right ),
   * \f]
   *
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void gradient( int const q,
                        real64 const (&invJ)[3][3],
                        real64 const (&var)[numNodes][3],
                        real64 ( &grad )[3][3] );


  /**
   * @brief Inner product of all basis function gradients and a rank-2
   *   symmetric tensor evaluated at a quadrature point.
   * @param q The 3d quadrature point index
   * @param invJ The inverse of the Jacobian transformation matrix.
   * @param var The rank-2 symmetric tensor at @p q.
   * @param R The vector resulting from the tensor contraction.
   *
   * More precisely, the operator is defined as:
   * \f[
   * R_i = \sum_a^{nSupport} \left( \frac{\partial N_a}{\partial X_j} var_{ij} \right),
   * \f]
   * where \f$\frac{\partial N_a}{\partial X_j}\f$ is the basis function gradient,
   *   \f$var_{ij}\f$ is the rank-2 symmetric tensor.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void plusGradNajAij( int const q,
                              real64 const (&invJ)[3][3],
                              real64 const (&var)[6],
                              real64 ( &R )[numNodes][3] );



  /**
   * @brief Calculates the isoparametric "Jacobian" transformation
   *   matrix/mapping from the parent space to the physical space.
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param qc The 1d quadrature point index in xi2 direction (0,1)
   * @param X Array containing the coordinates of the mesh support points.
   * @param J Array to store the Jacobian transformation.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void jacobianTransformation( int const qa,
                                      int const qb,
                                      int const qc,
                                      real64 const (&X)[8][3],
                                      real64 ( &J )[3][3] );

  /**
   * @brief Calculates the isoparametric "Jacobian" transformation
   *   matrix/mapping from the parent space to the physical space at a single point.
   * @param coords The parent coordinates at which to evaluate the shape function value
   * @param X Array containing the coordinates of the support points.
   * @param J Array to store the Jacobian transformation.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void jacobianTransformation( real64 const (&coords)[3],
                                      real64 const (&X)[numNodes][3],
                                      real64 ( &J )[3][3] );
  /**
   * @brief performs a trilinear interpolation to determine the real-world coordinates of a
   *   vertex
   * @param[in] alpha Interpolation coefficient in [0,1] for the first coordinate
   * @param[in] beta Interpolation coefficient in [0,1] for the second coordinate
   * @param[in] gamma Interpolation coefficient in [0,1] for the third coordinate
   * @param[in] X Real-world coordinates of the cell corners
   * @param[out] coords Real-world coordinates of the interpolated point
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void trilinearInterp( real64 const alpha,
                               real64 const beta,
                               real64 const gamma,
                               real64 const (&X)[8][3],
                               real64 ( &coords )[3] );

  /**
   * @brief computes the real-world coordinates of the support nodes
   * @param[in] Xmesh Array containing the coordinates of the corners of the mesh element
   * @param[out] X Array containing the coordinates of the support points.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void computeLocalCoords( real64 const (&Xmesh)[8][3],
                                  real64 const (&X)[numNodes][3] );


  /**
   * @brief Computes the superposition integral between Bernstein-Bézier functions indexed by 
   *  (i1, j1, k1, l1) and (i1, j2, k2, l2)
   * @param[in] i1
   * @param[in] j1
   * @param[in] k1
   * @param[in] l1
   * @param[in] i2
   * @param[in] j2
   * @param[in] k2
   * @param[in] l2
   * @return the superposition integral over the barycentric coordinates
   */
  
   constexpr static real64 computeSuperpositionIntegral(
       const int i1, const int j1, const int k1, const int l1, 
       const int i2, const int j2, const int k2, const int l2 )
   {
     std::cout << "test: " << integralTerm(i1+i2, i1, i2) << " " << integralTerm(j1+j2, j1, j2) << " " << integralTerm(k1+k2, k1, k2) << " " <<  integralTerm(l1+l2, l1, l2) << " " << integralTerm(i1+j1+k1+l1+i2+j2+k2+l2+3, i1+j1+k1+l1+3, i2+j2+k2+l2+3) << std::endl;
     return (integralTerm(i1+i2, i1, i2)*
             integralTerm(j1+j2, j1, j2)*
             integralTerm(k1+k2, k1, k2)*
             integralTerm(l1+l2, l1, l2))/
             integralTerm(i1+j1+k1+l1+i2+j2+k2+l2+3, 
                          i1+j1+k1+l1+3, i2+j2+k2+l2+3);
   }
  /**
   * @brief Computes a! / ( b! * c! ) with b + c >= a >= b, c
   * @param[in] a
   * @param[in] b
   * @param[in] c
   * @return a!/(b!*c!)
   */
   constexpr static real64 integralTerm(const int a, const int b, const int c)
   {
     real64 res = 1.0;
     int num = a;
     int den = c;
     for( int i = b; i > 0; i--)
     {
         res *= ( (real64) num ) /  i;
         num--;
     }
     for( int i = num; i > 0; i--)
     {
         res *= ( (real64) i ) /  den;
         den--;
     }
     for( int i = den; i > 0; i--)
     {
         res /= i;
     }
     
     return res;
  }

  /**
   * @brief computes the non-zero contributions of the d.o.f. indexd by q to the
   *   mass matrix M, i.e., the superposition matrix of the shape functions.
   * @param q The quadrature point index
   * @param X Array containing the coordinates of the mesh support points.
   * @return The diagonal mass term associated to q
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 computeMassTerm( localIndex const q,
                                 real64 const (&X)[4][3] )
  {
    real64 m[3][3] = {};
    for( int i = 0; i < 3; i++ )
    {
      for( int j = 0; j < 3; j++ )
      {
        m[ i ][ j ] = X[ i + 1 ][ j ] - X[ 0 ][ j ];
      }
    }
    real64 detJ = LvARray::math::abs( LvArray::tensorOps::determinant< 3 >( m ) );
    real64 J[3][3] = {{0}};
    jacobianTransformation( qa, qb, qc, X, J );
    return LvArray::math::abs( LvArray::tensorOps::determinant< 3 >( J ) )*w3D;
  }

  /**
   * @brief computes the non-zero contributions of the d.o.f. indexd by q to the
   *   damping matrix M, i.e., the superposition matrix of the shape functions
   *   integrated over a face.
   * @param q The quadrature point index
   * @param X Array containing the coordinates of the support points.
   * @return The diagonal damping term associated to q
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 computeDampingTerm( localIndex const q,
                                    real64 const (&X)[4][3] );

  /**
   * @brief computes the matrix B, defined as J^{-T}J^{-1}/det(J), where J is the Jacobian matrix,
   *   at the given Gauss-Lobatto point.
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param qc The 1d quadrature point index in xi2 direction (0,1)
   * @param X Array containing the coordinates of the support points.
   * @param J Array to store the Jacobian
   * @param B Array to store the matrix B, in Voigt notation
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void computeBMatrix( int const qa,
                              int const qb,
                              int const qc,
                              real64 const (&X)[8][3],
                              real64 ( &J )[3][3],
                              real64 ( &B )[6] );

  /**
   * @brief computes the non-zero contributions of the d.o.f. indexed by q to the
   *   stiffness matrix R, i.e., the superposition matrix of first derivatives
   *   of the shape functions.
   * @param q The quadrature point index
   * @param X Array containing the coordinates of the support points.
   * @param func Callback function accepting three parameters: i, j and R_ij
   */
  template< typename FUNC >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void computeStiffnessTerm( localIndex const q,
                                    real64 const (&X)[8][3],
                                    FUNC && func );

  /**
   * @brief computes the matrix B in the case of quasi-stiffness (e.g. for pseudo-acoustic case), defined as J^{-T}A_xy J^{-1}/det(J), where
   * J is the Jacobian matrix, and A_xy is a zero matrix except on A_xy(1,1) = 1 and A_xy(2,2) = 1.
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param qc The 1d quadrature point index in xi2 direction (0,1)
   * @param X Array containing the coordinates of the support points.
   * @param J Array to store the Jacobian
   * @param B Array to store the matrix B, in Voigt notation
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void computeBxyMatrix( int const qa,
                                int const qb,
                                int const qc,
                                real64 const (&X)[8][3],
                                real64 ( &J )[3][3],
                                real64 ( &B )[6] );

  /**
   * @brief computes the non-zero contributions of the d.o.f. indexed by q to the
   *   partial-stiffness matrix R, i.e., the superposition matrix of first derivatives in x and y
   *   of the shape functions. Warning, the matrix B is obtained by computeBxyMatrix instead of usual one.
   * @param q The quadrature point index
   * @param X Array containing the coordinates of the support points.
   * @param func Callback function accepting three parameters: i, j and R_ij
   */
  template< typename FUNC >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void computeStiffnessxyTerm( localIndex const q,
                                      real64 const (&X)[8][3],
                                      FUNC && func );

  /**
   * @brief computes the matrix B in the case of quasi-stiffness (e.g. for pseudo-acoustic case), defined as J^{-T}A_z J^{-1}/det(J), where
   * J is the Jacobian matrix, and A_z is a zero matrix except on A_z(3,3) = 1.
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param qc The 1d quadrature point index in xi2 direction (0,1)
   * @param X Array containing the coordinates of the support points.
   * @param J Array to store the Jacobian
   * @param B Array to store the matrix B, in Voigt notation
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void computeBzMatrix( int const qa,
                               int const qb,
                               int const qc,
                               real64 const (&X)[8][3],
                               real64 ( &J )[3][3],
                               real64 ( &B )[6] );

  /**
   * @brief computes the non-zero contributions of the d.o.f. indexed by q to the
   *   partial-stiffness matrix R, i.e., the superposition matrix of first derivatives in z only
   *   of the shape functions. Warning, the matrix B is obtained by computeBzMatrix instead of usual one.
   * @param q The quadrature point index
   * @param X Array containing the coordinates of the support points.
   * @param func Callback function accepting three parameters: i, j and R_ij
   */
  template< typename FUNC >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void computeStiffnesszTerm( localIndex const q,
                                     real64 const (&X)[8][3],
                                     FUNC && func );

/**
 * @brief Computes the "Grad(Phi)*B*Grad(Phi)" coefficient of the stiffness term. The matrix B must be provided and Phi denotes a basis
 * function.
 * @param qa The 1d quadrature point index in xi0 direction (0,1)
 * @param qb The 1d quadrature point index in xi1 direction (0,1)
 * @param qc The 1d quadrature point index in xi2 direction (0,1)
 * @param B Array of the B matrix, in Voigt notation
 * @param func Callback function accepting three parameters: i, j and R_ij
 */
  template< typename FUNC >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  computeGradPhiBGradPhi( int const qa,
                          int const qb,
                          int const qc,
                          real64 const (&B)[6],
                          FUNC && func );

  /**
   * @brief computes the non-zero contributions of the d.o.f. indexd by q to the
   *   x-part of the first order stiffness matrix R, i.e., the matrix composed of the
   *   the product of first derivatives of one shape function i and the shape function j itself.
   * @param q The quadrature point index
   * @param X Array containing the coordinates of the support points.
   * @param func Callback function accepting three parameters: i, j and R_ij
   */
  template< typename FUNC >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void computeFirstOrderStiffnessTermX( localIndex const q,
                                               real64 const (&X)[8][3],
                                               FUNC && func );
  /**
   * @brief computes the non-zero contributions of the d.o.f. indexd by q to the
   *   y-part of the first order stiffness matrix R, i.e., the matrix composed of the
   *   the product of first derivatives of one shape function i and the shape function j itself.
   * @param q The quadrature point index
   * @param X Array containing the coordinates of the support points.
   * @param func Callback function accepting three parameters: i, j and R_ij
   */
  template< typename FUNC >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void computeFirstOrderStiffnessTermY( localIndex const q,
                                               real64 const (&X)[8][3],
                                               FUNC && func );
  /**
   * @brief computes the non-zero contributions of the d.o.f. indexd by q to the
   *   z-part of the first order stiffness matrix R, i.e., the matrix composed of the
   *   the product of first derivatives of one shape function i and the shape function j itself.
   * @param q The quadrature point index
   * @param X Array containing the coordinates of the support points.
   * @param func Callback function accepting three parameters: i, j and R_ij
   */
  template< typename FUNC >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void computeFirstOrderStiffnessTermZ( localIndex const q,
                                               real64 const (&X)[8][3],
                                               FUNC && func );
  /**
   * @brief computes the non-zero contributions of the d.o.f. indexd by q to the
   *   stiffness matrix R for the elastic case, i.e., the superposition matrix of first derivatives
   *   of the shape functions. This callback returns the two indices indices i and j of matrix R and the value
   *   R[i][j] associated to those two indices.
   * @param q The quadrature point index
   * @param X Array containing the coordinates of the support points.
   * @param stiffnessVal Callback function accepting three parameters: i, j and R_ij
   */
  template< typename FUNC >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void computeFirstOrderStiffnessTerm( localIndex const q,
                                              real64 const (&X)[8][3],
                                              FUNC && stiffnessVal );


  /**
   * @brief Apply a Jacobian transformation matrix from the parent space to the
   *   physical space on the parent shape function derivatives, producing the
   *   shape function derivatives in the physical space.
   * @param q The quadrature point index
   * @param invJ The Jacobian transformation from parent->physical space.
   * @param gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void applyTransformationToParentGradients( int const q,
                                                    real64 const ( &invJ )[3][3],
                                                    real64 ( &gradN )[numNodes][3] );

  /**
   * @brief Apply a Jacobian transformation matrix from the parent space to the
   *   physical space on the parent shape function derivatives, producing the
   *   shape function derivatives in the physical space at a single point.
   * @param coords The parent coordinates at which to apply the transformation
   * @param invJ The Jacobian transformation from parent->physical space.
   * @param gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void applyTransformationToParentGradients( real64 const (&coords)[3],
                                                    real64 const ( &invJ )[3][3],
                                                    real64 ( &gradN )[numNodes][3] );


private:
  /// The length of one dimension of the parent element.
  constexpr static real64 parentLength = GL_BASIS::parentSupportCoord( 1 ) - GL_BASIS::parentSupportCoord( 0 );

  /// The volume of the element in the parent configuration.
  constexpr static real64 parentVolume = parentLength*parentLength*parentLength;
  /**
   * @brief Applies a function inside a generic loop in over the tensor product
   *   indices.
   * @tparam FUNC The type of function to call within the support loop.
   * @tparam PARAMS The parameter pack types to pass through to @p FUNC.
   * @param coords The parent coordinates at which to evaluate the shape function value
   * @param func The function to call within the support loop.
   * @param params The parameters to pass to @p func.
   */
  template< typename FUNC, typename ... PARAMS >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void supportLoop( real64 const (&coords)[3],
                           FUNC && func,
                           PARAMS &&... params );
  /**
   * @brief Applies a function inside a generic loop in over the tensor product
   *   indices.
   * @tparam FUNC The type of function to call within the support loop.
   * @tparam PARAMS The parameter pack types to pass through to @p FUNC.
   * @param q The quadrature node at which to evaluate the shape function value
   * @param func The function to call within the support loop.
   * @param params The parameters to pass to @p func.
   */
  template< typename FUNC, typename ... PARAMS >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void supportLoop( localIndex const q,
                           FUNC && func,
                           PARAMS &&... params );

};

/// @cond Doxygen_Suppress


template< typename GL_BASIS >
template< typename FUNC, typename ... PARAMS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BB_Tetrahedron< GL_BASIS >::supportLoop( real64 const (&coords)[3],
                                                              FUNC && func,
                                                              PARAMS &&... params )
{
  for( int c=0; c<num1dNodes; ++c )
  {
    for( int b=0; b<num1dNodes; ++b )
    {
      for( int a=0; a<num1dNodes; ++a )
      {
        real64 const dNdXi[3] = { GL_BASIS::gradient( a, coords[0] )*
                                  GL_BASIS::value( b, coords[1] )*
                                  GL_BASIS::value( c, coords[2] ),
                                  GL_BASIS::value( a, coords[0] )*
                                  GL_BASIS::gradient( b, coords[1] )*
                                  GL_BASIS::value( c, coords[2] ),
                                  GL_BASIS::value( a, coords[0] )*
                                  GL_BASIS::value( b, coords[1] )*
                                  GL_BASIS::gradient( c, coords[2] )};

        localIndex const nodeIndex = GL_BASIS::TensorProduct3D::linearIndex( a, b, c );

        func( dNdXi, nodeIndex, std::forward< PARAMS >( params )... );
      }
    }
  }
}

template< typename GL_BASIS >
template< typename FUNC, typename ... PARAMS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BB_Tetrahedron< GL_BASIS >::supportLoop( localIndex const q,
                                                              FUNC && func,
                                                              PARAMS &&... params )
{
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );
  for( int c=0; c<num1dNodes; ++c )
  {
    for( int b=0; b<num1dNodes; ++b )
    {
      for( int a=0; a<num1dNodes; ++a )
      {
        real64 const dNdXi[3] = { (b == qb && c == qc ) ? basisGradientAt( a, qa ) : 0,
                                  (a == qa && c == qc ) ? basisGradientAt( b, qb ) : 0,
                                  (a == qa && b == qb ) ? basisGradientAt( c, qc ) : 0 };

        localIndex const nodeIndex = GL_BASIS::TensorProduct3D::linearIndex( a, b, c );

        func( dNdXi, nodeIndex, std::forward< PARAMS >( params )... );
      }
    }
  }
}

//*************************************************************************************************

template< typename GL_BASIS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64
BB_Tetrahedron< GL_BASIS >::calcGradN( localIndex const q,
                                                            real64 const (&X)[numNodes][3],
                                                            real64 (& gradN)[numNodes][3] )
{
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );
  real64 Xmesh[8][3] = {{0}};
  for( int k = 0; k < 8; k++ )
  {
    const localIndex nodeIndex = meshIndexToLinearIndex3D( k );
    for( int i = 0; i < 3; i++ )
    {
      Xmesh[ k ][ i ] = X[ nodeIndex ][ i ];
    }
  }
  real64 J[3][3] = {{0}};

  jacobianTransformation( qa, qb, qc, Xmesh, J );

  real64 const detJ = LvArray::tensorOps::invert< 3 >( J );

  applyTransformationToParentGradients( q, J, gradN );

  return detJ;
}
//*************************************************************************************************
template< typename GL_BASIS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64
BB_Tetrahedron< GL_BASIS >::calcGradN( real64 const (&coords)[3],
                                                            real64 const (&X)[numNodes][3],
                                                            real64 (& gradN)[numNodes][3] )
{
  real64 J[3][3] = {{0}};

  jacobianTransformation( coords, X, J );

  real64 const detJ = LvArray::tensorOps::invert< 3 >( J );

  applyTransformationToParentGradients( coords, J, gradN );

  return detJ;
}
template< typename GL_BASIS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 BB_Tetrahedron< GL_BASIS >::
calcGradN( localIndex const q,
           real64 const (&X)[numNodes][3],
           StackVariables const & GEOS_UNUSED_PARAM( stack ),
           real64 ( & gradN )[numNodes][3] )
{
  return calcGradN( q, X, gradN );
}

//*************************************************************************************************
#if __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#endif

template< typename GL_BASIS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BB_Tetrahedron< GL_BASIS >::
jacobianTransformation( int const qa,
                        int const qb,
                        int const qc,
                        real64 const (&X)[8][3],
                        real64 ( & J )[3][3] )
{
  for( int k = 0; k < 8; k++ )
  {
    const int ka = k % 2;
    const int kb = ( k % 4 ) / 2;
    const int kc = k / 4;
    for( int j = 0; j < 3; j++ )
    {
      real64 jacCoeff = jacobianCoefficient1D( qa, 0, ka, j ) *
                        jacobianCoefficient1D( qb, 1, kb, j ) *
                        jacobianCoefficient1D( qc, 2, kc, j );
      for( int i = 0; i < 3; i++ )
      {
        J[i][j] +=  jacCoeff * X[k][i];
      }
    }
  }
}

template< typename GL_BASIS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BB_Tetrahedron< GL_BASIS >::
jacobianTransformation( real64 const (&coords)[3],
                        real64 const (&X)[numNodes][3],
                        real64 ( & J )[3][3] )
{
  supportLoop( coords, [] GEOS_HOST_DEVICE ( real64 const (&dNdXi)[3],
                                             int const nodeIndex,
                                             real64 const (&X)[numNodes][3],
                                             real64 (& J)[3][3] )
  {
    real64 const * const GEOS_RESTRICT Xnode = X[nodeIndex];
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        J[i][j] = J[i][j] + dNdXi[ j ] * Xnode[i];
      }
    }
  }, X, J );
}

template< typename GL_BASIS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BB_Tetrahedron< GL_BASIS >::
trilinearInterp( real64 const alpha,
                 real64 const beta,
                 real64 const gamma,
                 real64 const (&X)[8][3],
                 real64 (& coords)[3] )
{
  for( int i=0; i<3; i++ )
  {
    coords[i] = X[0][i]*( 1.0-alpha )*( 1.0-beta )*( 1.0-gamma )+
                X[1][i]*    alpha    *( 1.0-beta )*( 1.0-gamma )+
                X[2][i]*( 1.0-alpha )*    beta    *( 1.0-gamma )+
                X[3][i]*    alpha    *    beta    *( 1.0-gamma )+
                X[4][i]*( 1.0-alpha )*( 1.0-beta )*  gamma+
                X[5][i]*    alpha    *( 1.0-beta )*  gamma+
                X[6][i]*( 1.0-alpha )*    beta    *  gamma+
                X[7][i]*    alpha    *    beta    *  gamma;
  }
}


template< typename GL_BASIS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BB_Tetrahedron< GL_BASIS >::
computeLocalCoords( real64 const (&Xmesh)[8][3],
                    real64 const (&X)[numNodes][3] )
{
  int qa, qb, qc;
  for( int q=0; q<numNodes; q++ )
  {
    GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );
    real64 alpha = ( GL_BASIS::parentSupportCoord( qa ) + 1.0 ) / 2.0;
    real64 beta = ( GL_BASIS::parentSupportCoord( qb ) + 1.0 ) / 2.0;
    real64 gamma = ( GL_BASIS::parentSupportCoord( qc ) + 1.0 ) / 2.0;
    trilinearInterp( alpha, beta, gamma, Xmesh, X[q] );
  }
}

template< typename GL_BASIS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BB_Tetrahedron< GL_BASIS >::
jacobianTransformation2d( int const qa,
                          int const qb,
                          real64 const (&X)[4][3],
                          real64 ( & J )[3][2] )
{
  for( int k = 0; k < 4; k++ )
  {
    int ka = k % 2;
    int kb = k / 2;
    for( int j = 0; j < 2; j++ )
    {
      real64 jacCoeff = jacobianCoefficient1D( qa, 0, ka, j ) *
                        jacobianCoefficient1D( qb, 1, kb, j );
      for( int i = 0; i < 3; i++ )
      {
        J[i][j] +=  jacCoeff * X[k][i];
      }
    }
  }
}

template< typename GL_BASIS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64
BB_Tetrahedron< GL_BASIS >::
computeMassTerm( localIndex const q,
                 real64 const (&X)[8][3] )

template< typename GL_BASIS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64
BB_Tetrahedron< GL_BASIS >::
computeDampingTerm( localIndex const q,
                    real64 const (&X)[4][3] )
{
  int qa, qb;
  GL_BASIS::TensorProduct2D::multiIndex( q, qa, qb );
  const real64 w2D = GL_BASIS::weight( qa )*GL_BASIS::weight( qb );
  real64 B[3];
  real64 J[3][2] = {{0}};
  jacobianTransformation2d( qa, qb, X, J );
  // compute J^T.J, using Voigt notation for B
  B[0] = J[0][0]*J[0][0]+J[1][0]*J[1][0]+J[2][0]*J[2][0];
  B[1] = J[0][1]*J[0][1]+J[1][1]*J[1][1]+J[2][1]*J[2][1];
  B[2] = J[0][0]*J[0][1]+J[1][0]*J[1][1]+J[2][0]*J[2][1];
  return sqrt( LvArray::math::abs( LvArray::tensorOps::symDeterminant< 2 >( B ) ) )*w2D;
}

template< typename GL_BASIS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BB_Tetrahedron< GL_BASIS >::
computeBMatrix( int const qa,
                int const qb,
                int const qc,
                real64 const (&X)[8][3],
                real64 (& J)[3][3],
                real64 (& B)[6] )
{
  jacobianTransformation( qa, qb, qc, X, J );
  real64 const detJ = LvArray::tensorOps::determinant< 3 >( J );

  // compute J^T.J/det(J), using Voigt notation for B
  B[0] = (J[0][0]*J[0][0]+J[1][0]*J[1][0]+J[2][0]*J[2][0])/detJ;
  B[1] = (J[0][1]*J[0][1]+J[1][1]*J[1][1]+J[2][1]*J[2][1])/detJ;
  B[2] = (J[0][2]*J[0][2]+J[1][2]*J[1][2]+J[2][2]*J[2][2])/detJ;
  B[3] = (J[0][1]*J[0][2]+J[1][1]*J[1][2]+J[2][1]*J[2][2])/detJ;
  B[4] = (J[0][0]*J[0][2]+J[1][0]*J[1][2]+J[2][0]*J[2][2])/detJ;
  B[5] = (J[0][0]*J[0][1]+J[1][0]*J[1][1]+J[2][0]*J[2][1])/detJ;

  // compute detJ*J^{-1}J^{-T}
  LvArray::tensorOps::symInvert< 3 >( B );
}

template< typename GL_BASIS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BB_Tetrahedron< GL_BASIS >::
computeBzMatrix( int const qa,
                 int const qb,
                 int const qc,
                 real64 const (&X)[8][3],
                 real64 (& J)[3][3],
                 real64 (& B)[6] )
{
  jacobianTransformation( qa, qb, qc, X, J );
  real64 const detJ = LvArray::tensorOps::determinant< 3 >( J );

  real64 Jinv[3][3] = {{0}};
  LvArray::tensorOps::invert< 3 >( Jinv, J );

  // compute det(J)*J^{-1}Az*J^{-T}, using Voigt notation for B
  B[0] = detJ*(Jinv[0][2]*Jinv[0][2]);
  B[1] = detJ*(Jinv[1][2]*Jinv[1][2]);
  B[2] = detJ*(Jinv[2][2]*Jinv[2][2]);
  B[3] = detJ*(Jinv[1][2]*Jinv[2][2]);
  B[4] = detJ*(Jinv[0][2]*Jinv[2][2]);
  B[5] = detJ*(Jinv[0][2]*Jinv[1][2]);
}

template< typename GL_BASIS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BB_Tetrahedron< GL_BASIS >::
computeBxyMatrix( int const qa,
                  int const qb,
                  int const qc,
                  real64 const (&X)[8][3],
                  real64 (& J)[3][3],
                  real64 (& B)[6] )
{
  jacobianTransformation( qa, qb, qc, X, J );
  real64 const detJ = LvArray::tensorOps::determinant< 3 >( J );

  real64 Jinv[3][3] = {{0}};
  LvArray::tensorOps::invert< 3 >( Jinv, J );

  // compute det(J)*J^{-1}Axy*J^{-T}, using Voigt notation for B
  B[0] = detJ*(Jinv[0][0]*Jinv[0][0] + Jinv[0][1]*Jinv[0][1]);
  B[1] = detJ*(Jinv[1][1]*Jinv[1][1] + Jinv[1][0]*Jinv[1][0]);
  B[2] = detJ*(Jinv[2][0]*Jinv[2][0] + Jinv[2][1]*Jinv[2][1]);
  B[3] = detJ*(Jinv[1][0]*Jinv[2][0] + Jinv[1][1]*Jinv[2][1]);
  B[4] = detJ*(Jinv[0][0]*Jinv[2][0] + Jinv[0][1]*Jinv[2][1]);
  B[5] = detJ*(Jinv[0][0]*Jinv[1][0] + Jinv[0][1]*Jinv[1][1]);
}

template< typename GL_BASIS >
template< typename FUNC >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BB_Tetrahedron< GL_BASIS >::
computeGradPhiBGradPhi( int const qa,
                        int const qb,
                        int const qc,
                        real64 const (&B)[6],
                        FUNC && func )
{
  const real64 w = GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc );
  for( int i=0; i<num1dNodes; i++ )
  {
    const int ibc = GL_BASIS::TensorProduct3D::linearIndex( i, qb, qc );
    const int aic = GL_BASIS::TensorProduct3D::linearIndex( qa, i, qc );
    const int abi = GL_BASIS::TensorProduct3D::linearIndex( qa, qb, i );
    const real64 gia = basisGradientAt( i, qa );
    const real64 gib = basisGradientAt( i, qb );
    const real64 gic = basisGradientAt( i, qc );
    for( int j=0; j<num1dNodes; j++ )
    {
      const int jbc = GL_BASIS::TensorProduct3D::linearIndex( j, qb, qc );
      const int ajc = GL_BASIS::TensorProduct3D::linearIndex( qa, j, qc );
      const int abj = GL_BASIS::TensorProduct3D::linearIndex( qa, qb, j );
      const real64 gja = basisGradientAt( j, qa );
      const real64 gjb = basisGradientAt( j, qb );
      const real64 gjc = basisGradientAt( j, qc );
      // diagonal terms
      const real64 w0 = w * gia * gja;
      func( ibc, jbc, w0 * B[0] );
      const real64 w1 = w * gib * gjb;
      func( aic, ajc, w1 * B[1] );
      const real64 w2 = w * gic * gjc;
      func( abi, abj, w2 * B[2] );
      // off-diagonal terms
      const real64 w3 = w * gib * gjc;
      func( aic, abj, w3 * B[3] );
      func( abj, aic, w3 * B[3] );
      const real64 w4 = w * gia * gjc;
      func( ibc, abj, w4 * B[4] );
      func( abj, ibc, w4 * B[4] );
      const real64 w5 = w * gia * gjb;
      func( ibc, ajc, w5 * B[5] );
      func( ajc, ibc, w5 * B[5] );
    }
  }
}

template< typename GL_BASIS >
template< typename FUNC >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BB_Tetrahedron< GL_BASIS >::
computeStiffnessxyTerm( localIndex const q,
                        real64 const (&X)[8][3],
                        FUNC && func )
{
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );
  real64 B[6] = {0};
  real64 J[3][3] = {{0}};
  computeBxyMatrix( qa, qb, qc, X, J, B ); // The only change!
  computeGradPhiBGradPhi( qa, qb, qc, B, func );
}

template< typename GL_BASIS >
template< typename FUNC >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BB_Tetrahedron< GL_BASIS >::
computeStiffnesszTerm( localIndex const q,
                       real64 const (&X)[8][3],
                       FUNC && func )
{
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );
  real64 B[6] = {0};
  real64 J[3][3] = {{0}};
  computeBzMatrix( qa, qb, qc, X, J, B ); // The only change!
  computeGradPhiBGradPhi( qa, qb, qc, B, func );
}

template< typename GL_BASIS >
template< typename FUNC >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BB_Tetrahedron< GL_BASIS >::
computeStiffnessTerm( localIndex const q,
                      real64 const (&X)[8][3],
                      FUNC && func )
{
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );
  real64 B[6] = {0};
  real64 J[3][3] = {{0}};
  computeBMatrix( qa, qb, qc, X, J, B );
  computeGradPhiBGradPhi( qa, qb, qc, B, func );
}

template< typename GL_BASIS >
template< typename FUNC >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BB_Tetrahedron< GL_BASIS >::
computeFirstOrderStiffnessTerm( localIndex const q,
                                real64 const (&X)[8][3],
                                FUNC && func )
{
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );
  real64 J[3][3] = {{0}};
  jacobianTransformation( qa, qb, qc, X, J );
  real64 const detJ = LvArray::tensorOps::invert< 3 >( J );
  const real64 w = GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc );
  for( int i=0; i<num1dNodes; i++ )
  {
    const int ibc = GL_BASIS::TensorProduct3D::linearIndex( i, qb, qc );
    const int aic = GL_BASIS::TensorProduct3D::linearIndex( qa, i, qc );
    const int abi = GL_BASIS::TensorProduct3D::linearIndex( qa, qb, i );
    const real64 gia = basisGradientAt( i, qa );
    const real64 gib = basisGradientAt( i, qb );
    const real64 gic = basisGradientAt( i, qc );
    for( int j=0; j<num1dNodes; j++ )
    {
      const int jbc = GL_BASIS::TensorProduct3D::linearIndex( j, qb, qc );
      const int ajc = GL_BASIS::TensorProduct3D::linearIndex( qa, j, qc );
      const int abj = GL_BASIS::TensorProduct3D::linearIndex( qa, qb, j );
      const real64 gja = basisGradientAt( j, qa );
      const real64 gjb = basisGradientAt( j, qb );
      const real64 gjc = basisGradientAt( j, qc );
      // diagonal terms
      const real64 w00 = w * gia * gja;
      func( ibc, jbc, w00 * detJ, J, 0, 0 );
      const real64 w11 = w * gib * gjb;
      func( aic, ajc, w11 * detJ, J, 1, 1 );
      const real64 w22 = w * gic * gjc;
      func( abi, abj, w22 * detJ, J, 2, 2 );
      // off-diagonal terms
      const real64 w12 = w * gib * gjc;
      func( aic, abj, w12 * detJ, J, 1, 2 );
      func( abj, aic, w12 * detJ, J, 2, 1 );
      const real64 w02 = w * gia * gjc;
      func( ibc, abj, w02 * detJ, J, 0, 2 );
      func( abj, ibc, w02 * detJ, J, 2, 0 );
      const real64 w01 = w * gia * gjb;
      func( ibc, ajc, w01 * detJ, J, 0, 1 );
      func( ajc, ibc, w01 * detJ, J, 1, 0 );
    }
  }
}

template< typename GL_BASIS >
template< typename FUNC >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BB_Tetrahedron< GL_BASIS >::
computeFirstOrderStiffnessTermX( localIndex const q,
                                 real64 const (&X)[8][3],
                                 FUNC && func )
{
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );
  real64 J[3][3] = {{0}};
  jacobianTransformation( qa, qb, qc, X, J );
  const real64 detJ = LvArray::tensorOps::invert< 3 >( J );
  const real64 w = GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc );

  for( int i1 = 0; i1 < num1dNodes; ++i1 )
  {
    auto val = w * basisGradientAt( i1, qa );
    func( GL_BASIS::TensorProduct3D::linearIndex( i1, qb, qc ), q, detJ*J[0][0]*val, detJ*J[0][1]*val, detJ*J[0][2]*val );
  }

}

template< typename GL_BASIS >
template< typename FUNC >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BB_Tetrahedron< GL_BASIS >::
computeFirstOrderStiffnessTermY( localIndex const q,
                                 real64 const (&X)[8][3],
                                 FUNC && func )
{
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );
  real64 J[3][3] = {{0}};
  jacobianTransformation( qa, qb, qc, X, J );
  const real64 detJ = LvArray::tensorOps::invert< 3 >( J );
  const real64 w = GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc );

  for( int i2 = 0; i2 < num1dNodes; ++i2 )
  {
    auto val = w * basisGradientAt( i2, qb );
    func( GL_BASIS::TensorProduct3D::linearIndex( qa, i2, qc ), q, detJ*J[1][0]*val, detJ*J[1][1]*val, detJ*J[1][2]*val );
  }
}

template< typename GL_BASIS >
template< typename FUNC >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BB_Tetrahedron< GL_BASIS >::
computeFirstOrderStiffnessTermZ( localIndex const q,
                                 real64 const (&X)[8][3],
                                 FUNC && func )
{
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );
  real64 J[3][3] = {{0}};
  jacobianTransformation( qa, qb, qc, X, J );
  const real64 detJ = LvArray::tensorOps::invert< 3 >( J );
  const real64 w = GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc );

  for( int i3 = 0; i3 < num1dNodes; ++i3 )
  {
    auto val = w * basisGradientAt( i3, qc );
    func( GL_BASIS::TensorProduct3D::linearIndex( qa, qb, i3 ), q, detJ*J[2][0]*val, detJ*J[2][1]*val, detJ*J[2][2]*val );
  }
}

//*************************************************************************************************
template< typename GL_BASIS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BB_Tetrahedron< GL_BASIS >::
applyTransformationToParentGradients( int const q,
                                      real64 const ( &invJ )[3][3],
                                      real64 (& gradN)[numNodes][3] )
{
  supportLoop( q, [] GEOS_HOST_DEVICE ( real64 const (&dNdXi)[3],
                                        int const nodeIndex,
                                        real64 const (&invJ)[3][3],
                                        real64 (& gradN)[numNodes][3] )
  {
    // smaller register footprint by manually unrolling the for loops.
    gradN[nodeIndex][0] = dNdXi[0] * invJ[0][0] + dNdXi[1] * invJ[1][0] + dNdXi[2] * invJ[2][0];
    gradN[nodeIndex][1] = dNdXi[0] * invJ[0][1] + dNdXi[1] * invJ[1][1] + dNdXi[2] * invJ[2][1];
    gradN[nodeIndex][2] = dNdXi[0] * invJ[0][2] + dNdXi[1] * invJ[1][2] + dNdXi[2] * invJ[2][2];


  }, invJ, gradN );
}

//*************************************************************************************************
template< typename GL_BASIS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
BB_Tetrahedron< GL_BASIS >::
applyTransformationToParentGradients( real64 const (&coords)[3],
                                      real64 const ( &invJ )[3][3],
                                      real64 (& gradN)[numNodes][3] )
{
  supportLoop( coords, [] GEOS_HOST_DEVICE ( real64 const (&dNdXi)[3],
                                             int const nodeIndex,
                                             real64 const (&invJ)[3][3],
                                             real64 (& gradN)[numNodes][3] )
  {
    gradN[nodeIndex][0] = dNdXi[0] * invJ[0][0] + dNdXi[1] * invJ[1][0] + dNdXi[2] * invJ[2][0];
    gradN[nodeIndex][1] = dNdXi[0] * invJ[0][1] + dNdXi[1] * invJ[1][1] + dNdXi[2] * invJ[2][1];
    gradN[nodeIndex][2] = dNdXi[0] * invJ[0][2] + dNdXi[1] * invJ[1][2] + dNdXi[2] * invJ[2][2];
  }, invJ, gradN );
}

template< typename GL_BASIS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64
BB_Tetrahedron< GL_BASIS >::
transformedQuadratureWeight( localIndex const q,
                             real64 const (&X)[numNodes][3] )
{
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );
  real64 J[3][3] = {{0}};

  jacobianTransformation( qa, qb, qc, X, J );

  return LvArray::tensorOps::determinant< 3 >( J );
}



template< typename GL_BASIS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void BB_Tetrahedron< GL_BASIS >::
symmetricGradient( int const q,
                   real64 const (&invJ)[3][3],
                   real64 const (&var)[numNodes][3],
                   real64 (& grad)[6] )
{
  supportLoop( q, [] GEOS_HOST_DEVICE ( real64 const (&dNdXi)[3],
                                        int const nodeIndex,
                                        real64 const (&invJ)[3][3],
                                        real64 const (&var)[numNodes][3],
                                        real64 (& grad)[6] )
  {

    real64 gradN[3] = {0, 0, 0};
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        gradN[i] = gradN[i] + dNdXi[ j ] * invJ[j][i];
      }
    }

    grad[0] = grad[0] + gradN[0] * var[ nodeIndex ][0];
    grad[1] = grad[1] + gradN[1] * var[ nodeIndex ][1];
    grad[2] = grad[2] + gradN[2] * var[ nodeIndex ][2];
    grad[3] = grad[3] + gradN[2] * var[ nodeIndex ][1] + gradN[1] * var[ nodeIndex ][2];
    grad[4] = grad[4] + gradN[2] * var[ nodeIndex ][0] + gradN[0] * var[ nodeIndex ][2];
    grad[5] = grad[5] + gradN[1] * var[ nodeIndex ][0] + gradN[0] * var[ nodeIndex ][1];
  }, invJ, var, grad );
}

template< typename GL_BASIS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void BB_Tetrahedron< GL_BASIS >::
plusGradNajAij( int const q,
                real64 const (&invJ)[3][3],
                real64 const (&var)[6],
                real64 (& R)[numNodes][3] )
{
  supportLoop( q,
               [] GEOS_HOST_DEVICE
                 ( real64 const (&dNdXi)[3],
                 int const nodeIndex,
                 real64 const (&invJ)[3][3],
                 real64 const (&var)[6],
                 real64 (& R)[numNodes][3] )
  {

    real64 gradN[3] = {0, 0, 0};
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        gradN[i] = gradN[i] + dNdXi[ j ] * invJ[j][i];
      }
    }
    R[ nodeIndex ][ 0 ] = R[ nodeIndex ][ 0 ] - var[ 0 ] * gradN[ 0 ] - var[ 5 ] * gradN[ 1 ] - var[ 4 ] * gradN[ 2 ];
    R[ nodeIndex ][ 1 ] = R[ nodeIndex ][ 1 ] - var[ 5 ] * gradN[ 0 ] - var[ 1 ] * gradN[ 1 ] - var[ 3 ] * gradN[ 2 ];
    R[ nodeIndex ][ 2 ] = R[ nodeIndex ][ 2 ] - var[ 4 ] * gradN[ 0 ] - var[ 3 ] * gradN[ 1 ] - var[ 2 ] * gradN[ 2 ];
  }, invJ, var, R );
}



template< typename GL_BASIS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void BB_Tetrahedron< GL_BASIS >::
gradient( int const q,
          real64 const (&invJ)[3][3],
          real64 const (&var)[numNodes][3],
          real64 (& grad)[3][3] )
{
  supportLoop( q, [] GEOS_HOST_DEVICE ( real64 const (&dNdXi)[3],
                                        int const nodeIndex,
                                        real64 const (&invJ)[3][3],
                                        real64 const (&var)[numNodes][3],
                                        real64 (& grad)[3][3] )
  {
    for( int i = 0; i < 3; ++i )
    {
      real64 gradN=0.0;;
      for( int j = 0; j < 3; ++j )
      {
        gradN = gradN + dNdXi[ j ] * invJ[j][i];
      }
      for( int k = 0; k < 3; ++k )
      {
        grad[k][i] = grad[k][i] + gradN * var[ nodeIndex ][k];
      }
    }
  }, invJ, var, grad );
}
/**
 * Fixed-degree classes
 */
using BB1_Tetrahedron = BB_Tetrahedron< 1 >;
using BB2_Tetrahedron = BB_Tetrahedron< 2 >;
using BB3_Tetrahedron = BB_Tetrahedron< 3 >;
using BB4_Tetrahedron = BB_Tetrahedron< 4 >;
using BB5_Tetrahedron = BB_Tetrahedron< 5 >;
using BB6_Tetrahedron = BB_Tetrahedron< 6 >;
using BB7_Tetrahedron = BB_Tetrahedron< 7 >;
using BB8_Tetrahedron = BB_Tetrahedron< 8 >;
using BB9_Tetrahedron = BB_Tetrahedron< 9 >;
using BB10_Tetrahedron = BB_Tetrahedron< 10 >;

/// @endcond

#if __GNUC__
#pragma GCC diagnostic pop
#endif
#undef PARENT_GRADIENT_METHOD
}
}

#endif // GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_BBTETRAHEDRON_HPP_
