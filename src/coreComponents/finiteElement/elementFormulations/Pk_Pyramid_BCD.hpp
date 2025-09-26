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
 * @file Pk_Pyramid_BCD.hpp
 */

#ifndef GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_PKPYRAMIDBCD_HPP_
#define GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_PKPYRAMIDBCD_HPP_

#include "FiniteElementBase.hpp"
#include "denseLinearAlgebra/interfaces/blaslapack/BlasLapackLA.hpp"
#include <utility>



namespace geos
{
namespace finiteElement
{

namespace
{
template< int ORDER >
constexpr int Pk_Pyramid_BCD_NumNodes = ( ORDER + 1 ) * ( ORDER + 2 ) * ( 2 * ORDER + 3 ) / 6;
}


/**
 * This class is the basis class for the pyramid finite element cells with
 * shape functions defined using the Bergot-Cohen-Durufle method.
 * All the degree-specific versions (P1, P2, P3, ...) are defined at the end of this file.
 * For now only P1 is implemented, P2 and higher will be implemented later.
 */
template< int ORDER >
class Pk_Pyramid_BCD_impl : public FiniteElementBase_impl< Pk_Pyramid_BCD_NumNodes< ORDER >,
                                                           5,
                                                           Pk_Pyramid_BCD_NumNodes< ORDER > >
{
public:
  /// Convenience alias for the base class.
  using Base = FiniteElementBase_impl< Pk_Pyramid_BCD_NumNodes< ORDER >,
                                       5,
                                       Pk_Pyramid_BCD_NumNodes< ORDER > >;

  /// Stack variables for the element.
  using StackVariables = typename Base::StackVariables;

  /// Mesh data structure for the element.
  template< typename SubregionType >
  using MeshData = typename Base::template MeshData< SubregionType >;

  /// The number of shape functions per element.
  using Base::numNodes;

  /// The number of quadrature points per element.
  using Base::numQuadraturePoints;

  /// The maximum number of support points per element.
  using Base::maxSupportPoints;

  /// The order of the finite element.
  static constexpr int order = ORDER;

  /// The number of modal points per element.
  constexpr static localIndex numModes = numNodes;

  /// @copydoc FiniteElementBase_impl::getNumQuadraturePoints()
  GEOS_HOST_DEVICE
  static localIndex getNumQuadraturePoints()
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

  /// @copydoc FiniteElementBase_impl::getNumSupportPoints()
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static localIndex getNumSupportPoints()
  {
    return numNodes;
  }

  /// @copydoc FiniteElementBase_impl::getMaxSupportPoints()
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static localIndex getMaxSupportPoints()
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
   * @brief Get the number of modal points.
   * @param stack Object that holds stack variables.
   * @return The number of modal points.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static localIndex getNumModes( StackVariables const & stack )
  {
    GEOS_UNUSED_VAR( stack );
    return numNodes;
  }


  /**
   * @brief The linear index associated to the given one-dimensional indices in the three directions
   * @param qa The index in the first direction
   * @param qb The index in the second direction
   * @param qc The index in the third direction
   * @return The linear index in 3D
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static localIndex linearIndex3DVal( const localIndex qa, localIndex const qb, localIndex const qc )
  {
    localIndex index = 0;

    for( int l = 0; l < qc; ++l )
    {
      int n = order + 1 - l;
      index += n * n;
    }
    int n_k = order + 1 - qc;
    index += qa + qb * n_k;

    return index;
  }


  /**
   * @brief Generate the coordinates of the support points
   * @param coords The array to fill with the coordinates of the support points
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static constexpr void generatePointsCoordinates( real64 (& coords)[numNodes][3] )
  {
    // Generate the coordinates of the support points based on the order
    if constexpr (ORDER == 1)
    {
      coords[0][0] = -1.0; coords[0][1] = 1.0; coords[0][2] = 0.0;
      coords[1][0] =  -1.0; coords[1][1] = -1.0; coords[1][2] = 0.0;
      coords[2][0] =  1.0; coords[2][1] =  -1.0; coords[2][2] = 0.0;
      coords[3][0] = 1.0; coords[3][1] =  1.0; coords[3][2] = 0.0;
      coords[4][0] =  0.0; coords[4][1] =  0.0; coords[4][2] = 1.0;
    }
    else if constexpr (ORDER == 2)
    {
      //TODO
    }

  }


  /**
   * @brief Generate the indexes for the modal shape functions
   * @param func The function to call with the generated indexes
   */
  template< typename FUNC >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static constexpr void generateIndexes( FUNC && func )
  {

    for( localIndex i = 0; i <= ORDER; ++i )
    {
      for( localIndex j = 0; j <= ORDER; ++j )
      {
        localIndex maxIj = LvArray::math::max( i, j );
        for( localIndex k = 0; k <= ORDER - maxIj; ++k )
        {
          func( i, j, k );
        }
      }
    }
  }

  /**
   * @brief Generate the indexes for the modal shape functions
   * @param func The function to call with the generated indexes
   */
  template< typename FUNC >
  GEOS_FORCE_INLINE
  static constexpr void generateIndexesHost( FUNC && func )
  {

    for( localIndex i = 0; i <= ORDER; ++i )
    {
      for( localIndex j = 0; j <= ORDER; ++j )
      {
        localIndex maxIj = LvArray::math::max( i, j );
        for( localIndex k = 0; k <= ORDER - maxIj; ++k )
        {
          func( i, j, k );
        }
      }
    }
  }

  /**
   * @brief Returns the determinant of the Jacobian of the element
   * @param[in] X The coordinates of the tetrahedron
   * @return the (absolute value of the) determinant of the Jacobian on the element
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 jacobianDeterminant( real64 const (&X)[5][3] )

  {
    real64 m[3][3] = {};
    for( int i = 0; i < 4; i++ )
    {
      for( int j = 0; j < 3; j++ )
      {
        m[ i ][ j ] = X[ i + 1 ][ j ] - X[ 0 ][ j ];
      }
    }
    return LvArray::math::abs( LvArray::tensorOps::determinant< 3 >( m ) );
  }


  /**
   * @brief Calculate the determinant of the jacobian on the face opposite to the given vertex
   * @param[in] face The index of the vertex opposite to the desired face
   * @param[in] X The coordinates of the pyramid
   * @return the (absolute value of the) determinant of the Jacobian on the face
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 faceJacobianDeterminant( localIndex face,
                                         real64 const (&X)[5][3] )
  {
    int i1 = ( face + 1 ) % 4;
    int i2 = ( face + 2 ) % 4;
    int i3 = ( face + 3 ) % 4;

    real64 ab[3] = { X[ i2 ][ 0 ] - X[ i1 ][ 0 ],
                     X[ i2 ][ 1 ] - X[ i1 ][ 1 ],
                     X[ i2 ][ 2 ] - X[ i1 ][ 2 ] };
    real64 ac[3] = { X[ i3 ][ 0 ] - X[ i1 ][ 0 ],
                     X[ i3 ][ 1 ] - X[ i1 ][ 1 ],
                     X[ i3 ][ 2 ] - X[ i1 ][ 2 ] };
    real64 term1 = ab[1] * ac[2] - ab[2] * ac[1];
    real64 term2 = ab[2] * ac[0] - ab[0] * ac[2];
    real64 term3 = ab[0] * ac[1] - ab[1] * ac[0];
    return LvArray::math::sqrt( ( term1 * term1 + term2 * term2 + term3 * term3 ) );

  }

  /**
   * @brief Calculate shape functions values for each support point at a
   *   quadrature/nodal point.
   * @param q Index of the quadrature point.
   * @param stack Variables allocated on the stack as filled by @ref setupStack.
   * @param N An array to pass back the shape function values for each support
   *   point.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void calcN( localIndex const q,
                     StackVariables const & stack,
                     real64 ( & N )[numNodes] )
  {
    GEOS_UNUSED_VAR( stack );
    return calcN( q, N );
  }


  /**
   * @brief Helper to compute Jacobi Polynomials.
   * @param n The degree of the polynomial.
   * @param alpha The alpha parameter of the Jacobi polynomial.
   * @param beta The beta parameter of the Jacobi polynomial.
   * @param x The point at which to evaluate the polynomial.
   * @return The value of the Jacobi polynomial at point x.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static constexpr real64 EvaluateJacobiPolynomial( localIndex const n, real64 const alpha, real64 const beta, real64 const x )
  {

    if( n == 0 )
    {
      real64 val = 1.0;
      return val;
    }

    if( n == 1 )
    {
      real64 val = 0.5 * (alpha - beta + (alpha + beta + 2) * x);
      return val;
    }

    // Initial conditions for the recurrence relation
    real64 P_prev2 = 1.0;  // P_0
    real64 P_prev1 = 0.5 * (alpha - beta + (alpha + beta + 2) * x);  // P_1
    real64 P_current = 0.0;

    // Recurrence relation to compute the Jacobi polynomial
    for( int k = 2; k < n+1; ++k )
    {
      real64 a_k = 2 * k * (k + alpha + beta) * (2 * k + alpha + beta - 2);
      real64 b_k = (2 * k + alpha + beta - 1) * (alpha * alpha - beta * beta);
      real64 c_k = (2 * k + alpha + beta - 1) * (2 * k + alpha + beta) * (2 * k + alpha + beta - 2);
      real64 d_k = 2 * (k + alpha - 1) * (k + beta - 1) * (2 * k + alpha + beta);

      P_current = ((b_k + c_k * x) * P_prev1 - d_k * P_prev2) / a_k;

      // Mise à jour pour l'itération suivante
      P_prev2 = P_prev1;
      P_prev1 = P_current;
    }

    return P_current;
  }

/**
 * @brief Evaluate the derivative of the Jacobi polynomial at a point x.
 * @param n The degree of the polynomial.
 * @param alpha The alpha parameter of the Jacobi polynomial.
 * @param beta The beta parameter of the Jacobi polynomial.
 * @param x The point at which to evaluate the derivative.
 * @return The value of the derivative of the Jacobi polynomial at point x.
 */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static constexpr real64 EvaluateJacobiPolynomialDerivative( localIndex const n, real64 const alpha, real64 const beta, real64 const x )
  {
    // Particular case for n = 0
    if( n == 0 )
    {
      return 0.0;
    }

    // dP_n^{(alpha,beta)}(x)/dx = (n+alpha+beta+1)/2 * P_{n-1}^{(alpha+1,beta+1)}(x)
    real64 coeff = 0.5 * (n + alpha + beta + 1.0);
    real64 jacobi_nm1 = EvaluateJacobiPolynomial( n-1, alpha+1.0, beta+1.0, x );

    return coeff * jacobi_nm1;
  }

  /**
   * @brief Calculate modal base functions values for a modal point (i,j,k) at a
   *   point X.
   * @param[in] i,j,k Indexes of a modal point.
   * @param[in] X Coordinates in reference pyramid (array of size 3)
   * @return the modal function value
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static constexpr real64 calcModal( localIndex const i, localIndex const j, localIndex const k, real64 const (&X)[3] )
  {

    real64 const x = X[0];
    real64 const y = X[1];
    real64 const z = X[2];

    real64 const epsilon = 1e-10; // Small value to avoid compairing floating point numbers directly
    if( LvArray::math::abs( z-1.0 ) < epsilon )
    {
      if( i == 0 && j == 0 )
      {
        return (k+2)*(k+1)/2.0;
      }
      else
      {
        return 0.0;
      }
    }

    real64 xi = x / (1.0 - z);
    real64 eta = y / (1.0 - z);
    real64 chi = 2.0 * z - 1.0;
    localIndex max_ij = LvArray::math::max( i, j );
    real64 P_i = EvaluateJacobiPolynomial( i, 0.0, 0.0, xi );
    real64 P_j = EvaluateJacobiPolynomial( j, 0.0, 0.0, eta );
    real64 P_k = EvaluateJacobiPolynomial( k, 2.0 * max_ij+2.0, 0.0, chi );

    return P_i * P_j * std::pow( 1.0 - z, max_ij ) * P_k;

  }

  /**
   * @brief Calculate modal base functions derivative values for a modal point (i,j,k) at a
   *  point X.
   * @param[in] i,j,k Indexes of a modal point.
   * @param[in] X Coordinates in reference pyramid (array of size 3)
   * @param[out] gradPsiX A real to pass back the modal function derivative values
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static constexpr void calcGradModal( int i, int j, int k, real64 const (&X)[3],
                                       real64 (& gradPsiX)[3] )
  {
    real64 x = X[0];
    real64 y = X[1];
    real64 z = X[2];

    real64 const epsilon = 1e-10; // Small value to avoid compairing floating point numbers directly
    if( LvArray::math::abs( z-1.0 ) < epsilon )
    {
      gradPsiX[0] = 0.0;
      gradPsiX[1] = 0.0;
      gradPsiX[2] = 0.0;
      return;
    }


    real64 xi = x / (1.0 - z);
    real64 eta = y / (1.0 - z);
    real64 chi = 2.0 * z - 1.0;
    localIndex m = LvArray::math::max( i, j );
    real64 P_i   = EvaluateJacobiPolynomial( i, 0.0, 0.0, xi );
    real64 P_j   = EvaluateJacobiPolynomial( j, 0.0, 0.0, eta );
    real64 P_k   = EvaluateJacobiPolynomial( k, 2.0 * m + 2.0, 0.0, chi );
    real64 dP_i  = EvaluateJacobiPolynomialDerivative( i, 0.0, 0.0, xi );
    real64 dP_j  = EvaluateJacobiPolynomialDerivative( j, 0.0, 0.0, eta );
    real64 dP_k  = EvaluateJacobiPolynomialDerivative( k, 2.0 * m + 2.0, 0.0, chi );

    real64 f = pow( 1.0 - z, m );
    // Gradient en x
    gradPsiX[0] = (1.0 / (1.0 - z)) * dP_i * P_j * f * P_k;

    // Gradient en y
    gradPsiX[1] = (1.0 / (1.0 - z)) * dP_j * P_i * f * P_k;

    // Gradient en z
    real64 dxi_dz  = x / std::pow( 1.0 - z, 2 );
    real64 deta_dz = y / std::pow( 1.0 - z, 2 );
    real64 dchi_dz = 2.0;
    real64 df_dz   = -m * std::pow( 1.0 - z, m - 1 );
    gradPsiX[2] =
      dP_i * dxi_dz * P_j * f * P_k +
      P_i * dP_j * deta_dz * f * P_k +
      P_i * P_j * df_dz * P_k +
      P_i * P_j * f * dP_k * dchi_dz;


  }

  /**
   * @brief Compute the inverse of the VanDerMonde matrix
   * @return The inverse of the VanDerMonde matrix
   */
  GEOS_FORCE_INLINE
  static array2d< real64 > computeVanderMondeMatrixInverse( )
  {
    array2d< real64 > VDM;
    VDM.resize( numNodes, numNodes );
    real64 PsiX[numNodes] = {};
    real64 coords[numNodes][3] = {};

    generatePointsCoordinates( coords );

    for( int j = 0; j < numNodes; ++j )
    {
      localIndex count = 0;
      generateIndexesHost( [&]( localIndex const p, localIndex const r, localIndex const s )
      {
        PsiX[count] = calcModal( p, r, s, coords[j] );
        VDM[count][j] = PsiX[count];
        ++count;
      } );
    }
    array2d< real64 > VDM_inv;
    VDM_inv.resize( numNodes, numNodes );
    // Inversion of VanDerMonde matrix
    BlasLapackLA::matrixInverse( VDM, VDM_inv );
    return VDM_inv;
  }

  /**
   * @brief Evaluate shape functions of a linear pyramid (5-node) at a quadrature point.
   * @param[in] q A quadrature point index.
   * @param[in] VDM_inv The inverse of the VanDerMonde matrix
   * @param[out] N Array to store shape function values (array of size numNodes)
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static constexpr void calcN( localIndex q, arrayView2d< real64 const > VDM_inv, real64 (& N)[numNodes] )
  {

    real64 PsiX[numNodes] = {0.0};
    real64 coords[numNodes][3] = {{0.0}};

    //Generate the coordinates of the support points
    generatePointsCoordinates( coords );


    localIndex count = 0;
    //Compute the modal functions at the quadrature point
    generateIndexes( [&]( localIndex const p, localIndex const r, localIndex const s )
    {
      PsiX[count] = calcModal( p, r, s, coords[q] );
      ++count;
    } );

    //Compute the nodal functions at the quadrature point
    for( int i = 0; i < numNodes; ++i )
    {
      N[i] = 0.0;
      for( int j = 0; j < numNodes; ++j )
      {
        N[i] += VDM_inv[i][j] * PsiX[j];
      }
    }

  }
  /**
   * @brief Evaluate shape functions of a linear pyramid (5-node) at a quadrature point.
   * @param[in] X Coordinates in reference pyramid (array of size 3)
   * @param[in] VDM_inv The inverse of the VanDerMonde matrix
   * @param[out] N Array to store shape function values (array of size numNodes)
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static constexpr void calcN( real64 const (&X)[3], arrayView2d< real64 const > VDM_inv, real64 (& N)[numNodes] )
  {

    real64 PsiX[numNodes] = {};

    localIndex count = 0;
    //Compute the modal functions at the given point
    generateIndexes( [&]( localIndex const p, localIndex const r, localIndex const s )
    {
      PsiX[count] = calcModal( p, r, s, X );
      ++count;
    } );

    //Compute the nodal functions at the given points
    for( int i = 0; i < numNodes; ++i )
    {
      N[i] = 0.0;
      for( int j = 0; j < numNodes; ++j )
      {
        N[i] += VDM_inv[i][j] * PsiX[j];
      }
    }

  }

  /**
   * @brief Evaluate shape functions of a linear pyramid (5-node) at a quadrature point.
   * @param[in] q A quadrature point index.
   * @param[in] VDM_inv The inverse of the VanDerMonde matrix
   * @param[out] gradN Array to store shape function derivatives (array of size numNodes x 3)
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static constexpr void calcGradN( localIndex q, arrayView2d< real64 const > VDM_inv, real64 (& gradN)[numNodes][3] )
  {

    real64 gradModal[numNodes][3] = {{}};
    real64 coords[numNodes][3] = {};

    //Generate the coordinates of the support points
    generatePointsCoordinates( coords );

    localIndex count = 0;
    //Compute the modal functions derivatives at the quadrature point
    generateIndexes( [&]( localIndex const p, localIndex const r, localIndex const s )
    {
      calcGradModal( p, r, s, coords[q], gradModal[count] );
      ++count;
    } );


    //Compute the nodal functions derivatives at the quadrature point
    for( int i = 0; i < numNodes; ++i )
    {
      gradN[i][0] = 0.0;
      gradN[i][1] = 0.0;
      gradN[i][2] = 0.0;
      for( int j = 0; j < numNodes; ++j )
      {
        gradN[i][0] += VDM_inv[i][j] * gradModal[j][0];
        gradN[i][1] += VDM_inv[i][j] * gradModal[j][1];
        gradN[i][2] += VDM_inv[i][j] * gradModal[j][2];
      }
    }

  }

  /**
   * @brief Evaluate shape functions of a linear pyramid (5-node) at a given point.
   * @param[in] X Coordinates in reference pyramid (array of size 3)
   * @param[in] VDM_inv The inverse of the VanDerMonde matrix
   * @param[out] gradN Array to store shape function values (array of size numNodes)
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static constexpr void calcGradN( real64 const (&X)[3], arrayView2d< real64 const > VDM_inv, real64 (& gradN)[numNodes][3] )
  {

    real64 gradModal[numNodes][3] = {{}};

    localIndex count = 0;
    //Compute the modal functions derivatives at the given point
    generateIndexes( [&]( localIndex const p, localIndex const r, localIndex const s )
    {
      calcGradModal( p, r, s, X, gradModal[count] );
      ++count;
    } );

    //Compute the nodal functions derivatives at the given points
    for( int i = 0; i < numNodes; ++i )
    {
      gradN[i][0] = 0.0;
      gradN[i][1] = 0.0;
      gradN[i][2] = 0.0;
      for( int j = 0; j < numNodes; ++j )
      {
        gradN[i][0] += VDM_inv[i][j] * gradModal[j][0];
        gradN[i][1] += VDM_inv[i][j] * gradModal[j][1];
        gradN[i][2] += VDM_inv[i][j] * gradModal[j][2];
      }
    }

  }

  /**
   * @brief Computes the mass term Mij in the mass matrix
   * @param[in] VDM_inv The inverse of the VanDerMonde matrix
   * @param[in] func A callable object (e.g., a lambda function) that takes three parameters (localIndex a, localIndex b, real64 val)
   *   where 'a' and 'b' are the indices of the mass term in the mass matrix, and 'val' is the computed value of the mass term.
   *   The function should have the signature: void func(localIndex a, localIndex b, real64 val)
   */
  template< typename FUNC >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static constexpr void computeMassTerm( arrayView2d< real64 const > VDM_inv, FUNC && func )
  {

    // Gauss-Legendre points and weights for the quadrature
    constexpr real64 GLeQuadraturePoints[3] = { -0.7745966692, 0.0, 0.7745966692 };
    constexpr real64 GLeQuadratureWeights[3] = { 0.5555555556, 0.8888888889, 0.5555555556 };
    // Gauss-Jacobi points and weights for the quadrature
    constexpr real64 GJQuadraturePoints[3] = { 0.07299, 0.34700, 0.70500 };
    constexpr real64 GJQuadratureWeights[3] = { 0.15714, 0.14625, 0.02995 };



    real64 N[numNodes] = {0.0};
    for( localIndex a = 0; a < numNodes; ++a )
    {
      for( localIndex b = 0; b < numNodes; ++b )
      {
        real64 val = 0.0;
        for( int i = 0; i < 3; ++i )
        {
          for( int j = 0; j < 3; ++j )
          {
            for( int k = 0; k < 3; ++k )
            {

              real64 xi = GLeQuadraturePoints[i];
              real64 eta = GLeQuadraturePoints[j];
              real64 chi = GJQuadraturePoints[k];

              real64 weight = GLeQuadratureWeights[i] * GLeQuadratureWeights[j] * GJQuadratureWeights[k];
              real64 x_i = (1-chi)*xi;
              real64 y_j = (1-chi)*eta;
              real64 z_k = chi;
              real64 X[3] = {x_i, y_j, z_k};

              calcN( X, VDM_inv, N );

              val += weight * N[a] * N[b];

            }
          }
        }
        func( a, b, val );
      }
    }
  }

  /**
   * @brief Computes the stifness term Kij in the mass matrix
   * @param[in] VDM_inv The inverse of the VanDerMonde matrix
   * @param[in] func A callable object (e.g., a lambda function) that takes three parameters (localIndex a, localIndex b, real64 val)
   *   where 'a' and 'b' are the indices of the mass term in the stifness matrix, and 'val' is the computed value of the mass term.
   *   The function should have the signature: void func(localIndex a, localIndex b, real64 val)
   */
  template< typename FUNC >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void computeStiffnessTerm( arrayView2d< real64 const > VDM_inv, FUNC && func )
  {
    // Gauss-Legendre points and weights for the quadrature
    constexpr real64 GLeQuadraturePoints[3] = { -0.7745966692, 0.0, 0.7745966692 };
    constexpr real64 GLeQuadratureWeights[3] = { 0.5555555556, 0.8888888889, 0.5555555556 };
    // Gauss-Jacobi points and weights for the quadrature
    constexpr real64 GJQuadraturePoints[3] = { 0.07299, 0.34700, 0.70500 };
    constexpr real64 GJQuadratureWeights[3] = { 0.15714, 0.14625, 0.02995 };



    real64 gradN[numNodes][3] = {{0.0}};

    for( localIndex a = 0; a < numNodes; ++a )
    {
      for( localIndex b = 0; b < numNodes; ++b )
      {
        real64 val = 0.0;
        for( int i = 0; i < 3; ++i )
        {
          for( int j = 0; j < 3; ++j )
          {
            for( int k = 0; k < 3; ++k )
            {
              real64 xi = GLeQuadraturePoints[i];
              real64 eta = GLeQuadraturePoints[j];
              real64 chi = GJQuadraturePoints[k];

              real64 weight = GLeQuadratureWeights[i] * GLeQuadratureWeights[j] * GJQuadratureWeights[k];
              real64 x_i = (1-chi)*xi;
              real64 y_j = (1-chi)*eta;
              real64 z_k = chi;
              real64 X[3] = {x_i, y_j, z_k};

              calcGradN( X, VDM_inv, gradN );
              real64 dot = gradN[a][0] * gradN[b][0]
                           + gradN[a][1] * gradN[b][1]
                           + gradN[a][2] * gradN[b][2];
              val += weight * dot;
            }
          }
        }
        func( a, b, val );
      }
    }
  }

};

/// @copydoc Pk_Pyramid_BCD_impl
template< int ORDER >
class Pk_Pyramid_BCD final : public Pk_Pyramid_BCD_impl< ORDER >, public FiniteElementBase
{
public:

  /// The Implementation type
  using ImplType = Pk_Pyramid_BCD_impl< ORDER >;

  /// The number of nodes per element.
  using ImplType::numNodes;

  /// The max number of support points per element.
  using ImplType::maxSupportPoints;

  /// The number of quadrature points per element.
  using ImplType::numQuadraturePoints;

  Pk_Pyramid_BCD():
    FiniteElementBase( numNodes,
                       maxSupportPoints,
                       numQuadraturePoints )
  {}

#ifdef __CUDACC__
  #pragma diag_push
  #pragma nv_diag_suppress 20012
#endif
  GEOS_HOST_DEVICE
  virtual ~Pk_Pyramid_BCD() override final = default;
#ifdef __CUDACC__
  #pragma diag_pop
#endif

  /**
   * @brief Get the device-compatible implementation type.
   * @return A pointer to the device-compatible implementation type.
   */
  ImplType * getImpl()
  {
    return static_cast< ImplType * >(this);
  }

  /**
   * @brief Get the device-compatible implementation type.
   * @return A pointer to the device-compatible implementation type.
   */
  ImplType const * getImpl() const
  {
    return static_cast< ImplType const * >(this);
  }

};

/// Pyramid element with BCD basis functions of order 1.
using P1_Pyramid_BCD = Pk_Pyramid_BCD< 1 >;

/// Pyramid element with BCD basis functions of order 1.
using P1_Pyramid_BCD_impl = Pk_Pyramid_BCD_impl< 1 >;


}
}

#endif // GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_PKPYRAMIDBCD_HPP_
