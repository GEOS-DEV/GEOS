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



template< int ORDER >
class Pk_Pyramid_BCD final : public FiniteElementBase
{
public:

  /// The order of the finite element.
  static constexpr int order = ORDER;

  /// The number of shape functions per element.
  constexpr static localIndex numNodes = ( ORDER + 1 ) * ( ORDER + 2 ) * ( 2 * ORDER + 3 ) / 6;

  /// The number of faces points per element.
  constexpr static localIndex numFaces = 5;

  /// The maximum number of support points per element.
  constexpr static localIndex maxSupportPoints = numNodes;

  /// The number of quadrature points per element.
  constexpr static localIndex numQuadraturePoints = numNodes;

  /// The number of modal points per element.
  constexpr static localIndex numModes = numNodes;

  /** @cond Doxygen_Suppress */
  USING_FINITEELEMENTBASE
  /** @endcond Doxygen_Suppress */

  virtual ~Pk_Pyramid_BCD() = default;

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
   * @param N An array to pass back the shape function values for each support
   *   point.
   */
  //GEOS_HOST_DEVICE
  //GEOS_FORCE_INLINE
  //static void calcN( localIndex const q,
  //                   real64 (& N)[numNodes] )
  //{
  //  for( int a=0; a < numNodes; ++a )
  //  {
  //    N[ a ] = 0;
  //  }
  //  N[ q ] = 1.0;
  //}


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
   * @param[out] PsiX A real to pass back the modal function value
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
    else if( z > 1.0 || z < 0.0 )
    {
      GEOS_ERROR( "Invalid z coordinate for pyramid shape function calculation." );
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
    else if( z > 1.0 || z < 0.0 )
    {
      GEOS_ERROR( "Invalid z coordinate for pyramid shape function calculation." );
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
   * @brief Evaluate shape functions of a linear pyramid (5-node) at point X.
   * @param[in] X Coordinates in reference pyramid (array of size 3)
   * @param[out] N Array to store shape function values (array of size numNodes)
   */
  // GEOS_HOST_DEVICE
  // GEOS_FORCE_INLINE
  // static void calcN(real64 const (&X)[numNodes[3], real64 (&N)[numNodes])
  // {

  //   real64 VDM[numNodes][numNodes] = {};
  //   real64 PsiX[numNodes] = {};

  //   for (int j = 0; j < numNodes; ++j)
  //   {
  //     const real64 (&Mj)[3] = NODES[j];
  //     for (int i = 0; i < numModes; ++i) {
  //       auto [p, q, r] = MODES[i];
  //       calcModal(p, q, r, Mj, VDM[i][j]);
  //     }
  //   }

  //   for (int i = 0; i < numNodes; ++i) {
  //     auto [p, q, r] = MODES[i];
  //     calcModal(p, q, r, X, PsiX[i]);
  //   }

  //   for (int i = 0; i < numNodes; ++i) {
  //     N[i] = 0.0;
  //     for (int j = 0; j < numNodes; ++j) {
  //       N[i] += VDM_inv[i][j] * PsiX[j]; // VDM_inv A CALCULER
  //     }
  //   }

  // }


  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static constexpr void computeVanderMondeMatrix( array2d< real64 > & VDM )
  {
    VDM.resize( numNodes, numNodes );
    real64 PsiX[numNodes] = {};
    real64 coords[numNodes][3] = {};

    generatePointsCoordinates( coords );

    for( int j = 0; j < numNodes; ++j )
    {
      localIndex count = 0;
      generateIndexes( [&]( localIndex const p, localIndex const r, localIndex const s )
      {
        PsiX[count]=calcModal( p, r, s, coords[j] );
        VDM[count][j] = PsiX[count];
        ++count;
      } );
      //   }
    }
  }

  /// Other possible version
  /**
   * @brief Evaluate shape functions of a linear pyramid (5-node) at a quadrature point.
   * @param[in] q A quadrature point index.
   * @param[out] N Array to store shape function values (array of size numNodes)
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static constexpr void calcN( localIndex q, real64 (& N)[numNodes] )
  {

    array2d< real64 > VDM;
    VDM.resize( numModes, numNodes );
    real64 PsiX[numNodes] = {0.0};
    real64 coords[numNodes][3] = {{0.0}};

    generatePointsCoordinates( coords );

    //   for (int j = 0; j < numNodes; ++j)
    //   {
    //     localIndex count = 0;
    //    // for (int i = 0; i < numModes; ++i)
    //    // {
    //       //auto [p, q, r] = MODES[i];
    //       generateIndexes([&](localIndex const p, localIndex const r, localIndex const s)
    //       {

    //         PsiX[count]=calcModal(p, r, s, coords[j]);
    //         VDM[count][j] = PsiX[count];
    //         ++count;
    //       } );
    //  //   }
    //   }

    computeVanderMondeMatrix( VDM );
    array2d< real64 > VDM_inv;
    VDM_inv.resize( numNodes, numNodes );
    // Inversion of VanDerMonde matrix
    BlasLapackLA::matrixInverse( VDM, VDM_inv );


    localIndex count = 0;
    generateIndexes( [&]( localIndex const p, localIndex const r, localIndex const s )
    {
      PsiX[count] = calcModal( p, r, s, coords[q] );
      ++count;
    } );



    // }

    //real64 VDM_inv[numNodes][numNodes] = {{0}};
    // array2d<real64> VDM_inv;
    // VDM_inv.resize(numNodes, numNodes);
    // Inversion of VanDerMonde matrix
    // BlasLapackLA::matrixInverse( VDM, VDM_inv );

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
   * @param[out] N Array to store shape function values (array of size numNodes)
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static constexpr void calcN( real64 const (&X)[3], real64 (& N)[numNodes] )
  {

    array2d< real64 > VDM;
    VDM.resize( numModes, numNodes );
    real64 PsiX[numNodes] = {};

    //generatePointsCoordinates(coords);

    computeVanderMondeMatrix( VDM );

    //for (int i = 0; i < numNodes; ++i)
    //{
    // localIndex count = 0;
    localIndex count = 0;
    generateIndexes( [&]( localIndex const p, localIndex const r, localIndex const s )
    {
      PsiX[count] = calcModal( p, r, s, X );
      ++count;
    } );

    // }

    //real64 VDM_inv[numNodes][numNodes] = {{0}};
    array2d< real64 > VDM_inv;
    VDM_inv.resize( numNodes, numNodes );
    // Inversion of VanDerMonde matrix
    BlasLapackLA::matrixInverse( VDM, VDM_inv );

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
   * @param[out] gradN Array to store shape function derivatives (array of size numNodes x 3)
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static constexpr void calcGradN( localIndex q, real64 (& gradN)[numNodes][3] )
  {

    array2d< real64 > VDM;
    VDM.resize( numModes, numNodes );
    real64 gradModal[numNodes][3] = {{}};
    real64 coords[numNodes][3] = {};

    generatePointsCoordinates( coords );
    computeVanderMondeMatrix( VDM );
    array2d< real64 > VDM_inv;
    VDM_inv.resize( numNodes, numNodes );
    // Inversion of VanDerMonde matrix
    BlasLapackLA::matrixInverse( VDM, VDM_inv );



    localIndex count = 0;
    generateIndexes( [&]( localIndex const p, localIndex const r, localIndex const s )
    {
      calcGradModal( p, r, s, coords[q], gradModal[count] );
      ++count;
    } );


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
   * @param[in] q A quadrature point index.
   * @param[out] N Array to store shape function values (array of size numNodes)
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static constexpr void calcGradN( real64 const (&X)[3], real64 (& gradN)[numNodes][3] )
  {

    array2d< real64 > VDM;
    VDM.resize( numModes, numNodes );
    real64 gradModal[numNodes][3] = {{}};
    computeVanderMondeMatrix( VDM );
    array2d< real64 > VDM_inv;
    VDM_inv.resize( numNodes, numNodes );
    // Inversion of VanDerMonde matrix
    BlasLapackLA::matrixInverse( VDM, VDM_inv );

    localIndex count = 0;
    generateIndexes( [&]( localIndex const p, localIndex const r, localIndex const s )
    {
      calcGradModal( p, r, s, X, gradModal[count] );
      ++count;
    } );


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
   * @brief Calculate modal base functions derivative values for a modal point (i,j,k) at a
   *   point X.
   * @param[in] i,j,k Indexes of a modal point.
   * @param[in] X Coordinates in reference pyramid (array of size 3)
   * @param[out] GradPsiX A real to pass back the modal function derivative values
   */
//   static void calcGradModal( int i, int j, int k,
//                              real64 const (&X)[3],
//                              real64 (&gradPsiX)[3] )
//   {
//     real64 x = X[0];
//     real64 y = X[1];
//     real64 z = X[2];

//     if (z == 1.0) {
//         if (z == 1.0 || z < 0.0 || z > 1.0) return Eigen::Vector3d::Zero();

//         real64 xi = x / (1.0 - z);
//         real64 eta = y / (1.0 - z);
//         real64 chi = 2.0 * z - 1.0;

//         int m = std::max(i, j);

//         real64 P_i   = boost::math::jacobi(i, 0.0, 0.0, xi);
//         real64 P_j   = boost::math::jacobi(j, 0.0, 0.0, eta);
//         real64 P_k   = boost::math::jacobi(k, 2.0 * m + 2.0, 0.0, chi);

//         real64 dP_i  = boost::math::jacobi_derivative(i, 0.0, 0.0, xi, 1);
//         real64 dP_j  = boost::math::jacobi_derivative(j, 0.0, 0.0, eta, 1);
//         real64 dP_k  = boost::math::jacobi_derivative(k, 2.0 * m + 2.0, 0.0, chi, 1);

//         real64 f = std::pow(1.0 - z, m);

//         // Gradient en x
//         gradPsiX[0] = (1.0 / (1.0 - z)) * dP_i * P_j * f * P_k;

//         // Gradient en y
//         gradPsiX[1] = (1.0 / (1.0 - z)) * dP_j * P_i * f * P_k;

//         // Gradient en z
//         real64 dxi_dz  = x / std::pow(1.0 - z, 2);
//         real64 deta_dz = y / std::pow(1.0 - z, 2);
//         real64 dchi_dz = 2.0;
//         real64 df_dz   = -m * std::pow(1.0 - z, m - 1);

//         gradPsiX[2] =
//             dP_i * dxi_dz * P_j * f * P_k +
//             P_i * dP_j * deta_dz * f * P_k +
//             P_i * P_j * df_dz * P_k +
//             P_i * P_j * f * dP_k * dchi_dz;
//     }

//   }

//     /**
//    * @brief Evaluate shape functions derivatives of a linear pyramid (5-node) at point X.
//    * @param[in] X Coordinates in reference pyramid (array of size 3)
//    * @param[out] gradN Array to store shape function derivative values (array of size (numNodes = 5) * 3)
//    */
//   GEOS_HOST_DEVICE
//   GEOS_FORCE_INLINE
//   static void calcGradN(real64 const (&X)[3], real64 (&gradN)[numNodes][3])
//   {

//     real64 VDM[numNodes][numNodes] = {};
//     real64 gradPsiX[numNodes][3] = {};

//     for (int j = 0; j < numNodes; ++j) {
//       const real64 (&Mj)[3] = NODES[j];
//     }

//     for (int i = 0; i < numModes; ++i) {
//        auto [p, q, r] = MODES[i];
//        calcModal(p, q, r, Mj, VDM[i][j]);
//     }

//     for (int i = 0; i < numNodes; ++i) {
//       auto [p, q, r] = MODES[i];
//       calcgradModal(p, q, r, X, gradPsiX[i]);
//     }

//     for (int i = 0; i < numNodes; ++i) {
//       N[i] = 0.0;
//       for (int j = 0; j < numNodes; ++j) {
//         N[i] += VDM_inv[i][j] * gradPsiX[j]; // VDM_inv A CALCULER
//       }
//     }

//   }

  /**
   * @brief Evaluate shape functions derivatives of a linear pyramid (5-node) at point X.
   * @param[in] X Coordinates in reference pyramid (array of size 3)
   * @param[out] gradN Array to store shape function derivative values (array of size (numNodes = 5) * 3)
   */
  // GEOS_HOST_DEVICE
  // GEOS_FORCE_INLINE
  // static void calcGradN(real64 (&gradN)[numNodes])
  // {

  //   constexpr real64 coords[numNodes][3] = generatePointsCoordinates();

  //   real64 VDM[numNodes][numNodes] = {};
  //   real64 gradPsiX[numNodes][3] = {};

  //   for (localIndex j =0 ; j < numNodes ; ++j)
  //   {
  //     for (int i = 0; i < numModes; ++i)
  //     {
  //        generatesIndexes([&](localIndex const p, localIndex const q, localIndex const r)
  //        {
  //          calcModal(p, q, r, coords[j], VDM[i][j]);
  //        });
  //     }
  //  }

  //   for (int i = 0; i < numNodes; ++i)
  //   {
  //     generatesIndexes([&](localIndex const p, localIndex const q, localIndex const r)
  //     {
  //       calcgradModal(p, q, r, coords[i], gradPsiX[i]);
  //     } );
  //   }

  //   // Inversion of VanDerMonde matrix
  //   real64 VDM_inv[numNodes][numNodes] = {};
  //   LvArray::tensorOps::invert< 3 >( VDM_inv, VDM );


  //   //TODO: check if this is correct
  //   for (int i = 0; i < numNodes; ++i) {
  //     gradN[i] = 0.0;
  //     for (int j = 0; j < numNodes; ++j) {
  //       gradN[i][0] += VDM[i][j] * gradPsiX[j][0];
  //       gradN[i][1] += VDM[i][j] * gradPsiX[j][1];
  //       gradN[i][2] += VDM[i][j] * gradPsiX[j][2];
  //     }
  //   }

  // }


  /**
   * @brief Computes the mass term Mij in the mass matrix
   * @param[in] i,j Coordinates of the mass term in the mass matrix
   * @param[in] N Array containing the shape function values at the quadrature points
   * @param[out] Mij real to store the value of the mass term
   */
  template< typename FUNC >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static constexpr void computeMassTerm( FUNC && func )
  {

    //real64 Mij = 0;
    // Gauss-Legendre points and weights for the quadrature
    constexpr real64 GLeQuadraturePoints[3] = { -0.7745966692, 0.0, 0.7745966692 };
    constexpr real64 GLeQuadratureWeights[3] = { 0.5555555556, 0.8888888889, 0.5555555556 };
    // Gauss-Jacobi points and weights for the quadrature
    constexpr  real64 GJQuadraturePoints[3] = { 0.07299, 0.34700, 0.70500 };
    constexpr  real64 GJQuadratureWeights[3] = { 0.15714, 0.14625, 0.02995 };



    real64 N[numNodes] = {0.0};   // Initialize the shape function values array

    for( localIndex a = 0; a < numNodes; ++a )
    {
      for( localIndex b = 0; b < numNodes; ++b )
      {
        real64 val = 0.0;   // Initialize the value to accumulate
        for( int i = 0; i < 3; ++i )
        {
          for( int j = 0; j < 3; ++j )
          {
            for( int k = 0; k < 3; ++k )
            {
              //real64 xi =  QL[i];
              //real64 eta = QL[j];
              //real64 chi = QJ[k];

              real64 xi = GLeQuadraturePoints[i];
              real64 eta = GLeQuadraturePoints[j];
              real64 chi = GJQuadraturePoints[k];

              real64 weight = GLeQuadratureWeights[i] * GLeQuadratureWeights[j] * GJQuadratureWeights[k];
              //real64 weight = WL[i] * WL[j] * WJ[k];
              real64 x_i = (1-chi)*xi;
              real64 y_j = (1-chi)*eta;
              real64 z_k = chi;
              real64 X[3] = {x_i, y_j, z_k};

              calcN( X, N );
              val += weight * N[a] * N[b];

            }
          }
        }
        func( a, b, val ); // Call the function with the computed value
      }
    }
  }

//     /**
//    * @brief Computes the mass term Mij in the mass matrix
//    * @param[in] i,j Coordinates of the mass term in the mass matrix
//    * @param[in] N Array containing the shape function values at the quadrature points
//    * @param[out] Mij real to store the value of the mass term
//    */
//   GEOS_HOST_DEVICE
//   GEOS_FORCE_INLINE
//   static void computeMassTerm(int const i, int const j, real64 const (&N)[numNodes], real64 & Mij)
//   {
//       real64 QL[3] = {-0.7745966692,0.0,0.7745966692};
//       real64 WL[3] = {0.5555555556,0.8888888889,0.5555555556};
//       real64 QJ[3] = {0.07299,0.34700,0.70500};
//       real64 WJ[3] = {0.15714,0.14625,0.02995};

//       real64 Mij = 0;

//       for(int i = 0; i < 3; ++i) {
//           for(int j = 0; j < 3; ++j) {
//               for(int k = 0; k < 3; ++k) {
//                   real64 xi = QL[i];
//                   real64 eta = QL[j];
//                   real64 chi = QJ[k];

//                   real64 weight = WL[i] * WL[j] * WJ[k];

//                   real64 x_i = (1-chi)*xi;
//                   real64 y_j = (1-chi)*eta;
//                   real64 z_k = chi;

//                   real64 X[3] = {x_i, y_j, z_k};
//                   calcN(X,N);
//                   Mij+= weight * N[a] * N[b] ;

//               }
//           }
//       }
//   }



//   /**
//    * @brief Computes the stifness term Kij in the stifness matrix
//    * @param[in] i,j Coordinates of the stifness term in the stifness matrix
//    * @param[in] gradN Array containing the shape function derivative values at the quadrature points
//    * @param[out] Kij real to store the value of the stifness term
//    */
//   GEOS_HOST_DEVICE
//   GEOS_FORCE_INLINE
//   static void computeStifnessTerm(int const i, int const j, real64 const (&gradN)[numNodes][3], real64 & Kij)
//   {
//       real64 QL[3] = {-0.7745966692,0.0,0.7745966692};
//       real64 WL[3] = {0.5555555556,0.8888888889,0.5555555556};
//       real64 QJ[3] = {0.07299,0.34700,0.70500};
//       real64 WJ[3] = {0.15714,0.14625,0.02995};

//       real64 Kij = 0;

//       for(int i = 0; i < 3; ++i) {
//           for(int j = 0; j < 3; ++j) {
//               for(int k = 0; k < 3; ++k) {
//                   real64 xi = QL[i];
//                   real64 eta = QL[j];
//                   real64 chi = QJ[k];

//                   real64 weight = WL[i] * WL[j] * WJ[k];

//                   real64 x_i = (1-chi)*xi;
//                   real64 y_j = (1-chi)*eta;
//                   real64 z_k = chi;

//                   real64 X[3] = {x_i, y_j, z_k};
//                   calcGradN(X, gradN);

//                   Kij+= weight * gradN[a] * gradN[b] ; // PDT SCALAIRE A CODER
//               }
//           }
//       }
//   }

//   GEOS_HOST_DEVICE
//   GEOS_FORCE_INLINE
//   localIndex linearIndexForBase(localIndex i, localIndex j){
//      return i + (order + 1) * j;
//   }


//   /**
//    * @brief Computes the flow term Rij in the mass matrix
//    * @param[in] i,j Coordinates of the mass term in the mass matrix
//    * @param[in] N Array containing the shape function values at the quadrature points
//    * @param[out] Mij real to store the value of the mass term
//    */
//   GEOS_HOST_DEVICE
//   GEOS_FORCE_INLINE
//   static void computeBaseInterfaceTerm(real64 const (&R)[(order+1)*(order+1)]){

//     localIndex I;
//     std::vector<real64> x1D, w1D;
//     LagrangeBasis1D::GaussLobatto1D::getPointsAndWeights(order, x1D, w1D);

//     for(int i=0; i<ordre+1; ++i){
//         for(int j=0; j<ordre+1; ++j){
//             I = linearIndexForBase(i,j);
//             R[I]+= w1D[i] * w1D[j];
//         }
//     }
//   }

//   /**
//    * @brief Computes a! / ( b! * c! ) with b + c >= a >= b, c
//    * @param[in] a
//    * @param[in] b
//    * @param[in] c
//    * @return a!/(b!*c!)
//    */
//   GEOS_HOST_DEVICE
//   GEOS_FORCE_INLINE
//   constexpr static real64 integralTerm( const int a, const int b, const int c )
//   {
//     real64 res = 1.0;
//     int num = a;
//     int den = c;
//     for( int i = b; i > 0; i-- )
//     {
//       res *= ( (real64) num ) /  i;
//       num--;
//     }
//     for( int i = num; i > 0; i-- )
//     {
//       res *= ( (real64) i ) /  den;
//       den--;
//     }
//     for( int i = den; i > 0; i-- )
//     {
//       res /= i;
//     }

//     return res;
//   }

//   /**
//    * @brief Computes the correction factor for the superposition integral in case a function has been derived once.
//    *   The indices of the original function are ii1, j1, k1 and l1, and those of the once-derived one are i2, j2, k2 and l2
//    * @param[in] i1
//    * @param[in] j1
//    * @param[in] k1
//    * @param[in] l1
//    * @param[in] i2
//    * @param[in] j2
//    * @param[in] k2
//    * @param[in] l2
//    * @param[in] dim the dimension, 2 or 3
//    * @return the correction factor to be applied to the superposition integral
//    */
//   GEOS_HOST_DEVICE
//   GEOS_FORCE_INLINE
//   static constexpr real64 correctionFactorDerivative( int const i1, int const j1, int const k1, int const l1, int const i2, int const j2,
// int const k2, int const l2, int const dim )
//   {
//     return (i1+j1+k1+l1+dim)* (i2==0 ? 1 : i2) * (j2==0 ? 1 : j2)* (k2==0 ? 1 : k2)* (l2==0 ? 1 : l2);
//   }


//   /**
//    * @brief Computes the superposition integral between Bernstein-Bézier functions indexed by
//    *  (i1, j1, k1, l1) and (i1, j2, k2, l2)
//    * @param[in] i1
//    * @param[in] j1
//    * @param[in] k1
//    * @param[in] l1
//    * @param[in] i2
//    * @param[in] j2
//    * @param[in] k2
//    * @param[in] l2
//    * @return the superposition integral over the barycentric coordinates
//    */
//   GEOS_HOST_DEVICE
//   GEOS_FORCE_INLINE
//   static constexpr real64 computeSuperpositionIntegral( const int i1, const int j1, const int k1, const int l1,
//                                                         const int i2, const int j2, const int k2, const int l2 )
//   {
//     return (integralTerm( i1+i2, i1, i2 )*
//             integralTerm( j1+j2, j1, j2 )*
//             integralTerm( k1+k2, k1, k2 )*
//             integralTerm( l1+l2, l1, l2 ))/
//            integralTerm( i1+j1+k1+l1+i2+j2+k2+l2+3,
//                          i1+j1+k1+l1+3, i2+j2+k2+l2+3 );
//   }

//   /**
//    * @brief Computes the superposition integral over a face between Bernstein-Bézier functions whose indices are given by
//    *  (i1, j1, k1, l1=0 ) and (i1, j2, k2, l2=0)
//    * @param[in] i1
//    * @param[in] j1
//    * @param[in] k1
//    * @param[in] i2
//    * @param[in] j2
//    * @param[in] k2
//    * @return the superposition integral over the barycentric coordinates
//    */
//   GEOS_HOST_DEVICE
//   GEOS_FORCE_INLINE
//   constexpr static real64 computeFaceSuperpositionIntegral( const int i1, const int j1, const int k1,
//                                                             const int i2, const int j2, const int k2 )
//   {
//     return ((i1+k1+j1+3)*(i2+j2+k2+3)*
//             integralTerm( i1+i2, i1, i2 )*
//             integralTerm( j1+j2, j1, j2 )*
//             integralTerm( k1+k2, k1, k2 ))/
//            integralTerm( i1+j1+k1+i2+j2+k2+2,
//                          i1+j1+k1+2, i2+j2+k2+2 );
//   }


//   /**
//    * @brief Helper function for static for loop
//    * @tparam FUNC the callback function
//    * @tparam ...Is integer indices of the loop
//    * @param func the callback function to call for each index
//    * @param Is the integer indices of the loop
//    */
//   template< typename FUNC, int... Is >
//   GEOS_HOST_DEVICE
//   GEOS_FORCE_INLINE
//   static constexpr void loop( FUNC && func, std::integer_sequence< int, Is... > )
//   {
//     ( func( std::integral_constant< int, Is >{} ), ... );
//   }

//   /**
//    * @brief Helper function for loop over barycentric coordinates
//    * @tparam FUNC the callback function
//    * @param func the callback function to call for each index
//    */
//   template< typename FUNC >
//   GEOS_HOST_DEVICE
//   GEOS_FORCE_INLINE
//   static constexpr void barycentricCoordinateLoop( FUNC && func )
//   {
//     loop( [&func] ( auto const i ) {
//       func( std::integral_constant< int, i >{} );
//     }, std::make_integer_sequence< int, 4 >{} );
//   }

//   /**
//    * @brief Helper function for loop over tet basis functions
//    * @tparam FUNC the callback function
//    * @param func the callback function to call for each index
//    */
//   template< typename FUNC >
//   GEOS_HOST_DEVICE
//   GEOS_FORCE_INLINE
//   static constexpr void basisLoop( FUNC && func )
//   {
//     loop( [&func] ( auto const i )
//     {
//       constexpr int i1 = ORDER - i;
//       loop( [&func] ( auto const j )
//       {
//         constexpr int j1 = ORDER - j;
//         if constexpr ( j1 <= ORDER - i1 )
//         {
//           loop( [&func] ( auto const k )
//           {
//             constexpr int k1 = ORDER - k;
//             if constexpr ( k1 <= ORDER - i1 - j1 )
//             {
//               constexpr int l1 = ORDER - i1 - j1 - k1;
//               constexpr int c1 = dofIndex< i1, j1, k1 >();
//               func( std::integral_constant< int, c1 >{},
//                     std::integral_constant< int, i1 >{},
//                     std::integral_constant< int, j1 >{},
//                     std::integral_constant< int, k1 >{},
//                     std::integral_constant< int, l1 >{} );
//             }
//           }, std::make_integer_sequence< int, ORDER + 1 > {} );
//         }
//       }, std::make_integer_sequence< int, ORDER + 1 > {} );
//     }, std::make_integer_sequence< int, ORDER + 1 > {} );
//   }

//   /**
//    * @brief Helper function for loop over tet basis functions that have one index in a given set of indices.
//    *   If multiple indices are in the given list, the callback is called multiple times.
//    * @tparam FUNC the callback function
//    * @tparam Is the setindices
//    */
//   template< int... Is, typename FUNC >
//   GEOS_HOST_DEVICE
//   GEOS_FORCE_INLINE
//   static constexpr void conditionalBasisLoop( FUNC && func )
//   {
//     loop( [&func] ( auto const i )
//     {
//       constexpr int i1 = ORDER -  i;
//       loop( [&func, i1] ( auto const j )
//       {
//         constexpr int j1 = ORDER - j;
//         constexpr int ii1 = ORDER -i1;
//         if constexpr ( j1 <= ii1 )
//         {
//           loop( [&func, i1, j1] ( auto const k )
//           {
//             constexpr int k1 = ORDER - k;
//             constexpr int ji1 = ORDER - i1 - j1;
//             if constexpr ( k1 <= ji1 )
//             {
//               constexpr int l1 = ORDER - i1 - j1 - k1;
//               constexpr int c1 = dofIndex< i1, j1, k1 >();
//               ( ( (i1 == Is) &&
//                   ( void( func(
//                             std::integral_constant< int, 0 >{},
//                             std::integral_constant< int, i1 >{},
//                             std::integral_constant< int, c1 >{},
//                             std::integral_constant< int, j1 >{},
//                             std::integral_constant< int, k1 >{},
//                             std::integral_constant< int, l1 >{} ) ), 1 ) ) || ... );
//               ( ( (j1 == Is) &&
//                   ( void( func(
//                             std::integral_constant< int, 1 >{},
//                             std::integral_constant< int, j1 >{},
//                             std::integral_constant< int, c1 >{},
//                             std::integral_constant< int, i1 >{},
//                             std::integral_constant< int, k1 >{},
//                             std::integral_constant< int, l1 >{} ) ), 1 ) ) || ... );
//               ( ( (k1 == Is) &&
//                   ( void( func(
//                             std::integral_constant< int, 2 >{},
//                             std::integral_constant< int, k1 >{},
//                             std::integral_constant< int, c1 >{},
//                             std::integral_constant< int, i1 >{},
//                             std::integral_constant< int, j1 >{},
//                             std::integral_constant< int, l1 >{} ) ), 1 ) ) || ... );
//               ( ( (l1 == Is) &&
//                   ( void( func(
//                             std::integral_constant< int, 3 >{},
//                             std::integral_constant< int, l1 >{},
//                             std::integral_constant< int, c1 >{},
//                             std::integral_constant< int, i1 >{},
//                             std::integral_constant< int, j1 >{},
//                             std::integral_constant< int, k1 >{} ) ), 1 ) ) || ...);
//             }
//           }, std::make_integer_sequence< int, ORDER + 1 > {} );
//         }
//       }, std::make_integer_sequence< int, ORDER + 1 > {} );
//     }, std::make_integer_sequence< int, ORDER + 1 > {} );
//   }

//   /**
//    * @brief Helper function for loop over barycentric coordinates of a face.
//    * @tparam FUNC the callback function
//    */
//   template< typename FUNC >
//   GEOS_HOST_DEVICE
//   GEOS_FORCE_INLINE
//   static constexpr void faceBarycentricCoordinateLoop( FUNC && func )
//   {
//     loop( [&func] ( auto const i ) {
//       func( std::integral_constant< int, i >{} );
//     }, std::make_integer_sequence< int, 3 >{} );
//   }

//   /**
//    * @brief Computes the reference mass matrix, i.e., the superposition matrix of the shape functions
//    *   in barycentric coordinates. The real-world mass matrix can be obtained by using the multiplying
//    *   this matrix by the determinant of the Jacobian.
//    *
//    * @param[out] m The mass matrix
//    */
//   GEOS_HOST_DEVICE
//   GEOS_FORCE_INLINE
//   static
//   void
//   computeReferenceMassMatrix( arraySlice2d< real64 > const & m )
//   {
//     basisLoop( [ &m ] ( auto const cc1, auto const ii1, auto const jj1, auto const kk1, auto const ll1 )
//     {
//       constexpr int c1 = cc1;
//       constexpr int i1 = ii1;
//       constexpr int j1 = jj1;
//       constexpr int k1 = kk1;
//       constexpr int l1 = ll1;
//       // Needed for compilers that do not support constexpr lambdas
//       GEOS_UNUSED_VAR( c1, i1, j1, k1, l1 );
//       basisLoop( [ &m ] ( auto const c2, auto const i2, auto const j2, auto const k2, auto const l2 )
//       {
//         constexpr real64 val = computeSuperpositionIntegral( i1, j1, k1, l1, i2, j2, k2, l2 );
//         m[ c1 ][ c2 ] = val;
//       } );
//     } );
//   }

//   /**
//    * @brief Computes the reference damping matrix, i.e., the superposition matrix of the shape functions
//    *   in barycentric coordinates over faces. The real-world mass matrix can be obtained by using the multiplying
//    *   this matrix by the determinant of the face Jacobian.
//    *
//    * @param[out] d The damping matrix
//    * @param[in] face1Damped Whether the first face contributes to the damping term
//    * @param[in] face2Damped Whether the second face contributes to the damping term
//    * @param[in] face3Damped Whether the third face contributes to the damping term
//    * @param[in] face4Damped Whether the fourth face contributes to the damping term
//    */
//   template< typename DAMPING >
//   GEOS_HOST_DEVICE
//   GEOS_FORCE_INLINE
//   static
//   void
//   computeReferenceDampingMatrix( real64 (& d)[numNodes][numNodes], bool const face1Damped, bool const face2Damped, bool const
// face3Damped, bool const face4Damped )
//   {

//     for( int c1 = 0; c1 < numNodes; c1++ )
//     {
//       for( int c2 = 0; c2 < numNodes; c2++ )
//       {
//         d[ c1 ][ c2 ] = 0;
//       }
//     }
//     conditionalBasisLoop< 0 >( [&] ( auto const f1, auto const, auto const c1, auto const i1, auto const j1, auto const k1 )
//     {
//       conditionalBasisLoop< 0 >( [&] ( auto const f2, auto const, auto const c2, auto const i2, auto const j2, auto const k2 )
//       {
//         if constexpr ( f1 == f2 )
//         {
//           constexpr real64 val = computeFaceSuperpositionIntegral( i1, j1, k1, i2, j2, k2 );
//           if( ( f1 == 0 && face1Damped ) ||
//               ( f1 == 1 && face2Damped ) ||
//               ( f1 == 2 && face3Damped ) ||
//               ( f1 == 3 && face4Damped ) )
//           {
//             d[ c1 ][ c2 ] += val;
//           }
//         }
//       } );
//     } );
//   }

//   /**
//    * @brief Computes the non-zero contributions inside the element of the
//    *   mass matrix M, i.e., the superposition matrix of shape functions
//    * @param X Array containing the coordinates of the support points.
//    * @param func Callback function accepting three parameters: i, j (local d.o.f. inside the element) and M_ij
//    */
//   template< typename FUNC >
//   GEOS_HOST_DEVICE
//   GEOS_FORCE_INLINE
//   static
//   constexpr
//   void
//   computeMassTerm( real64 const (&X)[4][3],
//                    FUNC && func )
//   {
//     real64 detJ = LvArray::math::abs( jacobianDeterminant( X ));
//     basisLoop( [&func, &detJ] ( auto const cc1, auto const ii1, auto const jj1, auto const kk1, auto const ll1 )
//     {
//       constexpr int c1 = cc1;
//       constexpr int i1 = ii1;
//       constexpr int j1 = jj1;
//       constexpr int k1 = kk1;
//       constexpr int l1 = ll1;
//       //Needed for compilors that do not support constexpr lambdas
//       GEOS_UNUSED_VAR( c1, i1, j1, k1, l1 );
//       basisLoop( [&func, &detJ] ( auto const c2, auto const i2, auto const j2, auto const k2, auto const l2 )
//       {
//         constexpr real64 val = computeSuperpositionIntegral( i1, j1, k1, l1, i2, j2, k2, l2 );
//         func( c1, c2, val * detJ );
//       } );
//     } );
//   }

//   /**
//    * @brief Function to compute the factor for the flux derivative term
//    * @param X Array containing the coordinates of the support points.
//    * @param x1 Index of the first edge vertex
//    * @param x2 Index of the second edge vertex
//    * @param o1 Index of the first face vertex
//    * @param o2 Index of the second face vertex
//    */
//   GEOS_HOST_DEVICE
//   GEOS_FORCE_INLINE
//   static constexpr real64
//   computeFluxDerivativeFactor( real64 const (&X)[4][3], int x1, int x2, int o1, int o2 )    // Order x1, x2, o1, o2
//   {

//     real64 detJ = LvArray::math::abs( jacobianDeterminant( X ));

//     // height with respect to o2

//     real64 detJf1 = LvArray::math::abs( faceJacobianDeterminant( o1, X ));

//     real64 scal = 1.0;

//     if( o1 != o2 )
//     {

//       // scalar product

//       real64 detJf2 = LvArray::math::abs( faceJacobianDeterminant( o2, X ));

//       real64 el2[3][2] = { { edgeLength2( x1, x2, X ), edgeLength2( o1, o2, X ) },

//         { edgeLength2( x1, o1, X ), edgeLength2( x2, o2, X ) },

//         { edgeLength2( x1, o2, X ), edgeLength2( x2, o1, X ) } };

//       real64 h2 = (el2[1][0]+el2[1][1])-(el2[2][0]+el2[2][1]);

//       h2 = (4.0 * el2[0][0]*el2[0][1] - h2 * h2)/16.0;

//       scal = (4.0*h2-detJf1 * detJf1 - detJf2 * detJf2 ) / (2.0 * detJf1 * detJf2);

//     }

//     return scal * detJf1 / detJ;

//   }

//   /**
//    * @brief Computes the non-zero contributions inside the element of the
//    *   stiffness matrix R, i.e., the superposition matrix of first derivatives of shape functions
//    * @param X Array containing the coordinates of the support points.
//    * @param func Callback function accepting three parameters: i, j (local d.o.f. inside the element) and R_ij
//    */
//   template< typename FUNC >
//   GEOS_HOST_DEVICE
//   GEOS_FORCE_INLINE
//   static
//   constexpr
//   void
//   computeStiffnessTerm( real64 const (&X)[4][3],
//                         FUNC && func )
//   {
//     real64 detJ = LvArray::math::abs( jacobianDeterminant( X ));
//     real64 dLambdadX[4][3] = {};
//     for( int j = 0; j < 3; j++ )
//     {
//       dLambdadX[1][j] =
//         ( ( X[ 2 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ]) * ( X[ 3 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) - ( X[ 3 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ])
// * ( X[ 2 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) ) / detJ;
//       dLambdadX[2][j] =
//         ( ( X[ 3 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ]) * ( X[ 1 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) - ( X[ 1 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ])
// * ( X[ 3 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) ) / detJ;
//       dLambdadX[3][j] =
//         ( ( X[ 1 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ]) * ( X[ 2 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) - ( X[ 2 ][ (j+1)%3 ] - X[ 0 ][ (j+1)%3 ])
// * ( X[ 1 ][ (j+2)%3 ] - X[ 0 ][ (j+2)%3 ] ) ) / detJ;
//       dLambdadX[0][j] = -dLambdadX[1][j] - dLambdadX[2][j] - dLambdadX[3][j];
//     }
//     basisLoop( [&func, &dLambdadX, &detJ] ( auto const cc1, auto const ci1, auto const cj1, auto const ck1, auto const cl1 )
//     {

//       constexpr int c1 = cc1;
//       constexpr int i1 = ci1;
//       constexpr int j1 = cj1;
//       constexpr int k1 = ck1;
//       constexpr int l1 = cl1;
//       //Not used in some combinations, but needed for constexpr
//       GEOS_UNUSED_VAR( c1, i1, j1, k1, l1 );
//       basisLoop( [&func, &dLambdadX, &detJ] ( auto const cc2, auto const ci2, auto const cj2, auto const ck2, auto const cl2 )
//       {
//         constexpr int c2 = cc2;
//         constexpr int i2 = ci2;
//         constexpr int j2 = cj2;
//         constexpr int k2 = ck2;
//         constexpr int l2 = cl2;
//         //Not used in some combinations, but needed for constexpr
//         GEOS_UNUSED_VAR( c2, i2, j2, k2, l2 );
//         barycentricCoordinateLoop( [&func, &dLambdadX, &detJ] ( auto const cd1 )
//         {
//           constexpr int d1 = cd1;
//           barycentricCoordinateLoop( [&func, &dLambdadX, &detJ] ( auto const d2 )
//           {
//             constexpr int ii1 = i1 + ( d1 == 0 ) * ( -1 );
//             constexpr int ij1 = j1 + ( d1 == 1 ) * ( -1 );
//             constexpr int ik1 = k1 + ( d1 == 2 ) * ( -1 );
//             constexpr int il1 = l1 + ( d1 == 3 ) * ( -1 );
//             constexpr int ii2 = i2 + ( d2 == 0 ) * ( -1 );
//             constexpr int ij2 = j2 + ( d2 == 1 ) * ( -1 );
//             constexpr int ik2 = k2 + ( d2 == 2 ) * ( -1 );
//             constexpr int il2 = l2 + ( d2 == 3 ) * ( -1 );
//             constexpr real64 factor1 = correctionFactorDerivative( i1, j1, k1, l1, ii1, ij1, ik1, il1, 3 );
//             constexpr real64 factor2 = correctionFactorDerivative( i2, j2, k2, l2, ii2, ij2, ik2, il2, 3 );
//             if constexpr (ii1 >= 0 && ij1 >= 0 && ik1 >= 0 && il1 >= 0 &&
//                           ii2 >= 0 && ij2 >= 0 && ik2 >= 0 && il2 >= 0)
//             {
//               constexpr real64 val = computeSuperpositionIntegral( ii1, ij1, ik1, il1, ii2, ij2, ik2, il2 ) * factor1 * factor2;
//               func( c1, c2, val * detJ * ( dLambdadX[d1][0]*dLambdadX[d2][0] + dLambdadX[d1][1]*dLambdadX[d2][1] +
// dLambdadX[d1][2]*dLambdadX[d2][2] ) );
//             }
//           } );
//         } );
//       } );
//     } );
//   }

//   /**
//    * Compute th length of the edge between two vertices i1 and i2
//    * @param i1 Index of the first vertex
//    * @param i2 Index of the second vertex
//    * @param X Array containing the coordinates of the support points.
//    * @return The squared length of the edge between the two vertices
//    */
//   GEOS_HOST_DEVICE
//   GEOS_FORCE_INLINE
//   static real64 edgeLength2( localIndex i1,
//                              localIndex i2,
//                              real64 const (&X)[4][3] )
//   {
//     real64 ab[3] = { X[ i2 ][ 0 ] - X[ i1 ][ 0 ],
//                      X[ i2 ][ 1 ] - X[ i1 ][ 1 ],
//                      X[ i2 ][ 2 ] - X[ i1 ][ 2 ] };
//     return ab[0] * ab[0] + ab[1] * ab[1]+ ab[2]*ab[2];
//   }

//   /**
//    * @brief Computes the non-zero contributions inside the element of the surface terms, including the value of
//    *   the superposition integral of basis functions (used for the penalty and damping terms) and
//    *   the superposition integral of the derivative of a function with the value of the other, used for the flux terms.
//    * @param X Array containing the coordinates of the support points.
//    * @param funcP Callback function for non-zero penalty-type terms, accepting seven parameters:
//    *   c1, c2 (local d.o.f. inside the element), f (index of the face, i.e., index of the opposite vertex for this element),
//    *   i1, j1 and k1 (local indices for the first shape function) and value
//    *   i2, j2 and k2 (local indices for the second shape function) and value
//    * @param funcF Callback function for non-zero flux-type terms, accepting four parameters:
//    *   c2, c2 (local d.o.f. inside the element), f (index of the face, i.e., index of the opposite vertex for this element), and value
//    *   i1, j1 and k1 (local indices for the first shape function) and value
//    *   i2, j2 and k2 (local indices for the second shape function) and value
//    */
//   template< typename FUNCP, typename FUNCF >
//   GEOS_HOST_DEVICE
//   GEOS_FORCE_INLINE
//   static
//   constexpr
//   void
//   computeSurfaceTerms( real64 const (&X)[4][3],
//                        FUNCP && funcP,
//                        FUNCF && funcF )
//   {
//     real64 detJf[4] = { faceJacobianDeterminant( 0, X ), faceJacobianDeterminant( 1, X ),
//                         faceJacobianDeterminant( 2, X ), faceJacobianDeterminant( 3, X ) };
//     conditionalBasisLoop< 0, 1 >( [&funcP, &funcF, &detJf]  ( auto const cf1, auto const cd, auto const cc1, auto const ci1, auto const
// cj1, auto const ck1 )
//     {
//       constexpr int f1 = cf1;
//       constexpr int d1 = cd;
//       constexpr int c1 = cc1;
//       constexpr int i1 = ci1;
//       constexpr int j1 = cj1;
//       constexpr int k1 = ck1;
//       // Not used in some combinations, but needed for constexpr
//       GEOS_UNUSED_VAR( c1, i1, j1, k1 );
//       conditionalBasisLoop< 0 >( [&funcP, &funcF, &detJf]  ( auto const cf2, auto const, auto const cc2, auto const ci2, auto const cj2,
// auto const ck2 )
//       {
//         constexpr int f2 = cf2;
//         constexpr int c2 = cc2;
//         constexpr int i2 = ci2;
//         constexpr int j2 = cj2;
//         constexpr int k2 = ck2;
//         // Not used in some combinations, but needed for constexpr
//         GEOS_UNUSED_VAR( c2, i2, j2, k2 );

//         // The second function is nonzero on the face indexed by f2, so we integrate on this face.
//         if constexpr ( f1 == f2 )
//         {
//           // compute penalty term iff the other function is also nonzero on the same face (i.e., d1==0)
//           if constexpr ( d1 == 0 )
//           {
//             constexpr real64 val = computeFaceSuperpositionIntegral( i1, j1, k1, i2, j2, k2 );
//             funcP( c1, c2, f2, i1, j1, k1, i2, j2, k2, val * detJf[ f2 ] );
//           }
//           // Compute flux term. This is nonzero in two cases.
//           // first case: function has exponent 1 wrt to the same face. In this case, one can derive it once wrt to the
//           // corresponding lambda and it will obtain a nonzero function on the face.
//           if constexpr ( d1 == 1 )
//           {
//             constexpr real64 derFactor = ( i1 + j1 + k1 + 4 );
//             constexpr real64 val = computeFaceSuperpositionIntegral( i1, j1, k1, i2, j2, k2 ) * derFactor;
//             funcF( c1, c2, f2, -1, i1, j1, k1, i2, j2, k2, val * detJf[ f2 ] );
//           }
//           // second case: function has exponent zero wrt f2.
//           // In this case, one can derive it wrt to any other face.
//           else if constexpr ( d1 == 0 )
//           {
//             faceBarycentricCoordinateLoop( [ &funcF, &detJf ]( auto const cl )
//             {
//               constexpr int l = cl;
//               constexpr int ii1 = i1 + ( l == 0 ) * ( -1 );
//               constexpr int ij1 = j1 + ( l == 1 ) * ( -1 );
//               constexpr int ik1 = k1 + ( l == 2 ) * ( -1 );
//               if constexpr (ii1 >= 0 && ij1 >= 0 && ik1 >= 0)
//               {
//                 constexpr real64 derFactor = ( ii1 + ij1 + ik1 + 4 );
//                 constexpr real64 val = computeFaceSuperpositionIntegral( ii1, ij1, ik1, i2, j2, k2 ) * derFactor;
//                 constexpr int f = l >= f2 ? l + 1 : l;
//                 funcF( c1, c2, f2, f, i1, j1, k1, i2, j2, k2, val * detJf[f2] );
//               }
//             } );
//           }
//         }
//       } );
//     } );
//   }

};

// /**
//  *  Pyramidal element with Bergot-Cohen-Durufle basis functions of order 1.
//  */
// using P1_Pyramid_BCD = Pk_Pyramid_BCD< 1 >;
// /**
//  *  Pyramidal element with Bergot-Cohen-Durufle basis functions of order 2.
//  */
// using P2_Pyramid_BCD = Pk_Pyramid_BCD< 2 >;
// /**
//  *  Pyramidal element with Bergot-Cohen-Durufle basis functions of order 3.
//  */
// using P3_Pyramid_BCD = Pk_Pyramid_BCD< 3 >;
// //using P4_Pyramid_BCD = Pk_Pyramid_BCD< 4 >;
// //using P5_Pyramid_BCD = Pk_Pyramid_BCD< 5 >;

}
}

#endif // GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_PKPYRAMIDBCD_HPP_
