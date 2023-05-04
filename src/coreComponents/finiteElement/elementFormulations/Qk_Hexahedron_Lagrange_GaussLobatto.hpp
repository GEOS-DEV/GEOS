/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file Qk_Hexahedron_Lagrange_GaussLobatto.hpp
 */

#ifndef GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_Q1HEXAHEDRON
#define GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_Q1HEXAHEDRON

#include "FiniteElementBase.hpp"
#include "LagrangeBasis1.hpp"
#include "LagrangeBasis2.hpp"
#include "LagrangeBasis3GL.hpp"
#include "LagrangeBasis4GL.hpp"
#include "LagrangeBasis5GL.hpp"
#include <utility>



namespace geosx
{
namespace finiteElement
{

/**
 * This class is the basis class for the hexahedron finite element cells with
 * shape functions defined on Gauss-Lobatto quadrature points.
 * All the degree-specific versions (Q1, Q2, Q3, ...) are defined at the end of this file.
 */
template< typename GL_BASIS >
class Qk_Hexahedron_Lagrange_GaussLobatto final : public FiniteElementBase
{
public:

  /// The number of nodes/support points per element per dimension.
  constexpr static localIndex num1dNodes = GL_BASIS::numSupportPoints;

  /// The number of nodes/support points per element.
  constexpr static localIndex numNodes = GL_BASIS::TensorProduct3D::numSupportPoints;

  /// The number of nodes/support points per face
  constexpr static localIndex numNodesPerFace = GL_BASIS::TensorProduct2D::numSupportPoints;

  /// The maximum number of support points per element.
  constexpr static localIndex maxSupportPoints = numNodes;

  /// The number of quadrature points per element.
  constexpr static localIndex numQuadraturePoints = numNodes;

  /** @cond Doxygen_Suppress */
  USING_FINITEELEMENTBASE
  /** @endcond Doxygen_Suppress */

  virtual ~Qk_Hexahedron_Lagrange_GaussLobatto() override
  {}

  GEOSX_HOST_DEVICE
  virtual localIndex getNumQuadraturePoints() const override
  {
    return numQuadraturePoints;
  }

  /**
   * @brief Get the number of quadrature points.
   * @param stack Stack variables as filled by @ref setupStack.
   * @return The number of quadrature points.
   */
  GEOSX_HOST_DEVICE
  static localIndex getNumQuadraturePoints( StackVariables const & stack )
  {
    GEOSX_UNUSED_VAR( stack );
    return numQuadraturePoints;
  }

  GEOSX_HOST_DEVICE
  virtual localIndex getNumSupportPoints() const override
  {
    return numNodes;
  }

  GEOSX_HOST_DEVICE
  virtual localIndex getMaxSupportPoints() const override
  {
    return maxSupportPoints;
  }

  /**
   * @brief Get the number of support points.
   * @param stack Object that holds stack variables.
   * @return The number of support points.
   */
  GEOSX_HOST_DEVICE
  static localIndex getNumSupportPoints( StackVariables const & stack )
  {
    GEOSX_UNUSED_VAR( stack );
    return numNodes;
  }


  /**
   * @brief Calculate shape functions values at a single point.
   * @param[in] coords The parent coordinates at which to evaluate the shape function value
   * @param[out] N The shape function values.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void calcN( real64 const (&coords)[3],
                     real64 (& N)[numNodes] )
  {
    GL_BASIS::TensorProduct3D::value( coords, N );
  }


  /**
   * @brief Calculate shape functions values for each support point at a
   *   quadrature point.
   * @param q Index of the quadrature point.
   * @param N An array to pass back the shape function values for each support
   *   point.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void calcN( localIndex const q,
                     real64 (& N)[numNodes] )
  {
    int qa, qb, qc;
    GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );
    real64 const qCoords[3] = { GL_BASIS::parentSupportCoord( qa ),
                                GL_BASIS::parentSupportCoord( qb ),
                                GL_BASIS::parentSupportCoord( qc ) };

    GL_BASIS::TensorProduct3D::value( qCoords, N );
  }

  /**
   * @brief Calculate shape functions values for each support point at a
   *   quadrature point.
   * @param q Index of the quadrature point.
   * @param stack Variables allocated on the stack as filled by @ref setupStack.
   * @param N An array to pass back the shape function values for each support
   *   point.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void calcN( localIndex const q,
                     StackVariables const & stack,
                     real64 ( & N )[numNodes] )
  {
    GEOSX_UNUSED_VAR( stack );
    return calcN( q, N );
  }

  /**
   * @brief Calculate the shape functions derivatives wrt the physical
   *   coordinates.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the support points.
   * @param gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   * @return The determinant of the parent/physical transformation matrix.
   */
  GEOSX_HOST_DEVICE
  static real64 calcGradN( localIndex const q,
                           real64 const (&X)[numNodes][3],
                           real64 ( &gradN )[numNodes][3] );

  /**
   * @brief Calculate the shape functions derivatives wrt the physical
   *   coordinates at a single point.
   * @param[in] coords The parent coordinates at which to evaluate the shape function value
   * @param[in] X Array containing the coordinates of the support points.
   * @param[out] gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   * @return The determinant of the parent/physical transformation matrix.
   */
  GEOSX_HOST_DEVICE
  static real64 calcGradN( real64 const (&coords)[3],
                           real64 const (&X)[numNodes][3],
                           real64 ( &gradN )[numNodes][3] );

  /**
   * @brief Calculate the shape functions derivatives wrt the physical
   *   coordinates.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the support points.
   * @param stack Variables allocated on the stack as filled by @ref setupStack.
   * @param gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   * @return The determinant of the parent/physical transformation matrix.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
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
  GEOSX_HOST_DEVICE
  static real64 transformedQuadratureWeight( localIndex const q,
                                             real64 const (&X)[numNodes][3] );

  /**
   * @brief Calculates the isoparametric "Jacobian" transformation
   *   matrix/mapping from the parent space to the physical space on a 2D domain (face).
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param X Array containing the coordinates of the support points.
   * @param J Array to store the Jacobian transformation.
   */
  GEOSX_HOST_DEVICE
  static void jacobianTransformation2d( int const qa,
                                        int const qb,
                                        real64 const (&X)[numNodesPerFace][3],
                                        real64 ( &J )[3][2] );


  /**
   * @brief Calculates the isoparametric "Jacobian" transformation
   *   matrix/mapping from the parent space to the physical space.
   * @param q The quadrature point index in 3d space.
   * @param X Array containing the coordinates of the support points.
   * @param J Array to store the Jacobian transformation.
   * @return The determinant of the Jacobian transformation matrix.
   */
  GEOSX_HOST_DEVICE
  static real64 invJacobianTransformation( int const q,
                                           real64 const (&X)[numNodes][3],
                                           real64 ( & J )[3][3] )
  {
    int qa, qb, qc;
    GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );
    jacobianTransformation( qa, qb, qc, X, J );
    return LvArray::tensorOps::invert< 3 >( J );
  }


  /**
   * @brief Calculate the symmetric gradient of a vector valued support field
   *   at a quadrature point using the stored inverse of the Jacobian
   *   transformation matrix.
   * @param q The linear index of the quadrature point.
   * @param invJ The inverse of the Jacobian transformation matrix.
   * @param var The vector valued support field to apply the gradient
   *   operator on.
   * @param grad The symmetric gradient in Voigt notation.
   */
  GEOSX_HOST_DEVICE
  static void symmetricGradient( int const q,
                                 real64 const (&invJ)[3][3],
                                 real64 const (&var)[numNodes][3],
                                 real64 ( &grad )[6] );



  /**
   * @brief Calculate the gradient of a vector valued support field at a point
   *   using the stored basis function gradients for all support points.
   * @param q The linear index of the quadrature point.
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
  GEOSX_HOST_DEVICE
  static void gradient( int const q,
                        real64 const (&invJ)[3][3],
                        real64 const (&var)[numNodes][3],
                        real64 ( &grad )[3][3] );


  /**
   * @brief Inner product of all basis function gradients and a rank-2
   *   symmetric tensor evaluated at a quadrature point.
   * @param q The linear index of the quadrature point.
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
  GEOSX_HOST_DEVICE
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
   * @param X Array containing the coordinates of the support points.
   * @param J Array to store the Jacobian transformation.
   */
  GEOSX_HOST_DEVICE
  static void jacobianTransformation( int const qa,
                                      int const qb,
                                      int const qc,
                                      real64 const (&X)[numNodes][3],
                                      real64 ( &J )[3][3] );

  /**
   * @brief Calculates the isoparametric "Jacobian" transformation
   *   matrix/mapping from the parent space to the physical space at a single point.
   * @param coords The parent coordinates at which to evaluate the shape function value
   * @param X Array containing the coordinates of the support points.
   * @param J Array to store the Jacobian transformation.
   */
  GEOSX_HOST_DEVICE
  static void jacobianTransformation( real64 const (&coords)[3],
                                      real64 const (&X)[numNodes][3],
                                      real64 ( &J )[3][3] );

  /**
   * @brief computes the non-zero contributions of the d.o.f. indexd by q to the
   *   mass matrix M, i.e., the superposition matrix of the shape functions.
   * @param q The quadrature point index
   * @param X Array containing the coordinates of the support points.
   * @return The diagonal mass term associated to q
   */
  GEOSX_HOST_DEVICE
  static real64
  computeMassTerm( int q,
                   real64 const (&X)[numNodes][3] );

  /**
   * @brief computes the non-zero contributions of the d.o.f. indexd by q to the
   *   damping matrix M, i.e., the superposition matrix of the shape functions
   *   integrated over a face.
   * @param q The quadrature point index
   * @param X Array containing the coordinates of the support points.
   * @return The diagonal damping term associated to q
   */
  GEOSX_HOST_DEVICE
  static real64
  computeDampingTerm( int q,
                      real64 const (&X)[numNodesPerFace][3] );

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
  GEOSX_HOST_DEVICE
  static void
    computeBMatrix( int const qa,
                    int const qb,
                    int const qc,
                    real64 const (&X)[numNodes][3],
                    real64 ( &J )[3][3],
                    real64 ( &B )[6] );

  /**
   * @brief computes the non-zero contributions of the d.o.f. indexd by q to the
   *   stiffness matrix R, i.e., the superposition matrix of first derivatives
   *   of the shape functions.
   * @param q The quadrature point index
   * @param X Array containing the coordinates of the support points.
   * @param func Callback function accepting three parameters: i, j and R_ij
   */
  template< typename FUNC >
  GEOSX_HOST_DEVICE
  static void
  computeStiffnessTerm( int q,
                        real64 const (&X)[numNodes][3],
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
  GEOSX_HOST_DEVICE
  static void
  computeFirstOrderStiffnessTermX( int q,
                                   real64 const (&X)[numNodes][3],
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
  GEOSX_HOST_DEVICE
  static void
  computeFirstOrderStiffnessTermY( int q,
                                   real64 const (&X)[numNodes][3],
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
  GEOSX_HOST_DEVICE
  static void
  computeFirstOrderStiffnessTermZ( int q,
                                   real64 const (&X)[numNodes][3],
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
  GEOSX_HOST_DEVICE
  static void
  computeFirstOrderStiffnessTerm( int q,
                                  real64 const (&X)[numNodes][3],
                                  FUNC && stiffnessVal );


  /**
   * @brief Apply a Jacobian transformation matrix from the parent space to the
   *   physical space on the parent shape function derivatives, producing the
   *   shape function derivatives in the physical space.
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param qc The 1d quadrature point index in xi2 direction (0,1)
   * @param invJ The Jacobian transformation from parent->physical space.
   * @param gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   */
  GEOSX_HOST_DEVICE
  static void
    applyTransformationToParentGradients( int const qa,
                                          int const qb,
                                          int const qc,
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
  GEOSX_HOST_DEVICE
  static void
    applyTransformationToParentGradients( real64 const (&coords)[3],
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
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param qc The 1d quadrature point index in xi2 direction (0,1)
   * @param func The function to call within the support loop.
   * @param params The parameters to pass to @p func.
   */
  template< typename FUNC, typename ... PARAMS >
  GEOSX_HOST_DEVICE
  static void supportLoop( int const qa,
                           int const qb,
                           int const qc,
                           FUNC && func,
                           PARAMS &&... params );
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
  GEOSX_HOST_DEVICE
  static void supportLoop( real64 const (&coords)[3],
                           FUNC && func,
                           PARAMS &&... params );

  /**
   * @brief Applies a function inside a generic loop in over the tensor product
   *   indices, over a 2d domain.
   * @tparam FUNC The type of function to call within the support loop.
   * @tparam PARAMS The parameter pack types to pass through to @p FUNC.
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param func The function to call within the support loop.
   * @param params The parameters to pass to @p func.
   */
  template< typename FUNC, typename ... PARAMS >
  GEOSX_HOST_DEVICE
  static void supportLoop2d( int const qa,
                             int const qb,
                             FUNC && func,
                             PARAMS &&... params );
};

/// @cond Doxygen_Suppress

template< typename GL_BASIS >
template< typename FUNC, typename ... PARAMS >
GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE void
Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::supportLoop( int const qa,
                                                              int const qb,
                                                              int const qc,
                                                              FUNC && func,
                                                              PARAMS &&... params )
{

  real64 const qCoords[3] = { GL_BASIS::parentSupportCoord( qa ),
                              GL_BASIS::parentSupportCoord( qb ),
                              GL_BASIS::parentSupportCoord( qc ) };
  supportLoop( qCoords, std::forward< FUNC >( func ), std::forward< PARAMS >( params )... );
}


template< typename GL_BASIS >
template< typename FUNC, typename ... PARAMS >
GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE void
Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::supportLoop( real64 const (&coords)[3],
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
GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE void
Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::supportLoop2d( int const qa,
                                                                int const qb,
                                                                FUNC && func,
                                                                PARAMS &&... params )
{

  real64 const qCoords[2] = { GL_BASIS::parentSupportCoord( qa ),
                              GL_BASIS::parentSupportCoord( qb ) };

  for( int b=0; b<num1dNodes; ++b )
  {
    for( int a=0; a<num1dNodes; ++a )
    {
      real64 const dNdXi[2] = { GL_BASIS::gradient( a, qCoords[0] )*
                                GL_BASIS::value( b, qCoords[1] ),
                                GL_BASIS::value( a, qCoords[0] )*
                                GL_BASIS::gradient( b, qCoords[1] ) };

      localIndex nodeIndex = GL_BASIS::TensorProduct2D::linearIndex( a, b );
      func( dNdXi, nodeIndex, std::forward< PARAMS >( params )... );
    }
  }
}

//*************************************************************************************************
template< typename GL_BASIS >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::calcGradN( localIndex const q,
                                                            real64 const (&X)[numNodes][3],
                                                            real64 (& gradN)[numNodes][3] )
{
  real64 J[3][3] = {{0}};

  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );

  jacobianTransformation( qa, qb, qc, X, J );

  real64 const detJ = LvArray::tensorOps::invert< 3 >( J );

  applyTransformationToParentGradients( qa, qb, qc, J, gradN );

  return detJ;
}

//*************************************************************************************************
template< typename GL_BASIS >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::calcGradN( real64 const (&coords)[3],
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
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::
calcGradN( localIndex const q,
           real64 const (&X)[numNodes][3],
           StackVariables const & GEOSX_UNUSED_PARAM( stack ),
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
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::
jacobianTransformation( int const qa,
                        int const qb,
                        int const qc,
                        real64 const (&X)[numNodes][3],
                        real64 ( & J )[3][3] )
{
  supportLoop( qa, qb, qc, [] GEOSX_HOST_DEVICE ( real64 const (&dNdXi)[3],
                                                  int const nodeIndex,
                                                  real64 const (&X)[numNodes][3],
                                                  real64 (& J)[3][3] )
  {
    real64 const * const GEOSX_RESTRICT Xnode = X[nodeIndex];
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        J[i][j] = J[i][j] + dNdXi[ j ] * Xnode[i];
      }
    }

//    J[0][0] = J[0][0] + dNdXi[0] * Xnode[0];
//    J[0][1] = J[0][1] + dNdXi[1] * Xnode[0];
//    J[0][2] = J[0][2] + dNdXi[2] * Xnode[0];
//    J[1][0] = J[1][0] + dNdXi[0] * Xnode[1];
//    J[1][1] = J[1][1] + dNdXi[1] * Xnode[1];
//    J[1][2] = J[1][2] + dNdXi[2] * Xnode[1];
//    J[2][0] = J[2][0] + dNdXi[0] * Xnode[2];
//    J[2][1] = J[2][1] + dNdXi[1] * Xnode[2];
//    J[2][2] = J[2][2] + dNdXi[2] * Xnode[2];

  }, X, J );
}
template< typename GL_BASIS >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::
jacobianTransformation( real64 const (&coords)[3],
                        real64 const (&X)[numNodes][3],
                        real64 ( & J )[3][3] )
{
  supportLoop( coords, [] GEOSX_HOST_DEVICE ( real64 const (&dNdXi)[3],
                                              int const nodeIndex,
                                              real64 const (&X)[numNodes][3],
                                              real64 (& J)[3][3] )
  {
    real64 const * const GEOSX_RESTRICT Xnode = X[nodeIndex];
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
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::
jacobianTransformation2d( int const qa,
                          int const qb,
                          real64 const (&X)[numNodesPerFace][3],
                          real64 ( & J )[3][2] )
{
  supportLoop2d( qa, qb, [] GEOSX_HOST_DEVICE ( real64 const (&dNdXi)[2],
                                                int const nodeIndex,
                                                real64 const (&X)[numNodesPerFace][3],
                                                real64 ( & J)[3][2] )
  {
    real64 const * const GEOSX_RESTRICT Xnode = X[nodeIndex];
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 2; ++j )
      {
        J[i][j] = J[i][j] + dNdXi[ j ] * Xnode[i];
      }
    }
  }, X, J );
}

template< typename GL_BASIS >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::
computeMassTerm( int q,
                 real64 const (&X)[numNodes][3] )
{
  real64 J[3][3] = {{0}};
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );
  jacobianTransformation( qa, qb, qc, X, J );
  return LvArray::math::abs( LvArray::tensorOps::determinant< 3 >( J ) )*GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc );
}

template< typename GL_BASIS >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::
computeDampingTerm( int q,
                    real64 const (&X)[numNodesPerFace][3] )
{
  real64 B[3];
  real64 J[3][2] = {{0}};
  int qa, qb;
  GL_BASIS::TensorProduct2D::multiIndex( q, qa, qb );
  jacobianTransformation2d( qa, qb, X, J );
  // compute J^T.J, using Voigt notation for B
  B[0] = J[0][0]*J[0][0]+J[1][0]*J[1][0]+J[2][0]*J[2][0];
  B[1] = J[0][1]*J[0][1]+J[1][1]*J[1][1]+J[2][1]*J[2][1];
  B[2] = J[0][0]*J[0][1]+J[1][0]*J[1][1]+J[2][0]*J[2][1];
  return sqrt( LvArray::math::abs( LvArray::tensorOps::symDeterminant< 2 >( B ) ) )*GL_BASIS::weight( qa )*GL_BASIS::weight( qb );
}

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
template< typename GL_BASIS >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::
computeBMatrix( int const qa,
                int const qb,
                int const qc,
                real64 const (&X)[numNodes][3],
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
template< typename FUNC >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::
computeStiffnessTerm( int q,
                      real64 const (&X)[numNodes][3],
                      FUNC && func )
{
  real64 B[6] = {0};
  real64 J[3][3] = {{0}};
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );
  computeBMatrix( qa, qb, qc, X, J, B );
  // diagonal terms
  for( int i=0; i<num1dNodes; i++ )
  {
    for( int j=0; j<num1dNodes; j++ )
    {
      func( GL_BASIS::TensorProduct3D::linearIndex( i, qb, qc ),
            GL_BASIS::TensorProduct3D::linearIndex( j, qb, qc ),
            GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc )*B[0]*
            GL_BASIS::gradient( i, GL_BASIS::parentSupportCoord( qa ) )*
            GL_BASIS::gradient( j, GL_BASIS::parentSupportCoord( qa ) ) );
    }
  }
  for( int i=0; i<num1dNodes; i++ )
  {
    for( int j=0; j<num1dNodes; j++ )
    {
      func( GL_BASIS::TensorProduct3D::linearIndex( qa, i, qc ),
            GL_BASIS::TensorProduct3D::linearIndex( qa, j, qc ),
            GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc )*B[1]*
            GL_BASIS::gradient( i, GL_BASIS::parentSupportCoord( qb ) )*
            GL_BASIS::gradient( j, GL_BASIS::parentSupportCoord( qb ) ) );
    }
  }
  for( int i=0; i<num1dNodes; i++ )
  {
    for( int j=0; j<num1dNodes; j++ )
    {
      func( GL_BASIS::TensorProduct3D::linearIndex( qa, qb, i ),
            GL_BASIS::TensorProduct3D::linearIndex( qa, qb, j ),
            GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc )*B[2]*
            GL_BASIS::gradient( i, GL_BASIS::parentSupportCoord( qc ) )*
            GL_BASIS::gradient( j, GL_BASIS::parentSupportCoord( qc ) ) );
    }
  }
  // off-diagonal terms
  for( int i=0; i<num1dNodes; i++ )
  {
    for( int j=0; j<num1dNodes; j++ )
    {
      int ii = GL_BASIS::TensorProduct3D::linearIndex( qa, i, qc );
      int jj = GL_BASIS::TensorProduct3D::linearIndex( qa, qb, j );
      real64 val = GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc )*B[3]*
                   GL_BASIS::gradient( i, GL_BASIS::parentSupportCoord( qb ) )*
                   GL_BASIS::gradient( j, GL_BASIS::parentSupportCoord( qc ) );
      func( ii, jj, val );
      func( jj, ii, val );
    }
  }
  for( int i=0; i<num1dNodes; i++ )
  {
    for( int j=0; j<num1dNodes; j++ )
    {
      int ii = GL_BASIS::TensorProduct3D::linearIndex( i, qb, qc );
      int jj = GL_BASIS::TensorProduct3D::linearIndex( qa, qb, j );
      real64 val = GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc )*B[4]*
                   GL_BASIS::gradient( i, GL_BASIS::parentSupportCoord( qa ) )*
                   GL_BASIS::gradient( j, GL_BASIS::parentSupportCoord( qc ) );
      func( ii, jj, val );
      func( jj, ii, val );
    }
  }
  for( int i=0; i<num1dNodes; i++ )
  {
    for( int j=0; j<num1dNodes; j++ )
    {
      int ii = GL_BASIS::TensorProduct3D::linearIndex( i, qb, qc );
      int jj = GL_BASIS::TensorProduct3D::linearIndex( qa, j, qc );
      real64 val = GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc )*B[5]*
                   GL_BASIS::gradient( i, GL_BASIS::parentSupportCoord( qa ) )*
                   GL_BASIS::gradient( j, GL_BASIS::parentSupportCoord( qb ) );
      func( ii, jj, val );
      func( jj, ii, val );
    }
  }
}

template< typename GL_BASIS >
template< typename FUNC >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::
computeFirstOrderStiffnessTerm( int q,
                                real64 const (&X)[numNodes][3],
                                FUNC && func )
{
  real64 J[3][3] = {{0}};
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );
  jacobianTransformation( qa, qb, qc, X, J );
  real64 const detJ = LvArray::tensorOps::invert< 3 >( J );
  // diagonal terms
  for( int i=0; i<num1dNodes; i++ )
  {
    for( int j=0; j<num1dNodes; j++ )
    {
      func( GL_BASIS::TensorProduct3D::linearIndex( i, qb, qc ),
            GL_BASIS::TensorProduct3D::linearIndex( j, qb, qc ),
            detJ*GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc )*
            GL_BASIS::gradient( i, GL_BASIS::parentSupportCoord( qa ) )*
            GL_BASIS::gradient( j, GL_BASIS::parentSupportCoord( qa ) ),
            J,
            0,
            0 );
    }
  }
  for( int i=0; i<num1dNodes; i++ )
  {
    for( int j=0; j<num1dNodes; j++ )
    {
      func( GL_BASIS::TensorProduct3D::linearIndex( qa, i, qc ),
            GL_BASIS::TensorProduct3D::linearIndex( qa, j, qc ),
            detJ*GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc )*
            GL_BASIS::gradient( i, GL_BASIS::parentSupportCoord( qb ) )*
            GL_BASIS::gradient( j, GL_BASIS::parentSupportCoord( qb ) ),
            J,
            1,
            1 );
    }
  }
  for( int i=0; i<num1dNodes; i++ )
  {
    for( int j=0; j<num1dNodes; j++ )
    {
      func( GL_BASIS::TensorProduct3D::linearIndex( qa, qb, i ),
            GL_BASIS::TensorProduct3D::linearIndex( qa, qb, j ),
            detJ*GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc )*
            GL_BASIS::gradient( i, GL_BASIS::parentSupportCoord( qc ) )*
            GL_BASIS::gradient( j, GL_BASIS::parentSupportCoord( qc ) ),
            J,
            2,
            2 );
    }
  }
  // off-diagonal terms
  for( int i=0; i<num1dNodes; i++ )
  {
    for( int j=0; j<num1dNodes; j++ )
    {
      func( GL_BASIS::TensorProduct3D::linearIndex( qa, i, qc ),
            GL_BASIS::TensorProduct3D::linearIndex( qa, qb, j ),
            detJ*GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc )*
            GL_BASIS::gradient( i, GL_BASIS::parentSupportCoord( qb ) )*
            GL_BASIS::gradient( j, GL_BASIS::parentSupportCoord( qc ) ),
            J,
            1,
            2 );
    }
  }

  for( int i=0; i<num1dNodes; i++ )
  {
    for( int j=0; j<num1dNodes; j++ )
    {
      func( GL_BASIS::TensorProduct3D::linearIndex( qa, qb, j ),
            GL_BASIS::TensorProduct3D::linearIndex( qa, i, qc ),
            detJ*GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc )*
            GL_BASIS::gradient( i, GL_BASIS::parentSupportCoord( qb ) )*
            GL_BASIS::gradient( j, GL_BASIS::parentSupportCoord( qc ) ),
            J,
            2,
            1 );
    }
  }



  for( int i=0; i<num1dNodes; i++ )
  {
    for( int j=0; j<num1dNodes; j++ )
    {
      func( GL_BASIS::TensorProduct3D::linearIndex( i, qb, qc ),
            GL_BASIS::TensorProduct3D::linearIndex( qa, qb, j ),
            detJ*GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc )*
            GL_BASIS::gradient( i, GL_BASIS::parentSupportCoord( qa ) )*
            GL_BASIS::gradient( j, GL_BASIS::parentSupportCoord( qc ) ),
            J,
            0,
            2 );

    }
  }

  for( int i=0; i<num1dNodes; i++ )
  {
    for( int j=0; j<num1dNodes; j++ )
    {
      func( GL_BASIS::TensorProduct3D::linearIndex( qa, qb, j ),
            GL_BASIS::TensorProduct3D::linearIndex( i, qb, qc ),
            detJ*GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc )*
            GL_BASIS::gradient( i, GL_BASIS::parentSupportCoord( qa ) )*
            GL_BASIS::gradient( j, GL_BASIS::parentSupportCoord( qc ) ),
            J,
            2,
            0 );

    }
  }

  for( int i=0; i<num1dNodes; i++ )
  {
    for( int j=0; j<num1dNodes; j++ )
    {
      func( GL_BASIS::TensorProduct3D::linearIndex( i, qb, qc ),
            GL_BASIS::TensorProduct3D::linearIndex( qa, j, qc ),
            detJ*GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc )*
            GL_BASIS::gradient( i, GL_BASIS::parentSupportCoord( qa ) )*
            GL_BASIS::gradient( j, GL_BASIS::parentSupportCoord( qb ) ),
            J,
            0,
            1 );
    }
  }

  for( int i=0; i<num1dNodes; i++ )
  {
    for( int j=0; j<num1dNodes; j++ )
    {
      func( GL_BASIS::TensorProduct3D::linearIndex( qa, j, qc ),
            GL_BASIS::TensorProduct3D::linearIndex( i, qb, qc ),
            detJ*GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc )*
            GL_BASIS::gradient( i, GL_BASIS::parentSupportCoord( qa ) )*
            GL_BASIS::gradient( j, GL_BASIS::parentSupportCoord( qb ) ),
            J,
            1,
            0 );
    }
  }

}

template< typename GL_BASIS >
template< typename FUNC >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::
computeFirstOrderStiffnessTermX( int q,
                                 real64 const (&X)[numNodes][3],
                                 FUNC && func )
{
  real64 J[3][3] = {{0}};
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );
  jacobianTransformation( qa, qb, qc, X, J );
  real64 const detJ = LvArray::tensorOps::invert< 3 >( J );

  for( int i1 = 0; i1 < num1dNodes; ++i1 )
  {
    real64 val = GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc )*GL_BASIS::gradient( i1, GL_BASIS::parentSupportCoord( qa ) );
    func( GL_BASIS::TensorProduct3D::linearIndex( i1, qb, qc ),
          GL_BASIS::TensorProduct3D::linearIndex( qa, qb, qc ),
          detJ*J[0][0]*val, detJ*J[0][1]*val, detJ*J[0][2]*val );

  }

}

template< typename GL_BASIS >
template< typename FUNC >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::
computeFirstOrderStiffnessTermY( int q,
                                 real64 const (&X)[numNodes][3],
                                 FUNC && func )
{
  real64 J[3][3] = {{0}};
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );
  jacobianTransformation( qa, qb, qc, X, J );
  real64 const detJ = LvArray::tensorOps::invert< 3 >( J );

  for( int i2 = 0; i2 < num1dNodes; ++i2 )
  {
    real64 val = GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc )*GL_BASIS::gradient( i2, GL_BASIS::parentSupportCoord( qb ) );
    func( GL_BASIS::TensorProduct3D::linearIndex( qa, i2, qc ),
          GL_BASIS::TensorProduct3D::linearIndex( qa, qb, qc ),
          detJ*J[1][0]*val, detJ*J[1][1]*val, detJ*J[1][2]*val );

  }

}

template< typename GL_BASIS >
template< typename FUNC >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::
computeFirstOrderStiffnessTermZ( int q,
                                 real64 const (&X)[numNodes][3],
                                 FUNC && func )
{
  real64 J[3][3] = {{0}};
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );
  jacobianTransformation( qa, qb, qc, X, J );
  real64 const detJ = LvArray::tensorOps::invert< 3 >( J );

  for( int i3 = 0; i3 < num1dNodes; ++i3 )
  {
    real64 val = GL_BASIS::weight( qa )*GL_BASIS::weight( qb )*GL_BASIS::weight( qc )*GL_BASIS::gradient( i3, GL_BASIS::parentSupportCoord( qc ) );
    func( GL_BASIS::TensorProduct3D::linearIndex( qa, qb, i3 ),
          GL_BASIS::TensorProduct3D::linearIndex( qa, qb, qc ),
          detJ*J[2][0]*val, detJ*J[2][1]*val, detJ*J[2][2]*val );

  }

}

//*************************************************************************************************
template< typename GL_BASIS >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::
applyTransformationToParentGradients( int const qa,
                                      int const qb,
                                      int const qc,
                                      real64 const ( &invJ )[3][3],
                                      real64 (& gradN)[numNodes][3] )
{
  supportLoop( qa, qb, qc, [] GEOSX_HOST_DEVICE ( real64 const (&dNdXi)[3],
                                                  int const nodeIndex,
                                                  real64 const (&invJ)[3][3],
                                                  real64 (& gradN)[numNodes][3] )
  {
//    for( int i = 0; i < 3; ++i )
//    {
//      gradN[nodeIndex][i] = 0.0;
//      for( int j = 0; j < 3; ++j )
//      {
//        gradN[nodeIndex][i] = gradN[nodeIndex][i] + dNdXi[ j ] * invJ[j][i];
//      }
//    }
    // smaller register footprint by manually unrolling the for loops.
    gradN[nodeIndex][0] = dNdXi[0] * invJ[0][0] + dNdXi[1] * invJ[1][0] + dNdXi[2] * invJ[2][0];
    gradN[nodeIndex][1] = dNdXi[0] * invJ[0][1] + dNdXi[1] * invJ[1][1] + dNdXi[2] * invJ[2][1];
    gradN[nodeIndex][2] = dNdXi[0] * invJ[0][2] + dNdXi[1] * invJ[1][2] + dNdXi[2] * invJ[2][2];


  }, invJ, gradN );
}

//*************************************************************************************************
template< typename GL_BASIS >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::
applyTransformationToParentGradients( real64 const (&coords)[3],
                                      real64 const ( &invJ )[3][3],
                                      real64 (& gradN)[numNodes][3] )
{
  supportLoop( coords, [] GEOSX_HOST_DEVICE ( real64 const (&dNdXi)[3],
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
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::
transformedQuadratureWeight( localIndex const q,
                             real64 const (&X)[numNodes][3] )
{
  real64 J[3][3] = {{0}};

  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );

  jacobianTransformation( qa, qb, qc, X, J );

  return LvArray::tensorOps::determinant< 3 >( J );
}



template< typename GL_BASIS >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::
symmetricGradient( int const q,
                   real64 const (&invJ)[3][3],
                   real64 const (&var)[numNodes][3],
                   real64 (& grad)[6] )
{
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );

  supportLoop( qa, qb, qc, [] GEOSX_HOST_DEVICE ( real64 const (&dNdXi)[3],
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
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::
plusGradNajAij( int const q,
                real64 const (&invJ)[3][3],
                real64 const (&var)[6],
                real64 (& R)[numNodes][3] )
{
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );

  supportLoop( qa, qb, qc,
               [] GEOSX_HOST_DEVICE
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
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void Qk_Hexahedron_Lagrange_GaussLobatto< GL_BASIS >::
gradient( int const q,
          real64 const (&invJ)[3][3],
          real64 const (&var)[numNodes][3],
          real64 (& grad)[3][3] )
{
  int qa, qb, qc;
  GL_BASIS::TensorProduct3D::multiIndex( q, qa, qb, qc );

  supportLoop( qa, qb, qc, [] GEOSX_HOST_DEVICE ( real64 const (&dNdXi)[3],
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
 * This class contains the kernel accessible functions specific to the standard
 * Trilinear Hexahedron finite element with a Gaussian quadrature rule. It is
 * assumed that the indexing for the quadrature points mirrors that of the
 * nodes. Also note that the assumed node ordering is not the standard
 * right-hand-rule used in the literature. Here we use a Cartesian aligned
 * numbering in order to simplify the mapping to the parent coordinates and
 * tensor product indices.
 *
 *                  6                   7                       ____________________
 *                   o-----------------o                       |Node   xi0  xi1  xi2|
 *                  /.                /|                       |=====  ===  ===  ===|
 *                 / .               / |                       | 0     -1   -1   -1 |
 *              4 o-----------------o 5|                       | 1      1   -1   -1 |
 *                |  .              |  |                       | 2     -1    1   -1 |
 *                |  .              |  |                       | 3      1    1   -1 |
 *                |  .              |  |                       | 4     -1   -1    1 |
 *                |  .              |  |                       | 5      1   -1    1 |
 *                |2 o..............|..o 3       xi2           | 6     -1    1    1 |
 *                | ,               | /          |             | 7      1    1    1 |
 *                |,                |/           | / xi1       |____________________|
 *                o-----------------o            |/
 *               0                   1           ------ xi0
 *
 *
 */
using Q1_Hexahedron_Lagrange_GaussLobatto = Qk_Hexahedron_Lagrange_GaussLobatto< LagrangeBasis1 >;
/**
 * This class contains the kernel accessible functions specific to the standard
 * Trilinear Hexahedron finite element with a Gaussian quadrature rule. It is
 * assumed that the indexing for the quadrature points mirrors that of the
 * nodes. Also note that the assumed node ordering is not the standard
 * right-hand-rule used in the literature. Here we use a Cartesian aligned
 * numbering in order to simplify the mapping to the parent coordinates and
 * tensor product indices.
 *
 *                                                                  ____________________
 *                                                                 |Node   xi0  xi1  xi2|
 *                                                                 |=====  ===  ===  ===|
 *                                                                 |  0    -1   -1   -1 |
 *                                                                 |  1     0   -1   -1 |
 *                                                                 |  2     1   -1   -1 |
 *              24              25               26                |  3    -1    0   -1 |
 *                o--------------o--------------o                  |  4     0    0   -1 |
 *               /.                            /|                  |  5     1    0   -1 |
 *              / .                           / |                  |  6    -1    1   -1 |
 *          21 o  .           o 22        23 o  |                  |  7     0    1   -1 |
 *            /   .                         /   |                  |  8     1    1   -1 |
 *           /    .         19             /    |                  |  9    -1   -1    0 |
 *       18 o--------------o--------------o 20  |                  | 10     0   -1    0 |
 *          |     o              o        |     o                  | 11     1   -1    0 |
 *          |     .15             16      |     |17                | 12    -1    0    0 |
 *          |     .                       |     |                  | 13     0    0    0 |
 *          |  o  .           o           |  o  |                  | 14     1    0    0 |
 *          |   12.            13         |   14|                  | 15    -1    1    0 |
 *          |     .                       |     |                  | 16     0    1    0 |
 *        9 o     .        o 10           o 11  |                  | 17     1    1    0 |
 *          |     o..............o........|.....o                  | 18    -1   -1    1 |
 *          |    , 6              7       |    / 8                 | 19     0   -1    1 |
 *          |   ,                         |   /                    | 20     1   -1    1 |
 *          |  o              o           |  o         xi2         | 21    -1    0    1 |
 *          | , 3              4          | / 5        |           | 22     0    0    1 |
 *          |,                            |/           | / xi1     | 23     1    0    1 |
 *          o--------------o--------------o            |/          | 24    -1    1    1 |
 *         0                1              2           o----- xi0  | 25     0    1    1 |
 *                                                                 | 26     1    1    1 |
 *                                                                 |____________________|
 *
 */
using Q2_Hexahedron_Lagrange_GaussLobatto = Qk_Hexahedron_Lagrange_GaussLobatto< LagrangeBasis2 >;
/**
 * This class contains the kernel accessible functions specific to the standard
 * Trilinear Hexahedron finite element with a Gaussian quadrature rule. It is
 * assumed that the indexing for the quadrature points mirrors that of the
 * nodes. Also note that the assumed node ordering is not the standard
 * right-hand-rule used in the literature. Here we use a Cartesian aligned
 * numbering in order to simplify the mapping to the parent coordinates and
 * tensor product indices.
 *
 *
 *                                                                  _____________________________________
 *                                                                 |Node      xi0         xi1         xi2|
 *                                                                 |=====     ===         ===         ===|
 *                                                                 |  0       -1          -1          -1 |
 *                                                                 |  1   -1/sqrt(5)      -1          -1 |
 *                                                                 |  2    1/sqrt(5)      -1          -1 |
 *              60       61         62        63                   |  3        1          -1          -1 |
 *                o---------o---------o---------o                  |  4       -1      -1/sqrt(5)      -1 |
 *            56 /.     57        58        59 /|                  |  5   -1/sqrt(5)  -1/sqrt(5)      -1 |
 *              o .       o         o         o |                  |  6    1/sqrt(5)  -1/sqrt(5)      -1 |
 *          52 /  .   53        54        55 /  |                  |  7        1      -1/sqrt(5)      -1 |
 *            o   .     o         o         o   |                  |  8       -1       1/sqrt(5)      -1 |
 *        48 /    o 49      o 50      o 51 /    o                  |  9   -1/sqrt(5)   1/sqrt(5)      -1 |
 *          o---------o---------o---------o     |                  | 10    1/sqrt(5)   1/sqrt(5)      -1 |
 *          |   o .       o         o     |   o |                  | 11        1       1/sqrt(5)      -1 |
 *          |     .                       |     |                  | 12       -1           1          -1 |
 *          | o   o     o   o     o   o   | o   o                  | 13   -1/sqrt(5)       1          -1 |
 *          |     .                       |     |                  | 14    1/sqrt(5)       1          -1 |
 *          o   o .   o   o     o   o     o   o |                  | 15        1           1          -1 |
 *          |     .                       |     |                  | ..       ..          ..          .. |
 *          | o   .     o         o       | o   |                  | ..       ..          ..          .. |
 *          |     o.........o.........o...|.....o                  | 55        1      -1/sqrt(5)       1 |
 *          o    ,12  o     13  o     14  o    /15                 | 56       -1       1/sqrt(5)       1 |
 *          |   o         o         o     |   o                    | 57   -1/sqrt(5)   1/sqrt(5)       1 |
 *          |  ,8         9         10    |  /11       xi2         | 58    1/sqrt(5)   1/sqrt(5)       1 |
 *          | o         o         o       | o          |           | 59        1       1/sqrt(5)       1 |
 *          |,4         5         6       |/7          | / xi1     | 60       -1           1           1 |
 *          o---------o---------o---------o            |/          | 61   -1/sqrt(5)       1           1 |
 *         0         1         2         3             o----- xi0  | 62    1/sqrt(5)       1           1 |
 *                                                                 | 63        1           1           1 |
 *                                                                 |_____________________________________|
 *
 */
using Q3_Hexahedron_Lagrange_GaussLobatto = Qk_Hexahedron_Lagrange_GaussLobatto< LagrangeBasis3GL >;
/**
 * This class contains the kernel accessible functions specific to the standard
 * Trilinear Hexahedron finite element with a Gaussian quadrature rule. It is
 * assumed that the indexing for the quadrature points mirrors that of the
 * nodes. Also note that the assumed node ordering is not the standard
 * right-hand-rule used in the literature. Here we use a Cartesian aligned
 * numbering in order to simplify the mapping to the parent coordinates and
 * tensor product indices.
 *                                                                  _____________________________________
 *                120      121     122     123       124           |Node      xi0         xi1         xi2 |
 *                  o-------o-------o-------o-------o              |=====     ===         ===         === |
 *                 /.                              /|              |   0       -1          -1          -1 |
 *            115 o .  116o    117o    118o    119o |              |   1   -sqrt(3/7)      -1          -1 |
 *               /  o                            /  o              |   2        0          -1          -1 |
 *          110 o   .111o    112o    113o    114o   |              |   3    sqrt(3/7)      -1          -1 |
 *             /  o .                          /  o |              |   4        1          -1          -1 |
 *        105 o     . o106    o107    o108 109o     |              |   5       -1      -sqrt(3/7)      -1 |
 *           /  o   o      102     103    104/  o   o              |   6   -sqrt(3/7)  -sqrt(3/7)      -1 |
 *      100 o-------o-------o-------o-------o       |              |   7        0      -sqrt(3/7)      -1 |
 *          | o   o . 101                   | o   o |              |   8    sqrt(3/7)  -sqrt(3/7)      -1 |
 *          |       .                       |       |              |   9        1      -sqrt(3/7)      -1 |
 *          o   o   o       o       o       o   o   o              |  10       -1           0          -1 |
 *          |       .                       |       |              |  11   -sqrt(3/7)       0          -1 |
 *          | o   o .20     21      22    23| o   o |24            |  12        0           0          -1 |
 *          |       o.......o.......o.......|.......o              |  13    sqrt(3/7)       0          -1 |
 *          o   o  ,o       o       o       o   o  /               |  14        1           0          -1 |
 *          |     o       o       o       o |     o                |  ..       ..          ..          .. |
 *          | o  ,15      13      17      18| o  /19               |  ..       ..          ..          .. |
 *          |   o       o       o       o   |   o                  |  ..       ..          ..          .. |
 *          o  ,10  o   11  o   12  o   13  o  /14     xi2         | 121        -1          1           1 |
 *          | o       o       o       o     | o        |           | 122    -sqrt(3/7)      1           1 |
 *          |,5       6       7       8     |/9        | / xi1     | 123         0          1           1 |
 *          o-------o-------o-------o-------o          |/          | 124     sqrt(3/7)      1           1 |
 *         0        1       2       3        4         o----- xi0  | 125         1          1           1 |
 *                                                                 |______________________________________|
 *
 */
using Q4_Hexahedron_Lagrange_GaussLobatto = Qk_Hexahedron_Lagrange_GaussLobatto< LagrangeBasis4GL >;
/**
 * This class contains the kernel accessible functions specific to the standard
 * Trilinear Hexahedron finite element with a Gauss Lobatto quadrature rule. It is
 * assumed that the indexing for the quadrature points mirrors that of the
 * nodes. Also note that the assumed node ordering is not the standard
 * right-hand-rule used in the literature. Here we use a Cartesian aligned
 * numbering in order to simplify the mapping to the parent coordinates and
 * tensor product indices.
 *
 *
 *                210      211      212      213      214      215    _______________________________________________________
 *                   o--------o--------o--------o--------o--------o  |Node      xi0                        xi1            xi2|
 *                  /.                                           /|  |=====     ===                        ===            ===|
 *            204  / .  205      206      207      208      209 / |  |  0       -1                         -1             -1 |
 *                o  .     o        o        o        o        o  |  |  1   -sqrt(1/21(7+/sqrt(7))         -1             -1 |
 *               /   o                                        /   |  |  2    -sqrt(1/21(7-/sqrt(7))        -1             -1 |
 *         198  /    .199     200      201      202      203 /    o  |  3    sqrt(1/21(7-/sqrt(7))         -1             -1 |
 *             o     .  o        o        o        o        o     |  |  4    sqrt(1/21(7+/sqrt(7))         -1             -1 |
 *            /      .                                     /      |  |  5        1                         -1             -1 |
 *      192  /   193 o     194      195      196      197 /    o  |  |  6       -1                 -sqrt(1/21(7+/sqrt(7)) -1 |
 *          o        o        o        o        o        o        o  |  7   -sqrt(1/21(7+/sqrt(7)) -sqrt(1/21(7+/sqrt(7)) -1 |
 *         /         .                                  /         |  |  8   -sqrt(1/21(7-/sqrt(7)) -sqrt(1/21(7+/sqrt(7)) -1 |
 *    186 /    187   .  188      189      190      191 /    o     |  |  9    sqrt(1/21(7-/sqrt(7)) -sqrt(1/21(7+/sqrt(7)) -1 |
 *       o        o  o     o        o        o        o        o  |  | 10    sqrt(1/21(7+/sqrt(7)) -sqrt(1/21(7+/sqrt(7)) -1 |
 *      /            .                               /            o  | 11        1                 -sqrt(1/21(7+/sqrt(7)) -1 |
 * 180 /    181      . 182    183      184      185 /    o        |  | ..       ..                         ..             .. |
 *    o--------o--------o--------o--------o--------o        o     |  | ..       ..                         ..             .. |
 *    |           o  .                             |           o  |  | 204      -1                  sqrt(1/21(7+/sqrt(7)) 1  |
 *    |  o           o        o        o        o  |  o  o        o  | 205  -sqrt(1/21(7+/sqrt(7))  sqrt(1/21(7+/sqrt(7)) 1  |
 *    |     o        .                             |     o        |  | 206  -sqrt(1/21(7-/sqrt(7))  sqrt(1/21(7-/sqrt(7)) 1  |
 *    |        o     .                             |        o     |  | 207  sqrt(1/21(7+/sqrt(7))   sqrt(1/21(7-/sqrt(7)) 1  |
 *    o           o  .                             o           o  |  | 208  sqrt(1/21(7-/sqrt(7))   sqrt(1/21(7+/sqrt(7)) 1  |
 *    |  o           .                             |  o           |  | 209       1                  sqrt(1/21(7+/sqrt( *  1  |
 *    |     o        o--------o--------o--------o--|-----o--------o  | 210      -1                          1             1  |
 *    |        o    ,30       31      32        33 |     34 o    /35 | 211  -sqrt(1/21(7+/sqrt(7))          1             1  |
 *    o            ,                               o            /    | 212  -sqrt(1/21(7-/sqrt(7))          1             1  |
 *    |  o        o        o         o       o     |  o        o     | 213   sqrt(1/21(7-/sqrt(7))          1             1  |
 *    |     o    ,24       25        26      27    |  28 o    /29    | 214   sqrt(1/21(7+/sqrt(7))          1             1  |
 *    |         ,                                  |         /       | 215       1                          1             1  |
 *    o        o        o         o       o     22 o        o *      |_______________________________________________________|
 *    |  o    ,18       19        20      21       |  o    /23
 *    |      ,                                     |      /
 *    |     o        o         o       o        o  |     o
 *    o    ,12       13        14      15       16 o    /17
 *    |   ,                                        |   /
 *    |  o        o        o        o        o     |  o               xi2
 *    | ,6        7        8        9        10    | /11               |
 *    |,                                           |/                  | / xi1
 *    o--------o--------o--------o--------o--------o                   |/
 *    0        1        2        3        4        5                   o----- xi0
 */
using Q5_Hexahedron_Lagrange_GaussLobatto = Qk_Hexahedron_Lagrange_GaussLobatto< LagrangeBasis5GL >;

/// @endcond

#if __GNUC__
#pragma GCC diagnostic pop
#endif
#undef PARENT_GRADIENT_METHOD
}
}

#endif //FINITE_ELEMENT_SHAPE_KERNEL
