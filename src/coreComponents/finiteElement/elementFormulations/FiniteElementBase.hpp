/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FiniteElementBase.hpp
 */

#if defined(GEOSX_USE_CUDA)
#define CALC_FEM_SHAPE_IN_KERNEL
#endif


#ifndef GEOSX_FINITEELEMENT_FINITEELEMENTBASE_HPP_
#define GEOSX_FINITEELEMENT_FINITEELEMENTBASE_HPP_

#include "common/DataTypes.hpp"
#include "common/GeosxMacros.hpp"


namespace geosx
{
namespace finiteElement
{

/**
 * @class FiniteElementShapeFunctionKernelBase
 * @brief Base class for the finite element kernels.
 */
class FiniteElementBase
{
public:

  /**
   * @brief Virtual getter for the number of quadrature points per element.
   * @return The number of quadrature points per element.
   */
  virtual localIndex getNumQuadraturePoints() const = 0;

  /**
   * @brief Virtual getter for the number of support points per element.
   * @return The number of support points per element.
   */
  virtual localIndex getNumSupportPoints() const = 0;

  /**
   * @brief Destructor
   */
  virtual ~FiniteElementBase() = default;

  /**
   * @brief Computes the inverse of a 3x3 c-array.
   * @param J The array to invert...which is also used to store the inverse.
   * @return The determinant of @p J.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static real64 inverse( real64 (& J)[3][3] )
  {
    real64 const temp[3][3] =
    { { J[1][1]*J[2][2] - J[1][2]*J[2][1], J[0][2]*J[2][1] - J[0][1]*J[2][2], J[0][1]*J[1][2] - J[0][2]*J[1][1] },
      { J[1][2]*J[2][0] - J[1][0]*J[2][2], J[0][0]*J[2][2] - J[0][2]*J[2][0], J[0][2]*J[1][0] - J[0][0]*J[1][2] },
      { J[1][0]*J[2][1] - J[1][1]*J[2][0], J[0][1]*J[2][0] - J[0][0]*J[2][1], J[0][0]*J[1][1] - J[0][1]*J[1][0] } };

    real64 const det =  J[0][0] * temp[0][0] + J[1][0] * temp[0][1] + J[2][0] * temp[0][2];
    real64 const invDet = 1.0 / det;

    for( int i=0; i<3; ++i )
    {
      for( int j=0; j<3; ++j )
      {
        J[i][j] = temp[i][j] * invDet;
      }
    }
    return det;
  }

  /**
   * @brief Calculate the determinant of a 3x3 c-array.
   * @param J The input array.
   * @return The determinant of @p J
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static real64 detJ( real64 const (&J)[3][3] )
  {
    return J[0][0] * ( J[1][1]*J[2][2] - J[1][2]*J[2][1] ) +
           J[1][0] * ( J[0][2]*J[2][1] - J[0][1]*J[2][2] ) +
           J[2][0] * ( J[0][1]*J[1][2] - J[0][2]*J[1][1] );
  }

  template< typename LEAF >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 getGradN( localIndex const k,
                   localIndex const q,
                   real64 const (&X)[LEAF::numNodes][3],
                   real64 ( & gradN )[LEAF::numNodes][3] ) const
  {
    GEOSX_UNUSED_VAR( k );
    return LEAF::calcGradN( q, X, gradN );
  }


  template< typename LEAF >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 getGradN( localIndex const k,
                   localIndex const q,
                   int const,
                   real64 ( & gradN )[LEAF::numNodes][3] ) const
  {
    for( int a=0; a<LEAF::numNodes; ++a )
    {
      gradN[a][0] = m_viewGradN( k, q, a, 0 );
      gradN[a][1] = m_viewGradN( k, q, a, 1 );
      gradN[a][2] = m_viewGradN( k, q, a, 2 );
    }
    return m_viewDetJ( k, q );
  }


  /**
   * @brief Calculate the symmetric gradient of a vector valued support field
   *   at a point using the stored basis function gradients for all support
   *   points.
   * @param gradN The basis function gradients at a point in the element.
   * @param var The vector valued support field that the gradient operator will
   *  be applied to.
   * @param grad The symmetric gradient in Voigt notation.
   *
   * More precisely, the operator is defined as:
   * \f[
   * grad^s_{ij}  = \frac{1}{2} \sum_a^{nSupport} \left ( \frac{\partial N_a}{\partial X_j} var_{ai} + \frac{\partial N_a}{\partial X_i}
   * var_{aj}\right ),
   * \f]
   *
   */
  template< int NUM_SUPPORT_POINTS,
            typename GRADIENT_TYPE >
  GEOSX_HOST_DEVICE
  static void symmetricGradient( GRADIENT_TYPE const & gradN,
                                 real64 const (&var)[NUM_SUPPORT_POINTS][3],
                                 real64 ( &gradVar )[6] );



  /**
   * @brief Calculate the gradient of a vector valued support field at a point
   *   using the stored basis function gradients for all support points.
   * @param gradN The basis function gradients at a point in the element.
   * @param var The vector valued support field that the gradient operator will
   *  be applied to.
   * @param grad The  gradient.
   *
   * More precisely, the operator is defined as:
   * \f[
   * grad_{ij}  = \sum_a^{nSupport} \left ( \frac{\partial N_a}{\partial X_j} var_{ai}\right ),
   * \f]
   */
  template< int NUM_SUPPORT_POINTS,
            typename GRADIENT_TYPE >
  GEOSX_HOST_DEVICE
  static void gradient( GRADIENT_TYPE const & gradN,
                        real64 const (&var)[NUM_SUPPORT_POINTS][3],
                        real64 ( &gradVar )[3][3] );

  /**
   * @brief Inner product of all basis function gradients and a rank-2
   *   symmetric tensor.
   * @param gradN The basis function gradients at a point in the element.
   * @param var The rank-2 symmetric tensor at @p q.
   * @param R The vector resulting from the tensor contraction.
   *
   * More precisely, the operator is defined as:
   * \f[
   * R_i = \sum_a^{nSupport} \left ( \frac{\partial N_a}{\partial X_j} var_{ij}\right ),
   * \f]
   * where $\frac{\partial N_a}{\partial X_j}$ is the basis function gradient,
   *   $var_{ij}$ is the rank-2 symmetric tensor.
   */
  template< int NUM_SUPPORT_POINTS,
            typename GRADIENT_TYPE >
  GEOSX_HOST_DEVICE
  static void gradNajAij( GRADIENT_TYPE const & gradN,
                          real64 const (&var_detJxW)[6],
                          real64 ( &R )[NUM_SUPPORT_POINTS][3] );

  template< int NUM_SUPPORT_POINTS,
            typename GRADIENT_TYPE >
  GEOSX_HOST_DEVICE
  static void gradNajAij( GRADIENT_TYPE const & gradN,
                          real64 const (&var_detJxW)[3][3],
                          real64 ( &R )[NUM_SUPPORT_POINTS][3] );

  template< int NUM_SUPPORT_POINTS >
  GEOSX_HOST_DEVICE
  static void NaFi( real64 const (&N)[NUM_SUPPORT_POINTS],
                    real64 const (&forcingTerm_detJ)[3],
                    real64 ( &R )[NUM_SUPPORT_POINTS][3] );

  template< int NUM_SUPPORT_POINTS,
            typename GRADIENT_TYPE >
  GEOSX_HOST_DEVICE
  static void gradNajAij_plus_NaFi( GRADIENT_TYPE const & gradN,
                                    real64 const (&var_detJxW)[3][3],
                                    real64 const (&N)[NUM_SUPPORT_POINTS],
                                    real64 const (&forcingTerm_detJ)[3],
                                    real64 ( &R )[NUM_SUPPORT_POINTS][3] );

  template< int NUM_SUPPORT_POINTS,
            typename GRADIENT_TYPE >
  GEOSX_HOST_DEVICE
  static void gradNajAij_plus_NaFi( GRADIENT_TYPE const & gradN,
                                    real64 const (&var_detJxW)[6],
                                    real64 const (&N)[NUM_SUPPORT_POINTS],
                                    real64 const (&forcingTerm_detJ)[3],
                                    real64 ( &R )[NUM_SUPPORT_POINTS][3] );


  void setGradNView( arrayView4d< real64 const > const & source )
  {
    m_viewGradN = source;
  }

  void setDetJView( arrayView2d< real64 const > const & source )
  {
    m_viewDetJ = source;
  }

//  void setGradN( localIndex const k,
//                 localIndex const q )

protected:
  arrayView4d< real64 const > m_viewGradN;
  arrayView2d< real64 const > m_viewDetJ;
};

template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void FiniteElementBase::symmetricGradient( GRADIENT_TYPE const & gradN,
                                           real64 const (&var)[NUM_SUPPORT_POINTS][3],
                                           real64 (& gradVar)[6] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    gradVar[0] = gradVar[0] + gradN[a][0] * var[ a ][0];
    gradVar[1] = gradVar[1] + gradN[a][1] * var[ a ][1];
    gradVar[2] = gradVar[2] + gradN[a][2] * var[ a ][2];
    gradVar[3] = gradVar[3] + gradN[a][2] * var[ a ][1] + gradN[a][1] * var[ a ][2];
    gradVar[4] = gradVar[4] + gradN[a][2] * var[ a ][0] + gradN[a][0] * var[ a ][2];
    gradVar[5] = gradVar[5] + gradN[a][1] * var[ a ][0] + gradN[a][0] * var[ a ][1];
  }
}


template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void FiniteElementBase::gradient( GRADIENT_TYPE const & gradN,
                                  real64 const (&var)[NUM_SUPPORT_POINTS][3],
                                  real64 (& gradVar)[3][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        gradVar[i][j] = gradVar[i][j] + var[ a ][i] * gradN[a][j];
      }
    }
  }
}

template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void FiniteElementBase::gradNajAij( GRADIENT_TYPE const & gradN,
                                    real64 const (&var_detJxW)[6],
                                    real64 (& R)[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    R[a][0] = R[a][0] - var_detJxW[0] * gradN[a][0] - var_detJxW[5] * gradN[a][1] - var_detJxW[4] * gradN[a][2];
    R[a][1] = R[a][1] - var_detJxW[5] * gradN[a][0] - var_detJxW[1] * gradN[a][1] - var_detJxW[3] * gradN[a][2];
    R[a][2] = R[a][2] - var_detJxW[4] * gradN[a][0] - var_detJxW[3] * gradN[a][1] - var_detJxW[2] * gradN[a][2];
  }
}


template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void FiniteElementBase::gradNajAij( GRADIENT_TYPE const & gradN,
                                    real64 const (&var_detJxW)[3][3],
                                    real64 (& R)[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    R[a][0] = R[a][0] - var_detJxW[0][0] * gradN[a][0] - var_detJxW[0][1] * gradN[a][1] - var_detJxW[0][2] * gradN[a][2];
    R[a][1] = R[a][1] - var_detJxW[1][0] * gradN[a][0] - var_detJxW[1][1] * gradN[a][1] - var_detJxW[1][2] * gradN[a][2];
    R[a][2] = R[a][2] - var_detJxW[2][0] * gradN[a][0] - var_detJxW[2][1] * gradN[a][1] - var_detJxW[2][2] * gradN[a][2];
  }
}

template< int NUM_SUPPORT_POINTS >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void FiniteElementBase::NaFi( real64 const (&N)[NUM_SUPPORT_POINTS],
                              real64 const (&var_detJxW)[3],
                              real64 ( & R )[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    R[a][0] = R[a][0] + var_detJxW[0] * N[a];
    R[a][1] = R[a][1] + var_detJxW[1] * N[a];
    R[a][2] = R[a][2] + var_detJxW[2] * N[a];
  }
}


template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void FiniteElementBase::gradNajAij_plus_NaFi( GRADIENT_TYPE const & gradN,
                                              real64 const (&var_detJxW)[6],
                                              real64 const (&N)[NUM_SUPPORT_POINTS],
                                              real64 const (&forcingTerm_detJ)[3],
                                              real64 (& R)[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    R[a][0] = R[a][0] - var_detJxW[0] * gradN[a][0] - var_detJxW[5] * gradN[a][1] - var_detJxW[4] * gradN[a][2] + forcingTerm_detJ[0] * N[a];
    R[a][1] = R[a][1] - var_detJxW[5] * gradN[a][0] - var_detJxW[1] * gradN[a][1] - var_detJxW[3] * gradN[a][2] + forcingTerm_detJ[1] * N[a];
    R[a][2] = R[a][2] - var_detJxW[4] * gradN[a][0] - var_detJxW[3] * gradN[a][1] - var_detJxW[2] * gradN[a][2] + forcingTerm_detJ[2] * N[a];
  }
}

template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void FiniteElementBase::gradNajAij_plus_NaFi( GRADIENT_TYPE const & gradN,
                                              real64 const (&var_detJxW)[3][3],
                                              real64 const (&N)[NUM_SUPPORT_POINTS],
                                              real64 const (&forcingTerm_detJ)[3],
                                              real64 (& R)[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    R[a][0] = R[a][0] - var_detJxW[0][0] * gradN[a][0] - var_detJxW[0][1] * gradN[a][1] - var_detJxW[0][2] * gradN[a][2] + forcingTerm_detJ[0] * N[a];
    R[a][1] = R[a][1] - var_detJxW[1][0] * gradN[a][0] - var_detJxW[1][1] * gradN[a][1] - var_detJxW[1][2] * gradN[a][2] + forcingTerm_detJ[1] * N[a];
    R[a][2] = R[a][2] - var_detJxW[2][0] * gradN[a][0] - var_detJxW[2][1] * gradN[a][1] - var_detJxW[2][2] * gradN[a][2] + forcingTerm_detJ[2] * N[a];
  }
}
}
}

#endif //GEOSX_CORE_FINITEELEMENT_FINITEELEMENTBASE
