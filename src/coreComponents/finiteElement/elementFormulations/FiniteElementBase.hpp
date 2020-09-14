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


#ifndef GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_FINITEELEMENTBASE_HPP_
#define GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_FINITEELEMENTBASE_HPP_

#include "common/DataTypes.hpp"
#include "common/GeosxMacros.hpp"


namespace geosx
{
namespace finiteElement
{

/**
 * @brief Base class for FEM element implementations.
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
  static real64 inverse( real64 ( &J )[3][3] );

  /**
   * @brief Calculate the determinant of a 3x3 c-array.
   * @param J The input array.
   * @return The determinant of @p J
   */
  GEOSX_HOST_DEVICE
  static real64 detJ( real64 const (&J)[3][3] );

  /**
   * @brief Get the shape function gradients.
   * @tparam LEAF Type of the derived finite element implementation.
   * @param k The element index.
   * @param q The quadrature point index.
   * @param X Array of coordinates as the reference for the gradients.
   * @param gradN Return array of the shape function gradients.
   * @return The determinant of the Jacobian transformation matrix.
   *
   * This function calls the function to calculate shape function gradients.
   */
  template< typename LEAF >
  GEOSX_HOST_DEVICE
  real64 getGradN( localIndex const k,
                   localIndex const q,
                   real64 const (&X)[LEAF::numNodes][3],
                   real64 ( &gradN )[LEAF::numNodes][3] ) const;


  /**
   * @brief Get the shape function gradients.
   * @tparam LEAF Type of the derived finite element implementation.
   * @param k The element index.
   * @param q The quadrature point index.
   * @param X dummy variable.
   * @param gradN Return array of the shape function gradients.
   * @return The determinant of the Jacobian transformation matrix.
   *
   * This function returns pre-calculated shape function gradients.
   */
  template< typename LEAF >
  GEOSX_HOST_DEVICE
  real64 getGradN( localIndex const k,
                   localIndex const q,
                   int const X,
                   real64 ( &gradN )[LEAF::numNodes][3] ) const;

  /**
   * @name Value Operator Functions
   */
  ///@{

  /**
   * @brief Compute the interpolated value of a variable.
   * @tparam NUM_SUPPORT_POINTS The number of support points for the element.
   * @param N Array (for each support point) of shape function values at the
   *   coordinate the variable is to be interpolated.
   * @param var Array of variable values for each support point.
   * @param value The interpolated value of @p var.
   *
   * This is the standard finite element interpolation operator of a discrete
   * variable defined at the support points.
   * The operator is expressed as:
   * \f[
   * value  = \sum_a^{numSupport} \left ( N_a var_a \right ),
   * \f]

   * @note The shape function values @p N must be evaluated prior to calling this
   * function.
   *
   */
  template< int NUM_SUPPORT_POINTS >
  GEOSX_HOST_DEVICE
  static
  void value( real64 const (&N)[NUM_SUPPORT_POINTS],
              real64 const (&var)[NUM_SUPPORT_POINTS],
              real64 & value );

  /**
   * @brief Compute the interpolated value of a vector variable.
   * @tparam NUM_COMPONENTS Number of components for the vector variable.
   * @copydoc value
   */
  template< int NUM_SUPPORT_POINTS,
            int NUM_COMPONENTS >
  GEOSX_HOST_DEVICE
  static
  void value( real64 const (&N)[NUM_SUPPORT_POINTS],
              real64 const (&var)[NUM_SUPPORT_POINTS][NUM_COMPONENTS],
              real64 ( &value )[NUM_COMPONENTS] );

  ///@}

  /**
   * @name Gradient Operator Functions
   */
  ///@{

  /**
   * @brief Calculate the symmetric gradient of a vector valued support field
   *   at a point using the stored basis function gradients for all support
   *   points.
   * @tparam NUM_SUPPORT_POINTS The number of support points for the element.
   * @tparam GRADIENT_TYPE The type of the array object holding the shape
   * @param gradN The basis function gradients at a point in the element.
   * @param var The vector valued support field that the gradient operator will
   *  be applied to.
   * @param gradVar The symmetric gradient in Voigt notation.
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
   * @brief Calculate the gradient of a scalar valued support field at a point
   *   using the input basis function gradients.
   * @tparam NUM_SUPPORT_POINTS The number of support points for the element.
   * @tparam GRADIENT_TYPE The type of the array object holding the shape
   *   function gradients.
   * @param gradN The basis function gradients at a point in the element.
   * @param var The vector valued support field that the gradient operator will
   *  be applied to.
   * @param gradVar The  gradient.
   *
   * More precisely, the operator is defined as:
   * \f[
   * grad_{j}  = \sum_a^{nSupport} \left ( \frac{\partial N_a}{\partial X_j} var_{a}\right ),
   * \f]
   */
  template< int NUM_SUPPORT_POINTS,
            typename GRADIENT_TYPE >
  GEOSX_HOST_DEVICE
  static void gradient( GRADIENT_TYPE const & gradN,
                        real64 const (&var)[NUM_SUPPORT_POINTS],
                        real64 ( &gradVar )[3] );

  /**
   * @brief Calculate the gradient of a vector valued support field at a point
   *   using the input basis function gradients.
   * @copydoc gradient
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
  ///@}

  /**
   * @name Multi-Operator Functions
   */
  ///@{

  /**
   *
   * @tparam NUM_SUPPORT_POINTS The number of support points for the element.
   * @tparam GRADIENT_TYPE
   * @param N
   * @param gradN
   * @param var
   * @param value
   * @param gradVar
   */
  template< int NUM_SUPPORT_POINTS,
            typename GRADIENT_TYPE >
  GEOSX_HOST_DEVICE
  static void valueAndGradient( real64 const (&N)[NUM_SUPPORT_POINTS],
                                GRADIENT_TYPE const & gradN,
                                real64 const (&var)[NUM_SUPPORT_POINTS],
                                real64 & value,
                                real64 ( &gradVar )[3] );

  ///@}

  /**
   * @name Scattering Operator Functions
   *
   * These functions take quadrature data and map it to the support points
   * through some operator.
   */
  ///@{

  /**
   * @brief Inner product of each basis function gradient with a rank-2
   *   symmetric tensor.
   * @tparam NUM_SUPPORT_POINTS The number of support points for the element.
   * @tparam GRADIENT_TYPE The type of the array object holding the shape
   *   function gradients.
   * @param gradN The basis function gradients at a point in the element.
   * @param var_detJxW The rank-2 tensor at @p q scaled by J*W.
   * @param R The vector at each support point which will hold the result from
   *   the tensor contraction.
   *
   * More precisely, the operator is defined as:
   *
   * \f[
   * R_i = \sum_a^{nSupport} \left( \frac{\partial N_a}{\partial X_j} var_{ij}\right),
   * \f]
   *
   * where \f$\frac{\partial N_a}{\partial X_j}\f$ is the basis function gradient,
   *   \f$var_{ij}\f$ is the rank-2 symmetric tensor.
   */
  template< int NUM_SUPPORT_POINTS,
            typename GRADIENT_TYPE >
  GEOSX_HOST_DEVICE
  static void gradNajAij( GRADIENT_TYPE const & gradN,
                          real64 const (&var_detJxW)[6],
                          real64 ( &R )[NUM_SUPPORT_POINTS][3] );

  /**
   * @copydoc gradNajAij
   * @brief Inner product of each basis function gradient with a rank-2
   *   tensor.
   */
  template< int NUM_SUPPORT_POINTS,
            typename GRADIENT_TYPE >
  GEOSX_HOST_DEVICE
  static void gradNajAij( GRADIENT_TYPE const & gradN,
                          real64 const (&var_detJxW)[3][3],
                          real64 ( &R )[NUM_SUPPORT_POINTS][3] );

  /**
   * @brief Product of each shape function with a vector forcing term.
   * @tparam NUM_SUPPORT_POINTS The number of support points for the element.
   * @param N The shape function value at a predetermined coordinate in the element.
   * @param forcingTerm_detJxW A vector scaled by detJxW
   * @param R The vector at each support point which will hold the result from
   *   the tensor contraction.
   */
  template< int NUM_SUPPORT_POINTS >
  GEOSX_HOST_DEVICE
  static void NaFi( real64 const (&N)[NUM_SUPPORT_POINTS],
                    real64 const (&forcingTerm_detJxW)[3],
                    real64 ( &R )[NUM_SUPPORT_POINTS][3] );

  /**
   * @brief Inner product of each basis function gradient with a rank-2
   *   symmetric tensor added to the product each shape function with a vector.
   * @tparam NUM_SUPPORT_POINTS The number of support points for the element.
   * @tparam GRADIENT_TYPE The type of the array object holding the shape
   *   function gradients.
   * @param gradN The basis function gradients at a point in the element.
   * @param var_detJxW The rank-2 symmetric tensor at @p q scaled by J*W.
   * @param N The shape function value at a predetermined coordinate in the element.
   * @param forcingTerm_detJxW A vector scaled by detJxW
   * @param R The vector at each support point which will hold the result from
   *   the tensor contraction.
   *
   * \f[
   * R_i = \sum_a^{nSupport} \left ( \frac{\partial N_a}{\partial X_j} var_{ij} + N_a f_i \right ),
   * \f]
   */
  template< int NUM_SUPPORT_POINTS,
            typename GRADIENT_TYPE >
  GEOSX_HOST_DEVICE
  static void gradNajAij_plus_NaFi( GRADIENT_TYPE const & gradN,
                                    real64 const (&var_detJxW)[3][3],
                                    real64 const (&N)[NUM_SUPPORT_POINTS],
                                    real64 const (&forcingTerm_detJxW)[3],
                                    real64 ( &R )[NUM_SUPPORT_POINTS][3] );

  /**
   * @brief Inner product of each basis function gradient with a rank-2
   *   tensor added to the product each shape function with a vector.
   * @copydoc gradNajAij_plus_NaFi
   */
  template< int NUM_SUPPORT_POINTS,
            typename GRADIENT_TYPE >
  GEOSX_HOST_DEVICE
  static void gradNajAij_plus_NaFi( GRADIENT_TYPE const & gradN,
                                    real64 const (&var_detJxW)[6],
                                    real64 const (&N)[NUM_SUPPORT_POINTS],
                                    real64 const (&forcingTerm_detJxW)[3],
                                    real64 ( &R )[NUM_SUPPORT_POINTS][3] );


  /**
   * @brief Sets m_viewGradN equal to an input view.
   * @param source The view to assign to m_viewGradN.
   */
  void setGradNView( arrayView4d< real64 const > const & source )
  {
    m_viewGradN = source;
  }

  /**
   * @brief Sets m_viewDetJ equal to an input view.
   * @param source The view to assign to m_viewDetJ.
   */
  void setDetJView( arrayView2d< real64 const > const & source )
  {
    m_viewDetJ = source;
  }

protected:
  /// View to potentially hold pre-calculated shape function gradients.
  arrayView4d< real64 const > m_viewGradN;

  /// View to potentially hold pre-calculated weighted jacobian transformation
  /// determinants.
  arrayView2d< real64 const > m_viewDetJ;
};


/// @cond Doxygen_Suppress

//*************************************************************************************************
//***** Definitons ********************************************************************************
//*************************************************************************************************

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 FiniteElementBase::inverse( real64 (& J)[3][3] )
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

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 FiniteElementBase::detJ( real64 const (&J)[3][3] )
{
  return J[0][0] * ( J[1][1]*J[2][2] - J[1][2]*J[2][1] ) +
         J[1][0] * ( J[0][2]*J[2][1] - J[0][1]*J[2][2] ) +
         J[2][0] * ( J[0][1]*J[1][2] - J[0][2]*J[1][1] );
}

template< typename LEAF >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 FiniteElementBase::getGradN( localIndex const k,
                                    localIndex const q,
                                    real64 const (&X)[LEAF::numNodes][3],
                                    real64 (& gradN)[LEAF::numNodes][3] ) const
{
  GEOSX_UNUSED_VAR( k );
  return LEAF::calcGradN( q, X, gradN );
}

template< typename LEAF >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 FiniteElementBase::getGradN( localIndex const k,
                                    localIndex const q,
                                    int const X,
                                    real64 (& gradN)[LEAF::numNodes][3] ) const
{
  GEOSX_UNUSED_VAR( X );

  for( int a=0; a<LEAF::numNodes; ++a )
  {
    gradN[a][0] = m_viewGradN( k, q, a, 0 );
    gradN[a][1] = m_viewGradN( k, q, a, 1 );
    gradN[a][2] = m_viewGradN( k, q, a, 2 );
  }
  return m_viewDetJ( k, q );
}

//*************************************************************************************************
//***** Interpolated Value Functions **************************************************************
//*************************************************************************************************

template< int NUM_SUPPORT_POINTS >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void FiniteElementBase::value( real64 const (&N)[NUM_SUPPORT_POINTS],
                               real64 const (&var)[NUM_SUPPORT_POINTS],
                               real64 & value )
{
  value = N[0] * var[0];
  for( int a=1; a<NUM_SUPPORT_POINTS; ++a )
  {
    value = value + N[a] * var[a];
  }
}

template< int NUM_SUPPORT_POINTS,
          int NUM_COMPONENTS >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void FiniteElementBase::value( real64 const (&N)[NUM_SUPPORT_POINTS],
                               real64 const (&var)[NUM_SUPPORT_POINTS][NUM_COMPONENTS],
                               real64 (& value)[NUM_COMPONENTS] )
{

  for( int i=0; i<NUM_COMPONENTS; ++i )
  {
    value[i] = N[0] * var[0][i];
  }

  for( int a=1; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int i=0; i<NUM_COMPONENTS; ++i )
    {
      value[i] = value[i] + N[a] * var[a][i];
    }
  }
}


//*************************************************************************************************
//***** Variable Gradient Functions ***************************************************************
//*************************************************************************************************

template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void FiniteElementBase::symmetricGradient( GRADIENT_TYPE const & gradN,
                                           real64 const (&var)[NUM_SUPPORT_POINTS][3],
                                           real64 (& gradVar)[6] )
{
  gradVar[0] = gradN[0][0] * var[0][0];
  gradVar[1] = gradN[0][1] * var[0][1];
  gradVar[2] = gradN[0][2] * var[0][2];
  gradVar[3] = gradN[0][2] * var[0][1] + gradN[0][1] * var[0][2];
  gradVar[4] = gradN[0][2] * var[0][0] + gradN[0][0] * var[0][2];
  gradVar[5] = gradN[0][1] * var[0][0] + gradN[0][0] * var[0][1];

  for( int a=1; a<NUM_SUPPORT_POINTS; ++a )
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
                                  real64 const (&var)[NUM_SUPPORT_POINTS],
                                  real64 (& gradVar)[3] )
{
  for( int i = 0; i < 3; ++i )
  {
    gradVar[i] = var[0] * gradN[0][i];
  }
  for( int a=1; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int i = 0; i < 3; ++i )
    {
      gradVar[i] = gradVar[i] + var[a] * gradN[a][i];
    }
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
  for( int i = 0; i < 3; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      gradVar[i][j] = var[0][i] * gradN[0][j];
    }
  }
  for( int a=1; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        gradVar[i][j] = gradVar[i][j] + var[a][i] * gradN[a][j];
      }
    }
  }
}



template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void FiniteElementBase::valueAndGradient( real64 const (&N)[NUM_SUPPORT_POINTS],
                                          GRADIENT_TYPE const & gradN,
                                          real64 const (&var)[NUM_SUPPORT_POINTS],
                                          real64 & value,
                                          real64 (& gradVar)[3] )
{
  value = N[0] * var[0];
  for( int i = 0; i < 3; ++i )
  {
    gradVar[i] = var[0] * gradN[0][i];
  }

  for( int a=1; a<NUM_SUPPORT_POINTS; ++a )
  {
    value = value + N[a] * var[a];
    for( int i = 0; i < 3; ++i )
    {
      gradVar[i] = gradVar[i] + var[ a ] * gradN[a][i];
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
    R[a][0] = R[a][0] + var_detJxW[0] * gradN[a][0] + var_detJxW[5] * gradN[a][1] + var_detJxW[4] * gradN[a][2];
    R[a][1] = R[a][1] + var_detJxW[5] * gradN[a][0] + var_detJxW[1] * gradN[a][1] + var_detJxW[3] * gradN[a][2];
    R[a][2] = R[a][2] + var_detJxW[4] * gradN[a][0] + var_detJxW[3] * gradN[a][1] + var_detJxW[2] * gradN[a][2];
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
    R[a][0] = R[a][0] + var_detJxW[0][0] * gradN[a][0] + var_detJxW[0][1] * gradN[a][1] + var_detJxW[0][2] * gradN[a][2];
    R[a][1] = R[a][1] + var_detJxW[1][0] * gradN[a][0] + var_detJxW[1][1] * gradN[a][1] + var_detJxW[1][2] * gradN[a][2];
    R[a][2] = R[a][2] + var_detJxW[2][0] * gradN[a][0] + var_detJxW[2][1] * gradN[a][1] + var_detJxW[2][2] * gradN[a][2];
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
                                              real64 const (&forcingTerm_detJxW)[3],
                                              real64 (& R)[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    R[a][0] = R[a][0] + var_detJxW[0] * gradN[a][0] + var_detJxW[5] * gradN[a][1] + var_detJxW[4] * gradN[a][2] + forcingTerm_detJxW[0] * N[a];
    R[a][1] = R[a][1] + var_detJxW[5] * gradN[a][0] + var_detJxW[1] * gradN[a][1] + var_detJxW[3] * gradN[a][2] + forcingTerm_detJxW[1] * N[a];
    R[a][2] = R[a][2] + var_detJxW[4] * gradN[a][0] + var_detJxW[3] * gradN[a][1] + var_detJxW[2] * gradN[a][2] + forcingTerm_detJxW[2] * N[a];
  }
}

template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void FiniteElementBase::gradNajAij_plus_NaFi( GRADIENT_TYPE const & gradN,
                                              real64 const (&var_detJxW)[3][3],
                                              real64 const (&N)[NUM_SUPPORT_POINTS],
                                              real64 const (&forcingTerm_detJxW)[3],
                                              real64 (& R)[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    R[a][0] = R[a][0] + var_detJxW[0][0] * gradN[a][0] + var_detJxW[0][1] * gradN[a][1] + var_detJxW[0][2] * gradN[a][2] + forcingTerm_detJxW[0] * N[a];
    R[a][1] = R[a][1] + var_detJxW[1][0] * gradN[a][0] + var_detJxW[1][1] * gradN[a][1] + var_detJxW[1][2] * gradN[a][2] + forcingTerm_detJxW[1] * N[a];
    R[a][2] = R[a][2] + var_detJxW[2][0] * gradN[a][0] + var_detJxW[2][1] * gradN[a][1] + var_detJxW[2][2] * gradN[a][2] + forcingTerm_detJxW[2] * N[a];
  }
}
/// @endcond

}
}


/// Macro to simplify name resolution in derived classes.
#define USING_FINITEELEMENTBASE                               \
  /** @copydoc FiniteElementBase::value                  */ \
  using FiniteElementBase::value;                           \
  /** @copydoc FiniteElementBase::symmetricGradient      */ \
  using FiniteElementBase::symmetricGradient;               \
  /** @copydoc FiniteElementBase::gradient               */ \
  using FiniteElementBase::gradient;                        \
  /** @copydoc FiniteElementBase::valueAndGradient       */ \
  using FiniteElementBase::valueAndGradient;                \
  /** @copydoc FiniteElementBase::gradNajAij             */ \
  using FiniteElementBase::gradNajAij;                      \
  /** @copydoc FiniteElementBase::NaFi                   */ \
  using FiniteElementBase::NaFi;                            \
  /** @copydoc FiniteElementBase::gradNajAij_plus_NaFi   */ \
  using FiniteElementBase::gradNajAij_plus_NaFi;

#endif //GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_FINITEELEMENTBASE
