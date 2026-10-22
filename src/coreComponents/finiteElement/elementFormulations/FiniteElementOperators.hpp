#pragma once

#include "common/DataTypes.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{
namespace finiteElement
{
namespace feOps
{

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
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void value( real64 const (&N)[NUM_SUPPORT_POINTS],
            real64 const (&var)[NUM_SUPPORT_POINTS],
            real64 & value )
{
  value = LvArray::tensorOps::AiBi< NUM_SUPPORT_POINTS >( N, var );
}

/**
 * @brief Compute the interpolated value of a vector variable.
 * @tparam NUM_COMPONENTS Number of components for the vector variable.
 * @copydoc value
 */
template< int NUM_SUPPORT_POINTS,
          int NUM_COMPONENTS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void value( real64 const (&N)[NUM_SUPPORT_POINTS],
            real64 const (&var)[NUM_SUPPORT_POINTS][NUM_COMPONENTS],
            real64 (& value)[NUM_COMPONENTS] )
{
  LvArray::tensorOps::Ri_eq_AjiBj< NUM_COMPONENTS, NUM_SUPPORT_POINTS >( value, var, N );
}

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
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void symmetricGradient( GRADIENT_TYPE const & gradN,
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

/**
 * @brief Calculate the trace of the symmetric gradient of a vector valued support
 *   field (i.e. the volumetric strain for the displacement field) at a point using
 *   the stored basis function gradients for all support points.
 * @tparam NUM_SUPPORT_POINTS The number of support points for the element.
 * @tparam GRADIENT_TYPE The type of the array object holding the shape
 * @param gradN The basis function gradients at a point in the element.
 * @param var The vector valued support field that the gradient operator will
 *  be applied to.
 * @return The trace of the symetric gradient tensor.
 *
 */
template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 symmetricGradientTrace( GRADIENT_TYPE const & gradN,
                               real64 const (&var)[NUM_SUPPORT_POINTS][3] )
{
  real64 result = gradN[0][0] * var[0][0] + gradN[0][1] * var[0][1] + gradN[0][2] * var[0][2];

  for( int a=1; a<NUM_SUPPORT_POINTS; ++a )
  {
    result = result + gradN[a][0] * var[a][0] + gradN[a][1] * var[a][1] + gradN[a][2] * var[a][2];
  }
  return result;
}

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
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void gradient( GRADIENT_TYPE const & gradN,
               real64 const (&var)[NUM_SUPPORT_POINTS],
               real64 (& gradVar)[3] )
{
  LvArray::tensorOps::Ri_eq_AjiBj< 3, NUM_SUPPORT_POINTS >( gradVar, gradN, var );
}

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
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void gradient( GRADIENT_TYPE const & gradN,
               real64 const (&var)[NUM_SUPPORT_POINTS][3],
               real64 (& gradVar)[3][3] )
{
  LvArray::tensorOps::Rij_eq_AkiBkj< 3, 3, NUM_SUPPORT_POINTS >( gradVar, var, gradN );
}

///@}

/**
 * @name Multi-Operator Functions
 */
///@{

/**
 * @brief Calculate the value and gradient of a scalar valued support field
 *   at a point using the input basis function gradients.
 * @tparam NUM_SUPPORT_POINTS The number of support points for the element.
 * @tparam GRADIENT_TYPE The type of the array object holding the shape
 * @param N Array (for each support point) of shape function values at the
 *   coordinate the variable is to be interpolated.
 * @param gradN The basis function gradients at a point in the element.
 * @param var The vector valued support field that the gradient operator will
 *  be applied to.
 * @param value The value at the point for which N was specified.
 * @param gradVar The gradient at the point for which gradN was specified.
 */
template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void valueAndGradient( real64 const (&N)[NUM_SUPPORT_POINTS],
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
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void plusGradNajAij( GRADIENT_TYPE const & gradN,
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

/**
 * @copydoc plusGradNajAij
 * @brief Inner product of each basis function gradient with a rank-2
 *   tensor.
 */
template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void plusGradNajAij( GRADIENT_TYPE const & gradN,
                     real64 const (&var_detJxW)[3][3],
                     real64 (& R)[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    LvArray::tensorOps::Ri_add_AijBj< 3, 3 >( R[a], var_detJxW, gradN[a] );
  }
}

/**
 * @brief Product of each shape function with a vector forcing term.
 * @tparam NUM_SUPPORT_POINTS The number of support points for the element.
 * @param N The shape function value at a predetermined coordinate in the element.
 * @param forcingTerm_detJxW A vector scaled by detJxW
 * @param R The vector at each support point which will hold the result from
 *   the tensor contraction.
 */
template< int NUM_SUPPORT_POINTS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void plusNaFi( real64 const (&N)[NUM_SUPPORT_POINTS],
               real64 const (&var_detJxW)[3],
               real64 ( & R )[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    LvArray::tensorOps::scaledAdd< 3 >( R[a], var_detJxW, N[a] );
  }
}

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
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void plusGradNajAijPlusNaFi( GRADIENT_TYPE const & gradN,
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

/**
 * @brief Inner product of each basis function gradient with a rank-2
 *   tensor added to the product each shape function with a vector.
 * @copydoc plusGradNajAijPlusNaFi
 */
template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void plusGradNajAijPlusNaFi( GRADIENT_TYPE const & gradN,
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

///@}

} // namespace finiteElement::finiteElement::feOps
} // namespace geos::finiteElement
} // namespace geos
