

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



  /**
   * @brief Calculate the symmetric gradient of a vector valued support field
   *   at a point using the stored basis function gradients for all support
   *   points.
   * @param dNdX The basis function gradients at a point in the element.
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
            typename BASIS_GRAD_TYPE >
  GEOSX_HOST_DEVICE
  static void symmetricGradient( BASIS_GRAD_TYPE const & dNdX,
                                 real64 const (&var)[NUM_SUPPORT_POINTS][3],
                                 real64 ( &grad )[6] );



  /**
   * @brief Calculate the gradient of a vector valued support field at a point
   *   using the stored basis function gradients for all support points.
   * @param dNdX The basis function gradients at a point in the element.
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
            typename BASIS_GRAD_TYPE >
  GEOSX_HOST_DEVICE
  static void gradient( BASIS_GRAD_TYPE const & dNdX,
                        real64 const (&var)[NUM_SUPPORT_POINTS][3],
                        real64 ( &grad )[3][3] );

  /**
   * @brief Inner product of all basis function gradients and a rank-2
   *   symmetric tensor.
   * @param dNdX The basis function gradients at a point in the element.
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
            typename BASIS_GRAD_TYPE >
  GEOSX_HOST_DEVICE
  static void integrateBasisGradientInnerProduct( BASIS_GRAD_TYPE const & dNdX,
                                                  real64 const (&var_x_detJ_x_W)[6],
                                                  real64 ( &R )[NUM_SUPPORT_POINTS][3] );

  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRAD_TYPE >
  GEOSX_HOST_DEVICE
  static void integrateBasisGradientInnerProduct( BASIS_GRAD_TYPE const & dNdX,
                                                  real64 const (&var_x_detJ_x_W)[3][3],
                                                  real64 ( &R )[NUM_SUPPORT_POINTS][3] );


  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRAD_TYPE >
  GEOSX_HOST_DEVICE
  static void integrateBasisGradientInnerProductPlusForcing( BASIS_GRAD_TYPE const & dNdX,
                                                             real64 const (&var_x_detJ_x_W)[3][3],
                                                             real64 const (&N)[NUM_SUPPORT_POINTS],
                                                             real64 const (&forcingTerm_x_detJ)[3],
                                                             real64 ( &R )[NUM_SUPPORT_POINTS][3] );

  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRAD_TYPE >
  GEOSX_HOST_DEVICE
  static void integrateBasisGradientInnerProductPlusForcing( BASIS_GRAD_TYPE const & dNdX,
                                                             real64 const (&var_x_detJ_x_W)[6],
                                                             real64 const (&N)[NUM_SUPPORT_POINTS],
                                                             real64 const (&forcingTerm_x_detJ)[3],
                                                             real64 ( &R )[NUM_SUPPORT_POINTS][3] );

//TODO we want to keep views and provide interfaces to this data here for cases
//     where we pre-compute the shape function derivatives...maybe...tbd.
//private:
//  arrayView4d< real64 const > const m_dNdX;
//  arrayView2d< real64 const > const m_detJ;


};

template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRAD_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void FiniteElementBase::symmetricGradient( BASIS_GRAD_TYPE const & dNdX,
                                           real64 const (&var)[NUM_SUPPORT_POINTS][3],
                                           real64 (& grad)[6] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    grad[0] = grad[0] + dNdX[a][0] * var[ a ][0];
    grad[1] = grad[1] + dNdX[a][1] * var[ a ][1];
    grad[2] = grad[2] + dNdX[a][2] * var[ a ][2];
    grad[3] = grad[3] + dNdX[a][2] * var[ a ][1] + dNdX[a][1] * var[ a ][2];
    grad[4] = grad[4] + dNdX[a][2] * var[ a ][0] + dNdX[a][0] * var[ a ][2];
    grad[5] = grad[5] + dNdX[a][1] * var[ a ][0] + dNdX[a][0] * var[ a ][1];
  }
}


template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRAD_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void FiniteElementBase::gradient( BASIS_GRAD_TYPE const & dNdX,
                                  real64 const (&var)[NUM_SUPPORT_POINTS][3],
                                  real64 (& grad)[3][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        grad[i][j] = grad[i][j] + var[ a ][i] * dNdX[a][j];
      }
    }
  }
}

template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRAD_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void FiniteElementBase::integrateBasisGradientInnerProduct( BASIS_GRAD_TYPE const & dNdX,
                                                            real64 const (&var_x_detJ_x_W)[6],
                                                            real64 (& R)[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    R[a][0] = R[a][0] - var_x_detJ_x_W[0] * dNdX[a][0] - var_x_detJ_x_W[5] * dNdX[a][1] - var_x_detJ_x_W[4] * dNdX[a][2];
    R[a][1] = R[a][1] - var_x_detJ_x_W[5] * dNdX[a][0] - var_x_detJ_x_W[1] * dNdX[a][1] - var_x_detJ_x_W[3] * dNdX[a][2];
    R[a][2] = R[a][2] - var_x_detJ_x_W[4] * dNdX[a][0] - var_x_detJ_x_W[3] * dNdX[a][1] - var_x_detJ_x_W[2] * dNdX[a][2];
  }
}


template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRAD_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void FiniteElementBase::integrateBasisGradientInnerProduct( BASIS_GRAD_TYPE const & dNdX,
                                                            real64 const (&var_x_detJ_x_W)[3][3],
                                                            real64 (& R)[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    R[a][0] = R[a][0] - var_x_detJ_x_W[0][0] * dNdX[a][0] - var_x_detJ_x_W[0][1] * dNdX[a][1] - var_x_detJ_x_W[0][2] * dNdX[a][2];
    R[a][1] = R[a][1] - var_x_detJ_x_W[1][0] * dNdX[a][0] - var_x_detJ_x_W[1][1] * dNdX[a][1] - var_x_detJ_x_W[1][2] * dNdX[a][2];
    R[a][2] = R[a][2] - var_x_detJ_x_W[2][0] * dNdX[a][0] - var_x_detJ_x_W[2][1] * dNdX[a][1] - var_x_detJ_x_W[2][2] * dNdX[a][2];
  }
}


template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRAD_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void FiniteElementBase::integrateBasisGradientInnerProductPlusForcing( BASIS_GRAD_TYPE const & dNdX,
                                                                       real64 const (&var_x_detJ_x_W)[6],
                                                                       real64 const (&N)[NUM_SUPPORT_POINTS],
                                                                       real64 const (&forcingTerm_x_detJ)[3],
                                                                       real64 (& R)[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    R[a][0] = R[a][0] - var_x_detJ_x_W[0] * dNdX[a][0] - var_x_detJ_x_W[5] * dNdX[a][1] - var_x_detJ_x_W[4] * dNdX[a][2] + forcingTerm_x_detJ[0] * N[a];
    R[a][1] = R[a][1] - var_x_detJ_x_W[5] * dNdX[a][0] - var_x_detJ_x_W[1] * dNdX[a][1] - var_x_detJ_x_W[3] * dNdX[a][2] + forcingTerm_x_detJ[1] * N[a];
    R[a][2] = R[a][2] - var_x_detJ_x_W[4] * dNdX[a][0] - var_x_detJ_x_W[3] * dNdX[a][1] - var_x_detJ_x_W[2] * dNdX[a][2] + forcingTerm_x_detJ[2] * N[a];
  }
}

template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRAD_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void FiniteElementBase::integrateBasisGradientInnerProductPlusForcing( BASIS_GRAD_TYPE const & dNdX,
                                                                       real64 const (&var_x_detJ_x_W)[3][3],
                                                                       real64 const (&N)[NUM_SUPPORT_POINTS],
                                                                       real64 const (&forcingTerm_x_detJ)[3],
                                                                       real64 (& R)[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    R[a][0] = R[a][0] - var_x_detJ_x_W[0][0] * dNdX[a][0] - var_x_detJ_x_W[0][1] * dNdX[a][1] - var_x_detJ_x_W[0][2] * dNdX[a][2] + forcingTerm_x_detJ[0] * N[a];
    R[a][1] = R[a][1] - var_x_detJ_x_W[1][0] * dNdX[a][0] - var_x_detJ_x_W[1][1] * dNdX[a][1] - var_x_detJ_x_W[1][2] * dNdX[a][2] + forcingTerm_x_detJ[1] * N[a];
    R[a][2] = R[a][2] - var_x_detJ_x_W[2][0] * dNdX[a][0] - var_x_detJ_x_W[2][1] * dNdX[a][1] - var_x_detJ_x_W[2][2] * dNdX[a][2] + forcingTerm_x_detJ[2] * N[a];
  }
}
}
}

#endif //GEOSX_CORE_FINITEELEMENT_FINITEELEMENTBASE
