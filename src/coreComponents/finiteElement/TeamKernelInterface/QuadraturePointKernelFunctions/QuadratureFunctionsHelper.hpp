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
 * @file QuadratureFunctionsHelper.hpp
 */

#ifndef GEOSX_FINITEELEMENT_QUADRATUREFUNCTIONSHELPER_HPP_
#define GEOSX_FINITEELEMENT_QUADRATUREFUNCTIONSHELPER_HPP_

#include "common/DataTypes.hpp"
#include "common/GeosxMacros.hpp"

namespace geosx
{

/**
 * @namespace finiteElement Contains the finite element implementation.
 */
// namespace finiteElement
// {

/**
 * @brief Compute the determinant of a 3x3 matrix.
 * @param J The input 3x3 matrix.
 * @return The determinant of the input matrix.
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 determinant( real64 const (& J)[3][3] )
{
  real64 const detJ = J[0][0] * (J[1][1] * J[2][2] - J[2][1] * J[1][2])
                    - J[1][0] * (J[0][1] * J[2][2] - J[2][1] * J[0][2])
                    + J[2][0] * (J[0][1] * J[1][2] - J[1][1] * J[0][2]);
  return detJ;
}

/**
 * @brief Compute the adjugate matrix of a given matrix.
 * @param J The input matrix.
 * @param AdjJ The output matrix.
 * 
 * J : phys_dim x ref_dim, AdjJ: ref_dim x phys_dim
 * J^-1 := detJinv * AdjJ.
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void adjugate( real64 const (& J)[3][3], real64 (& AdjJ)[3][3] )
{
  AdjJ[0][0] = (J[1][1] * J[2][2]) - (J[1][2] * J[2][1]);
  AdjJ[0][1] = (J[2][1] * J[0][2]) - (J[0][1] * J[2][2]);
  AdjJ[0][2] = (J[0][1] * J[1][2]) - (J[1][1] * J[0][2]);
  AdjJ[1][0] = (J[2][0] * J[1][2]) - (J[1][0] * J[2][2]);
  AdjJ[1][1] = (J[0][0] * J[2][2]) - (J[0][2] * J[2][0]);
  AdjJ[1][2] = (J[1][0] * J[0][2]) - (J[0][0] * J[1][2]);
  AdjJ[2][0] = (J[1][0] * J[2][1]) - (J[2][0] * J[1][1]);
  AdjJ[2][1] = (J[2][0] * J[0][1]) - (J[0][0] * J[2][1]);
  AdjJ[2][2] = (J[0][0] * J[1][1]) - (J[0][1] * J[1][0]);
}

/**
 * @brief Compute the inverse and the determinant of a given 3x3 matrix.
 * @param J The input matrix.
 * @param Jinv The inverse matrix of the input matrix.
 * @return The determinant of the input matrix.
 * 
 * J : phys_dim x ref_dim, Jinv: ref_dim x phys_dim
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 computeInverseAndDeterminant( real64 const (& J)[3][3], real64 (& Jinv)[3][3] )
{
  real64 const detJ = J[0][0] * (J[1][1] * J[2][2] - J[2][1] * J[1][2])
                    - J[1][0] * (J[0][1] * J[2][2] - J[2][1] * J[0][2])
                    + J[2][0] * (J[0][1] * J[1][2] - J[1][1] * J[0][2]);
  real64 const detJinv = 1.0 / detJ;
  Jinv[0][0] = detJinv * ( (J[1][1] * J[2][2]) - (J[1][2] * J[2][1]) );
  Jinv[0][1] = detJinv * ( (J[2][1] * J[0][2]) - (J[0][1] * J[2][2]) );
  Jinv[0][2] = detJinv * ( (J[0][1] * J[1][2]) - (J[1][1] * J[0][2]) );
  Jinv[1][0] = detJinv * ( (J[2][0] * J[1][2]) - (J[1][0] * J[2][2]) );
  Jinv[1][1] = detJinv * ( (J[0][0] * J[2][2]) - (J[0][2] * J[2][0]) );
  Jinv[1][2] = detJinv * ( (J[1][0] * J[0][2]) - (J[0][0] * J[1][2]) );
  Jinv[2][0] = detJinv * ( (J[1][0] * J[2][1]) - (J[2][0] * J[1][1]) );
  Jinv[2][1] = detJinv * ( (J[2][0] * J[0][1]) - (J[0][0] * J[2][1]) );
  Jinv[2][2] = detJinv * ( (J[0][0] * J[1][1]) - (J[0][1] * J[1][0]) );
  return detJ;
}

/**
 * @brief Compute the gradient of a vectorial field in physical coordinates
 * from the gradient in reference coordinates. Using the formula:
 *   grad_phys := J^-1 * grad_ref.
 * @tparam ref_dim The dimension of the reference coordinate space.
 * @tparam phys_dim The dimension of the physical coordinate space.
 * @param Jinv The inverse of the mesh jacobian matrix.
 * @param grad_ref The gradient of the field in reference coordinates.
 * @param grad_phys The gradient of the field in physical coordinates.
 */
template < localIndex ref_dim,
           localIndex phys_dim >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void computePhysicalGradient( real64 const (& Jinv)[ref_dim][phys_dim],
                              real64 const (& grad_ref)[ref_dim],
                              real64 (& grad_phys)[phys_dim] )
{
  for (localIndex i = 0; i < phys_dim; i++)
  {
    real64 val = 0.0;
    for (localIndex j = 0; j < ref_dim; j++)
    {
      val = val + Jinv[ j ][ i ] * grad_ref[ j ];
    }
    grad_phys[ i ] = val;
  }
}

/**
 * @brief Compute the gradient of a vectorial field in physical coordinates
 * from the gradient in reference coordinates. Using the formula:
 *   grad_phys := J^-1 * grad_ref.
 * @tparam ref_dim The dimension of the reference coordinate space.
 * @tparam phys_dim The dimension of the physical coordinate space.
 * @tparam num_comp The number of components of the field.
 * @param Jinv The inverse of the mesh jacobian matrix.
 * @param grad_ref The gradient of the field in reference coordinates.
 * @param grad_phys The gradient of the field in physical coordinates.
 */
template < localIndex ref_dim,
           localIndex phys_dim,
           localIndex num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void computePhysicalGradient( real64 const (& Jinv)[ref_dim][phys_dim],
                              real64 const (& grad_ref)[num_comp][ref_dim],
                              real64 (& grad_phys)[num_comp][phys_dim] )
{
  for (localIndex c = 0; c < num_comp; c++)
  {
    for (localIndex i = 0; i < phys_dim; i++)
    {
      real64 val = 0.0;
      for (localIndex j = 0; j < ref_dim; j++)
      {
        val = val + Jinv[ j ][ i ] * grad_ref[ c ][ j ];
      }
      grad_phys[ c ][ i ] = val;
    }
  }
}

/**
 * @brief Compute the gradient of a vectorial field in physical coordinates
 * from the gradient in reference coordinates. Using the formula:
 *   grad_phys := J^-1 * grad,
 * with
 *   J^-1 := detJinv * AdjJ.
 * @tparam ref_dim The dimension of the reference coordinate space.
 * @tparam phys_dim The dimension of the physical coordinate space.
 * @param detJinv The inverse of the determinant of the mesh jacobian matrix.
 * @param AdjJ The adjugate of the mesh jacobian matrix.
 * @param grad_ref The gradient of the field in reference coordinates.
 * @param grad_phys The gradient of the field in physical coordinates.
 */ 
template < localIndex ref_dim,
           localIndex phys_dim >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void computePhysicalGradient( real64 const detJinv,
                              real64 const (& AdjJ)[ref_dim][phys_dim],
                              real64 const (& grad_ref)[ref_dim],
                              real64 (& grad_phys)[phys_dim] )
{
  for (localIndex i = 0; i < phys_dim; i++)
  {
    real64 val = 0.0;
    for (localIndex j = 0; j < ref_dim; j++)
    {
      val = val + AdjJ[ j ][ i ] * grad_ref[ j ];
    }
    grad_phys[ i ] = detJinv * val;
  }
}

/**
 * @brief Compute the gradient of a vectorial field in physical coordinates
 * from the gradient in reference coordinates. Using the formula:
 *   grad_phys := J^-1 * grad,
 * with
 *   J^-1 := detJinv * AdjJ.
 * @tparam ref_dim The dimension of the reference coordinate space.
 * @tparam phys_dim The dimension of the physical coordinate space.
 * @tparam num_comp The number of components of the field.
 * @param detJinv The inverse of the determinant of the mesh jacobian matrix.
 * @param AdjJ The adjugate of the mesh jacobian matrix.
 * @param grad_ref The gradient of the field in reference coordinates.
 * @param grad_phys The gradient of the field in physical coordinates.
 */ 
template < localIndex ref_dim,
           localIndex phys_dim,
           localIndex num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void computePhysicalGradient( real64 const detJinv,
                              real64 const (& AdjJ)[ref_dim][phys_dim],
                              real64 const (& grad_ref)[num_comp][ref_dim],
                              real64 (& grad_phys)[num_comp][phys_dim] )
{
  for (localIndex c = 0; c < num_comp; c++)
  {
    for (localIndex i = 0; i < phys_dim; i++)
    {
      real64 val = 0.0;
      for (localIndex j = 0; j < ref_dim; j++)
      {
        val = val + AdjJ[ j ][ i ] * grad_ref[ c ][ j ];
      }
      grad_phys[ c ][ i ] = detJinv * val;
    }
  }
}

/**
 * @brief Compute the gradient of a vectorial field in reference coordinates
 * from the gradient in physical coordinates. Using the formula:
 *   grad_ref := J^-T * grad_phys.
 * @tparam ref_dim The dimension of the reference coordinate space.
 * @tparam phys_dim The dimension of the physical coordinate space.
 * @param Jinv The inverse of the mesh jacobian matrix.
 * @param grad_phys The gradient of the field in physical coordinates.
 * @param grad_ref The gradient of the field in reference coordinates. 
 */
template < localIndex ref_dim,
           localIndex phys_dim >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void computeReferenceGradient( real64 const (& Jinv)[ref_dim][phys_dim],
                               real64 const (& grad_phys)[phys_dim],
                               real64 (& grad_ref)[ref_dim] )
{
  for (localIndex i = 0; i < ref_dim; i++)
  {
    real64 val = 0.0;
    for (localIndex j = 0; j < phys_dim; j++)
    {
      val = val + Jinv[ i ][ j ] * grad_phys[ j ];
    }
    grad_ref[ i ] = val;
  }
}

/**
 * @brief Compute the gradient of a vectorial field in reference coordinates
 * from the gradient in physical coordinates. Using the formula:
 *   grad_ref := J^-T * grad_phys.
 * @tparam ref_dim The dimension of the reference coordinate space.
 * @tparam phys_dim The dimension of the physical coordinate space.
 * @tparam num_comp The number of components of the field.
 * @param Jinv The inverse of the mesh jacobian matrix.
 * @param grad_phys The gradient of the field in physical coordinates.
 * @param grad_ref The gradient of the field in reference coordinates. 
 */
template < localIndex ref_dim,
           localIndex phys_dim,
           localIndex num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void computeReferenceGradient( real64 const (& Jinv)[ref_dim][phys_dim],
                               real64 const (& grad_phys)[num_comp][phys_dim],
                               real64 (& grad_ref)[num_comp][ref_dim] )
{
  for (localIndex c = 0; c < num_comp; c++)
  {
    for (localIndex i = 0; i < ref_dim; i++)
    {
      real64 val = 0.0;
      for (localIndex j = 0; j < phys_dim; j++)
      {
        val = val + Jinv[ i ][ j ] * grad_phys[ c ][ j ];
      }
      grad_ref[ c ][ i ] = val;
    }
  }
}

/**
 * @brief Compute the gradient of a vectorial field in reference coordinates
 * from the gradient in physical coordinates. Using the formula:
 *   grad_phys := J^-T * grad,
 * with
 *   J^-1 := detJinv * AdjJ.
 * @tparam ref_dim The dimension of the reference coordinate space.
 * @tparam phys_dim The dimension of the physical coordinate space.
 * @param detJinv The inverse of the determinant of the mesh jacobian matrix.
 * @param AdjJ The adjugate of the mesh jacobian matrix.
 * @param grad_phys The gradient of the field in physical coordinates.
 * @param grad_ref The gradient of the field in reference coordinates. 
 */
template < localIndex ref_dim,
           localIndex phys_dim >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void computeReferenceGradient( real64 const detJinv,
                               real64 const (& AdjJ)[ref_dim][phys_dim],
                               real64 const (& grad_phys)[phys_dim],
                               real64 (& grad_ref)[ref_dim] )
{
  for (localIndex i = 0; i < ref_dim; i++)
  {
    real64 val = 0.0;
    for (localIndex j = 0; j < phys_dim; j++)
    {
      val = val + AdjJ[ i ][ j ] * grad_phys[ j ];
    }
    grad_ref[ i ] = detJinv * val;
  }
}

/**
 * @brief Compute the gradient of a vectorial field in reference coordinates
 * from the gradient in physical coordinates. Using the formula:
 *   grad_phys := J^-T * grad,
 * with
 *   J^-1 := detJinv * AdjJ.
 * @tparam ref_dim The dimension of the reference coordinate space.
 * @tparam phys_dim The dimension of the physical coordinate space.
 * @tparam num_comp The number of components of the field.
 * @param detJinv The inverse of the determinant of the mesh jacobian matrix.
 * @param AdjJ The adjugate of the mesh jacobian matrix.
 * @param grad_phys The gradient of the field in physical coordinates.
 * @param grad_ref The gradient of the field in reference coordinates. 
 */
template < localIndex ref_dim,
           localIndex phys_dim,
           localIndex num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void computeReferenceGradient( real64 const detJinv,
                               real64 const (& AdjJ)[ref_dim][phys_dim],
                               real64 const (& grad_phys)[num_comp][phys_dim],
                               real64 (& grad_ref)[num_comp][ref_dim] )
{
  for (localIndex c = 0; c < num_comp; c++)
  {
    for (localIndex i = 0; i < ref_dim; i++)
    {
      real64 val = 0.0;
      for (localIndex j = 0; j < phys_dim; j++)
      {
        val = val + AdjJ[ i ][ j ] * grad_phys[ c ][ j ];
      }
      grad_ref[ c ][ i ] = detJinv * val;
    }
  }
}

/**
 * @brief Apply the quadrature weights for numerical integration.
 * 
 * @tparam dim The dimension of the space.
 * @param weight The quadrature point weight.
 * @param detJ The determinant of the mesh jacobian.
 * @param field The integrated field.
 */
template < localIndex dim >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void applyQuadratureWeights( real64 const weight, real64 const detJ, real64 (& field)[dim] )
{
  const real64 w = - weight * detJ;
  for (localIndex i = 0; i < dim; i++)
  {
    field[i] = w * field[i];
  }
}

/**
 * @brief Apply the quadrature weights for numerical integration.
 * 
 * @tparam num_comp The number of components of the field.
 * @tparam dim The diemsnion of the space.
 * @param weight The quadrature point weight.
 * @param detJ The determinant of the mesh jacobian.
 * @param field The integrated field.
 */
template < localIndex num_comp, localIndex dim >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void applyQuadratureWeights( real64 const weight, real64 const detJ, real64 (& field)[num_comp][dim] )
{
  const real64 w = - weight * detJ;
  for (localIndex j = 0; j < dim; j++)
  {
    for (localIndex i = 0; i < num_comp; i++)
    {
      field[i][j] = w * field[i][j];
    }    
  }  
}

/**
 * @brief Compute the symmetric gradient.
 * 
 * @param grad_phys The input gradient.
 * @param symm_grad The symmetric gradient.
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void computeSymmetricGradient( real64 const (& grad_phys)[3][3], real64 (& symm_grad)[6] )
{
  symm_grad[0] = grad_phys[0][0];
  symm_grad[1] = grad_phys[1][1];
  symm_grad[2] = grad_phys[2][2];
  symm_grad[5] = grad_phys[0][1] + grad_phys[1][0];
  symm_grad[4] = grad_phys[0][2] + grad_phys[2][0];
  symm_grad[3] = grad_phys[1][2] + grad_phys[2][1];
}

/**
 * @brief Compute the symmetric strain:
 *  symm_strain = 0.5 * ( grad_phys + grad_phys^T ).
 * 
 * @param grad_phys The gradient in physical coordinates.
 * @param symm_strain The strain.
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void computeStrain( real64 const (& grad_phys)[3][3], real64 (& symm_strain)[6] )
{
  symm_strain[0] = grad_phys[0][0];
  symm_strain[1] = grad_phys[1][1];
  symm_strain[2] = grad_phys[2][2];
  symm_strain[5] = 0.5 * ( grad_phys[0][1] + grad_phys[1][0] );
  symm_strain[4] = 0.5 * ( grad_phys[0][2] + grad_phys[2][0] );
  symm_strain[3] = 0.5 * ( grad_phys[1][2] + grad_phys[2][1] );
}

/**
 * @brief Compute the stress from a constitutive model.
 * 
 * @tparam Stack The type of the kernel stack variables.
 * @tparam ConstitutiveModel The type of the constitutive model.
 * @param stack The kernel stack variables.
 * @param constitutiveUpdate The constitutive model.
 * @param q The "linear" index of the quadrature point.
 * @param grad_phys The gradient in physical coordinates.
 * @param stress The stress at the given quadrature point.
 */
template < typename Stack, typename ConstitutiveModel >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void computeStress( Stack & stack,
                    ConstitutiveModel & constitutiveUpdate,
                    localIndex const q,
                    real64 const (& grad_phys)[ 3 ][ 3 ],
                    real64 (& stress)[ 3 ][ 3 ] )
{
  real64 symm_strain[ 6 ];
  computeSymmetricGradient( grad_phys, symm_strain );
  real64 symm_stress[ 6 ];
  constitutiveUpdate.smallStrainNoStateUpdate_StressOnly( stack.element_index, q, symm_strain, symm_stress );
  stress[ 0 ][ 0 ] = symm_stress[ 0 ];
  stress[ 1 ][ 1 ] = symm_stress[ 1 ];
  stress[ 2 ][ 2 ] = symm_stress[ 2 ];
  stress[ 0 ][ 1 ] = symm_stress[ 5 ];
  stress[ 0 ][ 2 ] = symm_stress[ 4 ];
  stress[ 1 ][ 2 ] = symm_stress[ 3 ];
  stress[ 1 ][ 0 ] = symm_stress[ 5 ];
  stress[ 2 ][ 0 ] = symm_stress[ 4 ];
  stress[ 2 ][ 1 ] = symm_stress[ 3 ];
}

/**
 * @brief Compute the stress from a constitutive model.
 * 
 * @tparam Stack The type of the kernel stack variables.
 * @tparam ConstitutiveModel The type of the constitutive model.
 * @param stack The kernel stack variables.
 * @param constitutiveUpdate The constitutive model.
 * @param q The "linear" index of the quadrature point.
 * @param Jinv The inverse of the mesh jacobian matrix.
 * @param grad_ref The gradient in reference coordinates.
 * @param stress The stress at the given quadrature point.
 */
template < typename Stack, typename ConstitutiveModel >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void computeStress( Stack & stack,
                    ConstitutiveModel & constitutiveUpdate,
                    localIndex const q,
                    real64 const (& Jinv)[ 3 ][ 3 ],
                    real64 const (& grad_ref)[ 3 ][ 3 ],
                    real64 (& stress)[ 3 ][ 3 ] )
{
  real64 grad_phys[ 3 ][ 3 ];
  computePhysicalGradient( Jinv, grad_ref, grad_phys );
  computeStress( stack, constitutiveUpdate, q, grad_phys, stress );
}

// } // namespace finiteElement
} // namespace geosx



#endif /* GEOSX_FINITEELEMENT_QUADRATUREFUNCTIONSHELPER_HPP_ */
