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

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 determinant( real64 const (& J)[3][3] )
{
  real64 const detJ = J[0][0] * (J[1][1] * J[2][2] - J[2][1] * J[1][2])
                    - J[1][0] * (J[0][1] * J[2][2] - J[2][1] * J[0][2])
                    + J[2][0] * (J[0][1] * J[1][2] - J[1][1] * J[0][2]);
  return detJ;
}

// J: phys_dim x ref_dim, AdjJ: ref_dim x phys_dim
// J^-1 := detJinv * AdjJ.
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

// J: phys_dim x ref_dim, Jinv: ref_dim x phys_dim
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 getInverseAndDeterminant( real64 const (& J)[3][3], real64 (& Jinv)[3][3] )
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

// Compute the gradient of a vectorial field in physical coordinates
// from the gradient in reference coordinates. Using the formula:
//   grad_phys := J^-1 * grad.
template < localIndex ref_dim,
           localIndex phys_dim,
           localIndex num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void computePhysicalGradient( real64 const (& Jinv)[ref_dim][phys_dim],
                              real64 const (& grad)[num_comp][ref_dim],
                              real64 (& grad_phys)[num_comp][phys_dim] )
{
  for (localIndex c = 0; c < num_comp; c++)
  {
    for (localIndex i = 0; i < phys_dim; i++)
    {
      real64 val = 0.0;
      for (localIndex j = 0; j < ref_dim; j++)
      {
        val = val + Jinv[ j ][ i ] * grad[ c ][ j ];
      }
      grad_phys[ c ][ i ] = val;
    }
  }
}

// Compute the gradient of a vectorial field in physical coordinates
// from the gradient in reference coordinates. Using the formula:
//   grad_phys := J^-1 * grad,
// with
//   J^-1 := detJinv * AdjJ.
template < localIndex ref_dim,
           localIndex phys_dim,
           localIndex num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void computePhysicalGradient( real64 const detJinv,
                              real64 const (& AdjJ)[ref_dim][phys_dim],
                              real64 const (& grad)[num_comp][ref_dim],
                              real64 (& grad_phys)[num_comp][phys_dim] )
{
  for (localIndex c = 0; c < num_comp; c++)
  {
    for (localIndex i = 0; i < phys_dim; i++)
    {
      real64 val = 0.0;
      for (localIndex j = 0; j < ref_dim; j++)
      {
        val = val + AdjJ[ j ][ i ] * grad[ c ][ j ];
      }
      grad_phys[ c ][ i ] = detJinv * val;
    }
  }
}

/// Compute symmetric gradient
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

/// Compute symm_strain = 0.5 * ( grad_phys + grad_phys^T )
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

// } // namespace finiteElement
} // namespace geosx



#endif /* GEOSX_FINITEELEMENT_QUADRATUREFUNCTIONSHELPER_HPP_ */
