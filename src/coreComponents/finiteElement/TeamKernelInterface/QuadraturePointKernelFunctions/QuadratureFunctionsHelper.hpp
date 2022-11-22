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

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void computePhysicalGradient( real64 const detJinv,
                              real64 const (& AdjJ)[3][3],
                              real64 const (& grad)[3][3],
                              real64 (& grad_phys)[3][3] )
{
  constexpr localIndex dim = 3;
  for (localIndex i = 0; i < dim; i++)
  {
    for (localIndex j = 0; j < dim; j++)
    {
      real64 val = 0.0;
      for (localIndex k = 0; k < dim; k++)
      {
        val = val + AdjJ[ i ][ k ] * grad[ k ][ j ];
      }
      grad_phys[ i ][ j ] = val * detJinv;
    }
  }
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void computeStrain( real64 const (& grad_phys)[3][3], real64 (& symm_strain)[6] )
{
  symm_strain[0] = grad_phys[0][0];
  symm_strain[1] = grad_phys[1][1];
  symm_strain[2] = grad_phys[2][2];
  symm_strain[3] = grad_phys[0][1] + grad_phys[1][0];
  symm_strain[4] = grad_phys[0][2] + grad_phys[2][0];
  symm_strain[5] = grad_phys[1][2] + grad_phys[2][1];
}

// } // namespace finiteElement
} // namespace geosx



#endif /* GEOSX_FINITEELEMENT_QUADRATUREFUNCTIONSHELPER_HPP_ */
