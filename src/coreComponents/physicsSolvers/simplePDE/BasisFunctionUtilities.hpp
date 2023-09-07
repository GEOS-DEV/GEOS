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
 * @file QuadratureUtilities.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SIMPLEPDE_BASISFUNCTIONUTILITIES_HPP_
#define GEOS_PHYSICSSOLVERS_SIMPLEPDE_BASISFUNCTIONUTILITIES_HPP_

#include "CellUtilities.hpp"

namespace geos
{

namespace BasisFunctionUtilities
{

/**
 * @brief Basis function type
 */
enum class BasisFunction : integer
{
  Lagrange
};

struct Gradient
{
  real64 val[3];
};

template< typename CELL_TYPE,
          BasisFunction BASIS_FUNCTION >
struct Helper
{};

template<>
struct Helper< HexahedronCell,
               BasisFunction::Lagrange >
{
  GEOS_HOST_DEVICE
  static Gradient getGradient( int const BasisIndex,
                                            real64 const Xiq[3] )
  {
    real64 dPhiLin[2] = { -1.0, 1.0 };

    int const a = BasisIndex & 1;
    int const b = ( BasisIndex & 2 ) >> 1;
    int const c = ( BasisIndex & 4 ) >> 2;

    Gradient gradient;
    gradient.val[0] = 0.125 * (       dPhiLin[a]          ) * ( 1.0 + dPhiLin[b] * Xiq[1] ) * ( 1.0 + dPhiLin[c] * Xiq[2] );
    gradient.val[1] = 0.125 * ( 1.0 + dPhiLin[a] * Xiq[0] ) * (       dPhiLin[b]          ) * ( 1.0 + dPhiLin[c] * Xiq[2] );
    gradient.val[2] = 0.125 * ( 1.0 + dPhiLin[a] * Xiq[0] ) * ( 1.0 + dPhiLin[b] * Xiq[1] ) * (       dPhiLin[c]          );
    return gradient;
  }
};
                  


// getQuadratureData< CELL_TYPE, INTEGRATION_RULE, INTEGRATION_ORDER >

template< typename CELL_TYPE,
          BasisFunction BASIS_FUNCTION >
GEOS_HOST_DEVICE
static Gradient getGradient( int const BasisIndex,
                             real64 const Xiq[3] )
{
  return Helper< CELL_TYPE, BASIS_FUNCTION >::getGradient( BasisIndex, Xiq );
}

/// Declare strings associated with enumeration values.
ENUM_STRINGS( BasisFunction,
              "Lagrange" );

} // namespace BasisFunctionUtilities

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_SIMPLEPDE_BASISFUNCTIONUTILITIES_HPP_