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
  real64 data[3];
};

template< CellUtilities::ReferenceCell REFERENCE_CELL_TYPE,
          BasisFunction BASIS_FUNCTION,
          int NUM_SUPPORT_POINTS >
struct Helper
{};

template<>
struct Helper< CellUtilities::ReferenceCell::Cube,
               BasisFunction::Lagrange,
               8 >
{
  GEOS_HOST_DEVICE
  static Gradient getParentGradient( int const BasisIndex,
                                     real64 const (&Xiq)[3] )
  {
    constexpr real64 dPhiLin[2] = { -1.0, 1.0 };

    int const a = BasisIndex & 1;
    int const b = ( BasisIndex & 2 ) >> 1;
    int const c = ( BasisIndex & 4 ) >> 2;

    Gradient gradient;
    gradient.data[0] = 0.125 * (       dPhiLin[a]          ) * ( 1.0 + dPhiLin[b] * Xiq[1] ) * ( 1.0 + dPhiLin[c] * Xiq[2] );
    gradient.data[1] = 0.125 * ( 1.0 + dPhiLin[a] * Xiq[0] ) * (       dPhiLin[b]          ) * ( 1.0 + dPhiLin[c] * Xiq[2] );
    gradient.data[2] = 0.125 * ( 1.0 + dPhiLin[a] * Xiq[0] ) * ( 1.0 + dPhiLin[b] * Xiq[1] ) * (       dPhiLin[c]          );
    return gradient;
  }

  GEOS_HOST_DEVICE
  static void getGradient( real64 const (&Xiq)[3],
                           typename CuboidCell::JacobianType const & Jinv,
                           real64 (& dNdX)[8][3] )
  {
    for( int i = 0; i < 8; ++i )
    {
      // ... ... Parent space
      Gradient dNdXi = getParentGradient( i, Xiq );

      // ... ... Physical space
      Jinv.leftMultiplyTranspose( dNdXi.data, dNdX[i] );
    }
  }

  GEOS_HOST_DEVICE
  static void getGradient( real64 const (&Xiq)[3],
                           typename CubeCell::JacobianType const & Jinv,
                           real64 (& dNdX)[8][3] )
  {
    for( int i = 0; i < 8; ++i )
    {
      // ... ... Parent space
      Gradient dNdXi = getParentGradient( i, Xiq );

      // ... ... Physical space
      Jinv.leftMultiplyTranspose( dNdXi.data, dNdX[i] );
    }
  }
};

template<>
struct Helper< CellUtilities::ReferenceCell::Wedge,//WedgeCell,
               BasisFunction::Lagrange,
               6 >
{
  GEOS_HOST_DEVICE
  static Gradient getParentGradient( int const BasisIndex,
                                     real64 const (&Xiq)[3] )
  {
    constexpr real64 dpsiTRI[2][3] = { { -1.0, 1.0, 0.0 }, { -1.0, 0.0, 1.0 } };
    constexpr real64 dpsiLIN[2] = { -1.0, 1.0 };

    int const a = BasisIndex / 2;
    int const b = BasisIndex & 1;


    Gradient gradient;
    gradient.data[0] = dpsiTRI[0][a] * 0.5 * ( 1.0 + dpsiLIN[b] * Xiq[2] );
    gradient.data[1] = dpsiTRI[1][a] * 0.5 * ( 1.0 + dpsiLIN[b] * Xiq[2] );
    gradient.data[2] = ( ( BasisIndex < 2 ) +  dpsiTRI[0][a] * Xiq[0] + dpsiTRI[1][a] * Xiq[1] )
                       * dpsiLIN[b];
    return gradient;
  }

  GEOS_HOST_DEVICE
  static void getGradient( real64 const (&Xiq)[3],
                           typename WedgeCell::JacobianType const & Jinv,
                           real64 (& dNdX)[6][3] )
  {
    for( int i = 0; i < 6; ++i )
    {
      // ... ... Parent space
      Gradient dNdXi = getParentGradient( i, Xiq );

      // ... ... Physical space
      Jinv.leftMultiplyTranspose( dNdXi.data, dNdX[i] );
    }
  }
};

template<>
struct Helper< CellUtilities::ReferenceCell::Tetrahedron,//TetrahedronCell,
               BasisFunction::Lagrange,
               4 >
{
  GEOS_HOST_DEVICE
  static Gradient getParentGradient( int const BasisIndex,
                                     real64 const (&Xiq)[3] )
  {
    GEOS_UNUSED_VAR( Xiq );
    Gradient gradient;
    gradient.data[0] = -1.0 + (BasisIndex > 0) + (BasisIndex == 1 );
    gradient.data[1] = -1.0 + (BasisIndex > 0) + (BasisIndex == 2 );
    gradient.data[2] = -1.0 + (BasisIndex > 0) + (BasisIndex == 3 );
    return gradient;
  }

  GEOS_HOST_DEVICE
  static void getGradient( real64 const (&Xiq)[3],
                           typename WedgeCell::JacobianType const & Jinv,
                           real64 (& dNdX)[4][3] )
  {
    GEOS_UNUSED_VAR( Xiq );
    for( int i = 0; i < 4; ++i )
    {
      // ... ... Parent space
      Gradient dNdXi = getParentGradient( i, Xiq );

      // ... ... Physical space
      Jinv.leftMultiplyTranspose( dNdXi.data, dNdX[i] );
    }
  }
};

template<>
struct Helper< CellUtilities::ReferenceCell::Pyramid,//PyramidCell,
               BasisFunction::Lagrange,
               5 >
{
  GEOS_HOST_DEVICE
  static Gradient getParentGradient( int const BasisIndex,
                                     real64 const (&Xiq)[3] )
  {
    constexpr real64 dPhiLin[2] = { -1.0, 1.0 };

    int const a = BasisIndex & 1;
    int const b = ( BasisIndex & 2 ) >> 1;
    int const c = ( BasisIndex & 4 ) >> 2;

    Gradient gradient;
    gradient.data[0] = (1 - c ) * 0.125 * (       dPhiLin[a]          ) * ( 1.0 + dPhiLin[b] * Xiq[1] ) * ( 1.0 + dPhiLin[c] * Xiq[2] );
    gradient.data[1] = (1 - c ) * 0.125 * ( 1.0 + dPhiLin[a] * Xiq[0] ) * (       dPhiLin[b]          ) * ( 1.0 + dPhiLin[c] * Xiq[2] );
    gradient.data[2] = (1 - c ) * 0.125 * ( 1.0 + dPhiLin[a] * Xiq[0] ) * ( 1.0 + dPhiLin[b] * Xiq[1] ) * (       dPhiLin[c]          ) + c * 0.5;
    return gradient;
  }

  GEOS_HOST_DEVICE
  static void getGradient( real64 const (&Xiq)[3],
                           typename PyramidCell::JacobianType const & Jinv,
                           real64 (& dNdX)[5][3] )
  {
    for( int i = 0; i < 5; ++i )
    {
      // ... ... Parent space
      Gradient dNdXi = getParentGradient( i, Xiq );

      // ... ... Physical space
      Jinv.leftMultiplyTranspose( dNdXi.data, dNdX[i] );
    }
  }
};

// template< typename CELL_TYPE,
//           BasisFunction BASIS_FUNCTION >
// GEOS_HOST_DEVICE
// static Gradient getParentGradient( int const BasisIndex,
//                                    real64 const (&Xiq)[3] )
// {
//   return Helper< CELL_TYPE, BASIS_FUNCTION >::getParentGradient( BasisIndex, Xiq );
// }

template< typename CELL_TYPE,
          BasisFunction BASIS_FUNCTION,
          int NUM_SUPPORT_POINTS >
GEOS_HOST_DEVICE
static void getGradient( real64 const (&Xiq)[3],
                         typename CELL_TYPE::JacobianType const & Jinv,
                         real64 (& dNdX)[NUM_SUPPORT_POINTS][3] )
{
  // for( int i = 0; i < NUM_SUPPORT_POINTS; ++i )
  // {
  //   // ... ... Parent space
  //   Gradient dNdXi = getParentGradient< CELL_TYPE, BASIS_FUNCTION >( i, Xiq );

  //   // ... ... Physical space
  //   Jinv.leftMultiplyTranspose( dNdXi.data, dNdX[i] );
  // }
  Helper< CellUtilities::ParentCell< CELL_TYPE >::value, BASIS_FUNCTION, NUM_SUPPORT_POINTS >::getGradient( Xiq, Jinv, dNdX );
}



/// Declare strings associated with enumeration values.
ENUM_STRINGS( BasisFunction,
              "Lagrange" );

} // namespace BasisFunctionUtilities

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_SIMPLEPDE_BASISFUNCTIONUTILITIES_HPP_
