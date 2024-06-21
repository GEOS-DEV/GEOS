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

#ifndef GEOS_PHYSICSSOLVERS_SIMPLEPDE_QUADRATUREUTILITIES_HPP_
#define GEOS_PHYSICSSOLVERS_SIMPLEPDE_QUADRATUREUTILITIES_HPP_

#include "CellUtilities.hpp"

namespace geos
{

namespace QuadratureUtilities
{

/**
 * @brief Integration rule type
 */
enum class Rule : integer
{
  Gauss
};

struct Data
{
  real64 wq;
  real64 Xiq[3];
};

template< CellUtilities::ReferenceCell REFERENCE_CELL_TYPE,
          Rule RULE,
          int N_INTEGRATION_POINTS >
struct Helper
{};

template<>
struct Helper< CellUtilities::ReferenceCell::Cube,
               Rule::Gauss,
               8 >
{
  GEOS_HOST_DEVICE
  static Data getData( int q )
  {
    real64 const val = 0.5773502691896257645092;
    int const a = q & 1;
    int const b = ( q & 2 ) >> 1;
    int const c = ( q & 4 ) >> 2;

    Data data;
    data.wq = 1.0; // weight
    data.Xiq[0] = ( 2 * a - 1 ) * val;
    data.Xiq[1] = ( 2 * b - 1 ) * val;
    data.Xiq[2] = ( 2 * c - 1 ) * val;
    return data;
  }
};

template<>
struct Helper< CellUtilities::ReferenceCell::Wedge,
               Rule::Gauss,
               6 >
{
  GEOS_HOST_DEVICE
  inline
  constexpr static real64 parentCoords0( int const a )
  {
    return 0.5 * ( a & 2 );
  }

  GEOS_HOST_DEVICE
  inline
  constexpr static real64 parentCoords1( int const a )
  {
    return 0.25 * ( a & 4 );
  }

  GEOS_HOST_DEVICE
  inline
  constexpr static real64 parentCoords2( int const a )
  {
    return -1.0 + 2 * ( a & 1 );
  }

  GEOS_HOST_DEVICE
  static Data getData( int q )
  {
    constexpr static real64 quadratureCrossSectionCoord = 1.0 / 6.0;
    constexpr static real64 quadratureLongitudinalCoord = 1.0 / 1.732050807568877293528;

    Data data;
    data.wq = 1.0 / 6.0; // weight
    data.Xiq[0] = quadratureCrossSectionCoord + 0.5 * parentCoords0( q );
    data.Xiq[1] = quadratureCrossSectionCoord + 0.5 * parentCoords1( q );
    data.Xiq[2] = quadratureLongitudinalCoord * parentCoords2( q );
    return data;
  }
};

template<>
struct Helper< CellUtilities::ReferenceCell::Tetrahedron,
               Rule::Gauss,
               1 >
{
  GEOS_HOST_DEVICE
  static Data getData( int q )
  {
    GEOS_UNUSED_VAR( q );

    Data data;
    data.wq = 1.0 / 1.6;
    data.Xiq[0] = 0.25;
    data.Xiq[1] = 0.25;
    data.Xiq[2] = 0.25;
    return data;
  }
};

template<>
struct Helper< CellUtilities::ReferenceCell::Pyramid,
               Rule::Gauss,
               5 >
{
  constexpr static real64 quadratureCrossSectionCoord = 0.584237394672177;
  constexpr static real64 quadratureLongitudinalCoordNeg = -2.0 / 3.0;
  constexpr static real64 quadratureLongitudinalCoordDelta = 16.0 / 15.0;

  constexpr static real64 weight = 81.0 / 100.0;
  constexpr static real64 weightDelta  = 125.0 / 27.0 - weight;

  GEOS_HOST_DEVICE
  inline
  constexpr static real64 parentCoords0( int const a )
  {
    return -1.0 + 2.0 * ( a & 1 ) + 0.25 * ( a & 4 );
  }

  GEOS_HOST_DEVICE
  inline
  constexpr static real64 parentCoords1( int const a )
  {
    return -1.0 + ( a & 2 ) + 0.25 * ( a & 4 );
  }

  GEOS_HOST_DEVICE
  inline
  constexpr static real64 parentCoords2( int const a )
  {
    return -1.0 + 0.5 * ( a & 4 );
  }

  GEOS_HOST_DEVICE
  inline
  constexpr static real64 quadratureParentCoords0( localIndex const q )
  {
    return parentCoords0( q ) * quadratureCrossSectionCoord;
  }

  GEOS_HOST_DEVICE
  inline
  constexpr static real64 quadratureParentCoords1( localIndex const q )
  {
    return parentCoords1( q ) * quadratureCrossSectionCoord;
  }

  GEOS_HOST_DEVICE
  inline
  constexpr static real64 quadratureParentCoords2( localIndex const q )
  {
    return quadratureLongitudinalCoordNeg + 0.5 * ( 1 + parentCoords2( q ) ) * quadratureLongitudinalCoordDelta;
  }

  GEOS_HOST_DEVICE
  static Data getData( int q )
  {
    Data data;
    data.wq = weight + 0.5 * ( 1 + parentCoords2( q ) ) * weightDelta;
    data.Xiq[0] = quadratureParentCoords0( q );
    data.Xiq[1] = quadratureParentCoords1( q );
    data.Xiq[2] = quadratureParentCoords2( q );
    return data;
  }
};

// getQuadratureData< CELL_TYPE, INTEGRATION_RULE, INTEGRATION_ORDER >

template< typename CELL_TYPE,
          Rule RULE_TYPE,
          int N_INTEGRATION_POINTS >
GEOS_HOST_DEVICE
static Data getData( int quadraturePointIndex )
{
  return Helper< CellUtilities::ParentCell< CELL_TYPE >::value, RULE_TYPE, N_INTEGRATION_POINTS >::getData( quadraturePointIndex );
}

/// Declare strings associated with enumeration values.
ENUM_STRINGS( Rule,
              "Gauss" );

} // namespace QuadratureUtilities

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_SIMPLEPDE_QUADRATUREUTILITIES_HPP_
