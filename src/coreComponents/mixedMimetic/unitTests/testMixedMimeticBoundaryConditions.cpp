/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "mixedMimetic/MixedMimeticBoundaryConditions.hpp"
#include "mainInterface/initialization.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geos;
using namespace geos::mixedMimeticBoundary;

namespace
{

real64 constexpr area = 2.5;
real64 constexpr tol = 1.0e-13;

// outward flux of a cell of value xK through a non-condensed boundary face: r F - xK + trace( F ) = 0
real64 liveRowFlux( FaceBoundaryCondition const & bc, real64 const halfResistance, real64 const xK )
{
  real64 dTrace = 0.0;
  real64 const trace0 = bc.trace( 1.0, area, 0.0, dTrace );
  return ( xK - trace0 ) / ( halfResistance + dTrace );
}

// the same flux from the condensed two-point closure: ( r + resistance ) F = xK - g
real64 condensedFlux( FaceBoundaryCondition const & bc, real64 const halfResistance, real64 const xK )
{
  return ( xK - bc.value ) / ( halfResistance + bc.resistance( area ) );
}

}

TEST( MixedMimeticBoundaryConditions, interiorByDefault )
{
  FaceBoundaryCondition const bc;
  EXPECT_EQ( bc.type, BoundaryType::interior );
  EXPECT_FALSE( bc.isEssential() );
  EXPECT_FALSE( bc.isCondensable() );
}

TEST( MixedMimeticBoundaryConditions, dirichlet )
{
  FaceBoundaryCondition const bc{ BoundaryType::dirichlet, 3.0e5, 0.0 };
  EXPECT_FALSE( bc.isEssential() );
  EXPECT_TRUE( bc.isCondensable() );

  // x_f = g independently of the flux
  real64 dTrace = 1.0;
  EXPECT_DOUBLE_EQ( bc.trace( -1.0, area, 7.0, dTrace ), 3.0e5 );
  EXPECT_DOUBLE_EQ( dTrace, 0.0 );
  EXPECT_DOUBLE_EQ( bc.resistance( area ), 0.0 );
}

TEST( MixedMimeticBoundaryConditions, neumann )
{
  FaceBoundaryCondition const bc{ BoundaryType::neumann, -4.0, 0.0 };
  EXPECT_TRUE( bc.isEssential() );
  EXPECT_FALSE( bc.isCondensable() );

  // the face dof is sigma g |f|, so that the outward flux sigma f is g |f| for both orientations
  for( real64 const sigma : { 1.0, -1.0 } )
  {
    EXPECT_DOUBLE_EQ( sigma * bc.essentialFlux( sigma, area ), -4.0 * area );
  }

  // no-flow
  FaceBoundaryCondition const noFlow{ BoundaryType::neumann, 0.0, 0.0 };
  EXPECT_DOUBLE_EQ( noFlow.essentialFlux( -1.0, area ), 0.0 );
}

TEST( MixedMimeticBoundaryConditions, robin )
{
  real64 const alpha = 25.0;
  real64 const g = 290.0;
  FaceBoundaryCondition const bc{ BoundaryType::robin, g, alpha };
  EXPECT_FALSE( bc.isEssential() );
  EXPECT_TRUE( bc.isCondensable() );

  // the trace satisfies the law sigma f = alpha |f| ( x_f - g )
  real64 const outwardFlux = 12.0;
  real64 dTrace = 0.0;
  real64 const trace = bc.trace( 1.0, area, outwardFlux, dTrace );
  EXPECT_NEAR( alpha * area * ( trace - g ), outwardFlux, tol * outwardFlux );
  EXPECT_DOUBLE_EQ( dTrace, 1.0 / ( alpha * area ) );
  EXPECT_DOUBLE_EQ( bc.resistance( area ), 1.0 / ( alpha * area ) );
}

TEST( MixedMimeticBoundaryConditions, liveAndCondensedClosuresAgree )
{
  real64 const halfResistance = 0.4;
  real64 const xK = 350.0;
  FaceBoundaryCondition const conditions[2] = { { BoundaryType::dirichlet, 300.0, 0.0 },
    { BoundaryType::robin, 300.0, 5.0 } };
  for( FaceBoundaryCondition const & bc : conditions )
  {
    real64 const flux = liveRowFlux( bc, halfResistance, xK );
    EXPECT_NEAR( flux, condensedFlux( bc, halfResistance, xK ), tol * flux );
  }
}

TEST( MixedMimeticBoundaryConditions, robinLimits )
{
  real64 const halfResistance = 0.4;
  real64 const xK = 350.0;
  real64 const g = 300.0;
  FaceBoundaryCondition const dirichlet{ BoundaryType::dirichlet, g, 0.0 };
  real64 const dirichletFlux = liveRowFlux( dirichlet, halfResistance, xK );

  // alpha -> infinity: Dirichlet
  FaceBoundaryCondition const stiff{ BoundaryType::robin, g, 1.0e14 };
  EXPECT_NEAR( liveRowFlux( stiff, halfResistance, xK ), dirichletFlux, 1.0e-12 * dirichletFlux );

  // alpha -> 0: the flux vanishes as alpha |f| ( x_K - g )
  real64 const alpha = 1.0e-10;
  FaceBoundaryCondition const soft{ BoundaryType::robin, g, alpha };
  EXPECT_NEAR( liveRowFlux( soft, halfResistance, xK ), alpha * area * ( xK - g ), 1.0e-9 * alpha * area * ( xK - g ) );
}

TEST( MixedMimeticBoundaryConditions, faceView )
{
  array1d< integer > type( 3 );
  array1d< real64 > value( 3 );
  array1d< real64 > coefficient( 3 );
  type[0] = BoundaryType::interior;
  type[1] = BoundaryType::neumann; value[1] = 2.0;
  type[2] = BoundaryType::robin; value[2] = 290.0; coefficient[2] = 25.0;

  FaceBoundaryView const view{ type.toViewConst(), value.toViewConst(), coefficient.toViewConst() };
  EXPECT_EQ( view[0].type, BoundaryType::interior );
  EXPECT_TRUE( view[1].isEssential() );
  EXPECT_DOUBLE_EQ( view[1].value, 2.0 );
  EXPECT_DOUBLE_EQ( view[2].resistance( area ), 1.0 / ( 25.0 * area ) );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
