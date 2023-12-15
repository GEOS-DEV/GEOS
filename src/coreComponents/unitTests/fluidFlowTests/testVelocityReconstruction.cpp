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

// Source includes
#include "common/Logger.hpp"
#include "mainInterface/initialization.hpp"
#include "finiteVolume/CellElementStencilTPFA.hpp"
#include "testFlowKernelHelpers.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geos;

void setReferences( const real64 flux, const real64 dx, const real64 dy, real64 (& expectedVel)[3] )
{

  expectedVel[0] = flux/dx;
  expectedVel[1] = flux/dy;
  expectedVel[2] = 0.0;

}

TEST( testAligned2D, Velocity_aligned2D )
{
  constexpr int nfaces = 4;    //2-hexa problem
  constexpr int numElemtInStencil = 2;

  real64 const faceNormal[nfaces][3] = {{1., 0., 0.},
    {1., 0., 0.},
    {0., 1., 0.},
    {0., 1., 0.}};
  real64 const cellToFaceVec[nfaces][numElemtInStencil][3] = { {{1., 0., 0.}, {-1., 0., 0}},
    {{1., 0., 0.}, {-1., 0., 0}},
    {{0., 1., 0.}, {0., -1., 0}},
    {{0., 1., 0.}, {0., -1., 0}}
  };
  real64 const transMult[nfaces] = {1., 1., 1., 1.};
  real64 weight[] = {1, -1};
  real64 const geomStabSum = 1.;    //irrelevant for now

  localIndex elementRegionIndices[] = {0, 0};
  localIndex elementSubRegionIndices[] = {0, 0};
  localIndex ei_[4][numElemtInStencil] = {
    {0, 1},       //face #1
    {2, 3},       //face #2
    {0, 2},       //face #3
    {1, 3}       //face #4
  };

  const real64 flux = 0.001;
  const int numElem = 4;
  constexpr localIndex nPhases = 2;
  constexpr localIndex nDirs = 3;
  const real64 phaseVelocity[numElem][nPhases][nDirs] = {
    {{0., 0., 0.}, {0., 0., 0.}},
    {{0., 0., 0.}, {0., 0., 0.}},
    {{0., 0., 0.}, {0., 0., 0.}},
    {{0., 0., 0.}, {0., 0., 0.}}
  };

  CellElementStencilTPFA tpfa;
  for( int kf = 0; kf < nfaces; ++kf )
  {

    localIndex elementIndices[] = {ei_[kf][0], ei_[kf][1]};
    tpfa.add( 2, elementRegionIndices, elementSubRegionIndices, elementIndices, weight, kf );
    tpfa.addVectors( transMult[kf], geomStabSum, faceNormal[kf], cellToFaceVec[kf] );
  }

  CellElementStencilTPFA::KernelWrapper wrapper = tpfa.createKernelWrapper();
  array2d< real64 > globalCellDim( 2, 3 );
  globalCellDim[0][0] = 2.; globalCellDim[1][0] = 2.;
  globalCellDim[0][1] = 2.; globalCellDim[1][1] = 2.;
  globalCellDim[0][2] = 0.; globalCellDim[1][2] = 0.;
  arrayView2d< const real64 > const globalCellDimView = globalCellDim.toViewConst();

  CellElementStencilTPFA::IndexContainerViewConstType const & seri = tpfa.getElementRegionIndices();
  CellElementStencilTPFA::IndexContainerViewConstType const & sesri = tpfa.getElementSubRegionIndices();
  CellElementStencilTPFA::IndexContainerViewConstType const & sei = tpfa.getElementIndices();

  ElementRegionManager::ElementViewAccessor< array3d< real64 > > phaseVelocityView = AccessorHelper< true >::makeElementAccessor< 3 >( phaseVelocity[0][0],
                                                                                                                                       tpfa.stencilSize( 0 ),
                                                                                                                                       seri[nfaces-1],
                                                                                                                                       sesri[nfaces-1],
                                                                                                                                       sei[nfaces-1],
                                                                                                                                       nPhases,
                                                                                                                                       nDirs );
  for( int iconn = 0; iconn < nfaces; ++iconn )
  {

    wrapper.computeVelocity( iconn /*iconn*/, 0 /*ip*/, flux, {globalCellDimView[0], globalCellDimView[1]}, phaseVelocityView.toNestedView());
    wrapper.computeVelocity( iconn /*iconn*/, 1 /*ip*/, 100*flux, {globalCellDimView[0], globalCellDimView[1]}, phaseVelocityView.toNestedView());
  }

  for( int ip = 0; ip < 2; ++ip )
  {
    auto fip = (99*ip+1)*flux;    //shorthand to get the 1:100 ratio between phases as imposed above
    for( int ib = 0; ib < numElem; ++ib )
    {
      real64 expectedVel[3];
      setReferences( fip, globalCellDim[0][0], globalCellDim[0][1], expectedVel );
      EXPECT_EQ( phaseVelocityView[0][0][ib][ip][0], expectedVel[0] );
      EXPECT_EQ( phaseVelocityView[0][0][ib][ip][1], expectedVel[1] );
      EXPECT_EQ( phaseVelocityView[0][0][ib][ip][2], expectedVel[2] );
    }
  }
}

int main( int argc, char * *argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
