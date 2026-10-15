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

#include "physicsSolvers/solidMechanics/kernels/ExplicitMPM.hpp"

#include <gtest/gtest.h>

namespace geos
{
namespace
{

struct MockCohesiveZoneUpdates
{
  explicit MockCohesiveZoneUpdates( arrayView1d< integer > const calls )
    : m_calls( calls )
  {}

  GEOS_HOST_DEVICE
  void jumpDisplacementUpdate( localIndex const pairIndex,
                               real64 const & normalDisplacement,
                               real64 const & tangentialDisplacement,
                               real64 & normalStress,
                               real64 & shearStress ) const
  {
    GEOS_UNUSED_VAR( normalDisplacement );
    GEOS_UNUSED_VAR( tangentialDisplacement );
    ++m_calls[pairIndex];
    normalStress = -2.0;
    shearStress = 0.0;
  }

  GEOS_HOST_DEVICE
  void saveConvergedState( localIndex const pairIndex,
                           localIndex const quadraturePoint ) const
  {
    GEOS_UNUSED_VAR( pairIndex );
    GEOS_UNUSED_VAR( quadraturePoint );
  }

  arrayView1d< integer > const m_calls;
};

TEST( CohesiveZoneStateUpdateKernel, AccumulatesThreeFieldPairsIntoSharedSlots )
{
  localIndex constexpr numPairs = 3;
  localIndex constexpr numFieldSlots = 3;

  array2d< localIndex > pairToFieldSlot( numPairs, 2 );
  pairToFieldSlot[0][0] = 0;
  pairToFieldSlot[0][1] = 1;
  pairToFieldSlot[1][0] = 0;
  pairToFieldSlot[1][1] = 2;
  pairToFieldSlot[2][0] = 1;
  pairToFieldSlot[2][1] = 2;

  array1d< real64 > fieldSlotMass( numFieldSlots );
  array2d< real64 > fieldSlotDisplacement( numFieldSlots, 3 );
  array2d< real64 > fieldSlotParticleSurfaceNormal( numFieldSlots, 3 );
  array3d< real64 > fieldSlotDeformationGradientCofactor( numFieldSlots, 3, 3 );
  array2d< real64 > fieldSlotCohesiveForce( numFieldSlots, 3 );
  array2d< real64 > fieldSlotReferenceSurfaceNormal( numFieldSlots, 3 );
  array1d< real64 > fieldSlotReferenceArea( numFieldSlots );
  array1d< integer > constitutiveCalls( numPairs );

  real64 const normals[numFieldSlots][3] = {
    { 1.0, 0.0, 0.0 },
    { -1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 }
  };

  for( localIndex slot = 0; slot < numFieldSlots; ++slot )
  {
    fieldSlotMass[slot] = 1.0;
    fieldSlotReferenceArea[slot] = 1.0;
    for( localIndex i = 0; i < 3; ++i )
    {
      fieldSlotDisplacement[slot][i] = 0.0;
      fieldSlotCohesiveForce[slot][i] = 0.0;
      fieldSlotParticleSurfaceNormal[slot][i] = normals[slot][i];
      fieldSlotReferenceSurfaceNormal[slot][i] = normals[slot][i];
      for( localIndex j = 0; j < 3; ++j )
      {
        fieldSlotDeformationGradientCofactor[slot][i][j] = i == j ? 1.0 : 0.0;
      }
    }
  }
  for( localIndex pair = 0; pair < numPairs; ++pair )
  {
    constitutiveCalls[pair] = 0;
  }

  MockCohesiveZoneUpdates const constitutiveWrapper( constitutiveCalls );
  solidMechanicsMPMKernels::CohesiveZoneStateUpdateKernel::launch< parallelDevicePolicy<> >(
    numPairs,
    constitutiveWrapper,
    1.0,
    0,
    1.0e-12,
    0,
    pairToFieldSlot,
    0,
    0,
    0,
    1.0,
    1.0,
    1.0,
    fieldSlotMass,
    fieldSlotDisplacement,
    fieldSlotParticleSurfaceNormal,
    fieldSlotDeformationGradientCofactor,
    fieldSlotCohesiveForce,
    fieldSlotReferenceSurfaceNormal,
    fieldSlotReferenceArea );
  parallelDeviceSync();
  fieldSlotCohesiveForce.move( hostMemorySpace, true );
  constitutiveCalls.move( hostMemorySpace, true );

  real64 netForce[3] = {};
  real64 absoluteForce = 0.0;
  for( localIndex slot = 0; slot < numFieldSlots; ++slot )
  {
    for( localIndex i = 0; i < 3; ++i )
    {
      netForce[i] += fieldSlotCohesiveForce[slot][i];
      absoluteForce += LvArray::math::abs( fieldSlotCohesiveForce[slot][i] );
    }
  }

  EXPECT_GT( absoluteForce, 0.0 );
  EXPECT_NEAR( netForce[0], 0.0, 1.0e-12 );
  EXPECT_NEAR( netForce[1], 0.0, 1.0e-12 );
  EXPECT_NEAR( netForce[2], 0.0, 1.0e-12 );
  for( localIndex pair = 0; pair < numPairs; ++pair )
  {
    EXPECT_EQ( constitutiveCalls[pair], 1 );
  }

  // Capturing diagnostics must preserve one constitutive update per pair and the
  // original slot forces, even when multiple pairs share a slot.
  using D = solidMechanicsMPMKernels::CohesivePairDebug;
  array2d< real64 > dumpedForce( numFieldSlots, 3 );
  array1d< integer > dumpedCalls( numPairs );
  array2d< real64 > debugValues( numPairs, D::numColumns );
  array1d< int > debugStatus( numPairs );
  for( localIndex slot = 0; slot < numFieldSlots; ++slot )
  {
    for( int i = 0; i < 3; ++i )
    {
      dumpedForce[slot][i] = 0.0;
    }
  }
  for( localIndex pair = 0; pair < numPairs; ++pair )
  {
    dumpedCalls[pair] = 0;
  }
  MockCohesiveZoneUpdates const dumpedWrapper( dumpedCalls );
  solidMechanicsMPMKernels::CohesiveZoneStateUpdateKernel::launch<
    parallelDevicePolicy<>, MockCohesiveZoneUpdates, true >(
    numPairs, dumpedWrapper, 1.0, 0, 1.0e-12, 0, pairToFieldSlot,
    0, 0, 0, 1.0, 1.0, 1.0, fieldSlotMass, fieldSlotDisplacement,
    fieldSlotParticleSurfaceNormal, fieldSlotDeformationGradientCofactor,
    dumpedForce, fieldSlotReferenceSurfaceNormal, fieldSlotReferenceArea,
    debugValues, debugStatus );
  parallelDeviceSync();
  dumpedForce.move( hostMemorySpace, true );
  dumpedCalls.move( hostMemorySpace, true );
  debugValues.move( hostMemorySpace, true );
  debugStatus.move( hostMemorySpace, true );
  pairToFieldSlot.move( hostMemorySpace, false );

  real64 contributions[numFieldSlots][3] = {};
  for( localIndex pair = 0; pair < numPairs; ++pair )
  {
    EXPECT_EQ( dumpedCalls[pair], 1 );
    EXPECT_EQ( debugStatus[pair], D::evaluated );
    EXPECT_DOUBLE_EQ( debugValues[pair][D::returnedNormalStress], -2.0 );
    EXPECT_DOUBLE_EQ( debugValues[pair][D::appliedNormalStress], -2.0 );
    for( int i = 0; i < 3; ++i )
    {
      real64 const forceA = debugValues[pair][D::pairForceA + i];
      real64 const forceB = debugValues[pair][D::pairForceB + i];
      EXPECT_NEAR( forceA + forceB, 0.0, 1.0e-12 );
      EXPECT_DOUBLE_EQ( forceA, debugValues[pair][D::tractionA + i] * debugValues[pair][D::pairArea] );
      contributions[pairToFieldSlot[pair][0]][i] += forceA;
      contributions[pairToFieldSlot[pair][1]][i] += forceB;
    }
  }
  // Opposed normals with unit reference areas give area 1, not the sum of both sides.
  EXPECT_DOUBLE_EQ( debugValues[0][D::pairArea], 1.0 );
  for( localIndex slot = 0; slot < numFieldSlots; ++slot )
  {
    for( int i = 0; i < 3; ++i )
    {
      EXPECT_NEAR( dumpedForce[slot][i], fieldSlotCohesiveForce[slot][i], 1.0e-12 );
      EXPECT_NEAR( contributions[slot][i], dumpedForce[slot][i], 1.0e-12 );
    }
  }
}

TEST( CohesiveZoneStateUpdateKernel, CapturesEachSkipReasonWithoutAddingForces )
{
  using D = solidMechanicsMPMKernels::CohesivePairDebug;
  int constexpr numPairs = 5;
  int constexpr numSlots = 2 * numPairs;
  array2d< localIndex > pairs( numPairs, 2 );
  array1d< real64 > mass( numSlots );
  array2d< real64 > displacement( numSlots, 3 );
  array2d< real64 > normal( numSlots, 3 );
  array3d< real64 > cofactor( numSlots, 3, 3 );
  array2d< real64 > force( numSlots, 3 );
  array2d< real64 > referenceNormal( numSlots, 3 );
  array1d< real64 > area( numSlots );
  array1d< integer > calls( numPairs );
  array2d< real64 > values( numPairs, D::numColumns );
  array1d< int > status( numPairs );
  for( int pair = 0; pair < numPairs; ++pair )
  {
    pairs[pair][0] = 2 * pair;
    pairs[pair][1] = 2 * pair + 1;
    calls[pair] = 0;
  }
  for( int slot = 0; slot < numSlots; ++slot )
  {
    mass[slot] = 1.0;
    area[slot] = 1.0;
    for( int i = 0; i < 3; ++i )
    {
      displacement[slot][i] = force[slot][i] = 0.0;
      normal[slot][i] = referenceNormal[slot][i] = i == 0 ? ( slot % 2 == 0 ? 1.0 : -1.0 ) : 0.0;
      for( int j = 0; j < 3; ++j )
      {
        cofactor[slot][i][j] = i == j ? 1.0 : 0.0;
      }
    }
  }
  mass[0] = 0.0;
  normal[2][0] = 0.0;
  mass[5] = 0.0;
  normal[7][0] = 0.0;
  normal[9][0] = 1.0; // Equal normals cancel in the mass-weighted pair normal.

  MockCohesiveZoneUpdates const wrapper( calls );
  solidMechanicsMPMKernels::CohesiveZoneStateUpdateKernel::launch<
    parallelDevicePolicy<>, MockCohesiveZoneUpdates, true >(
    numPairs, wrapper, 1.0, 0, 1.0e-12, 0, pairs,
    0, 0, 0, 1.0, 1.0, 1.0, mass, displacement, normal, cofactor,
    force, referenceNormal, area, values, status );
  parallelDeviceSync();
  force.move( hostMemorySpace, true );
  calls.move( hostMemorySpace, true );
  status.move( hostMemorySpace, true );
  values.move( hostMemorySpace, true );

  int const expected[numPairs] = { D::inactive, D::inactiveNormalA, D::inactiveMassB,
                                   D::inactiveNormalB, D::degenerateNormal };
  for( int pair = 0; pair < numPairs; ++pair )
  {
    EXPECT_EQ( status[pair], expected[pair] );
    EXPECT_EQ( calls[pair], 0 );
    for( int i = 0; i < 3; ++i )
    {
      EXPECT_DOUBLE_EQ( values[pair][D::pairForceA + i], 0.0 );
      EXPECT_DOUBLE_EQ( values[pair][D::pairForceB + i], 0.0 );
    }
  }
  for( int slot = 0; slot < numSlots; ++slot )
  {
    for( int i = 0; i < 3; ++i )
    {
      EXPECT_DOUBLE_EQ( force[slot][i], 0.0 );
    }
  }
}

} // namespace
} // namespace geos

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
