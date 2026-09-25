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
}

} // namespace
} // namespace geos

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
