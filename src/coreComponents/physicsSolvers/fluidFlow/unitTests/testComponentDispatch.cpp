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

#include "physicsSolvers/fluidFlow/kernels/compositional/CompositionalMultiphaseHybridFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/ReactiveCompositionalMultiphaseOBLKernels.hpp"
#include "mainInterface/initialization.hpp"

#include <gtest/gtest.h>

using namespace geos;

namespace
{

TEST( ComponentDispatch, FluidStorageCapacity )
{
  EXPECT_EQ( constitutive::MultiFluidConstants::MAX_NUM_COMPONENTS, std::max( 9, GEOS_MAX_COMPONENTS ) );
}

struct HybridDispatchProbe
{
  template< integer NF, integer NC, integer NP, typename IP >
  static void launch( integer & faces, integer & components, integer & phases )
  {
    static_assert( NC <= GEOS_MAX_COMPONENTS, "Disabled component counts must not be instantiated" );
    faces = NF;
    components = NC;
    phases = NP;
  }
};

TEST( ComponentDispatch, ConfiguredRange )
{
  for( integer components = 1; components <= GEOS_MAX_COMPONENTS; ++components )
  {
    integer calls = 0;
    geos::internal::kernelLaunchSelectorCompSwitch( components, [&] ( auto NC )
    {
      static_assert( NC() <= GEOS_MAX_COMPONENTS, "Disabled component counts must not be instantiated" );
      EXPECT_EQ( NC(), components );
      ++calls;
    } );
    EXPECT_EQ( calls, 1 );

    for( bool thermal : { false, true } )
    {
      geos::internal::kernelLaunchSelectorCompThermSwitch( components, thermal, [&] ( auto NC, auto THERMAL )
      {
        EXPECT_EQ( NC(), components );
        EXPECT_EQ( THERMAL(), thermal );
      } );
    }

    for( integer phases = 1; phases <= 3; ++phases )
    {
      geos::internal::kernelLaunchSelectorCompPhaseSwitch( components, phases, [&] ( auto NC, auto NP )
      {
        EXPECT_EQ( NC(), components );
        EXPECT_EQ( NP(), phases );
      } );
    }
  }
}

TEST( ComponentDispatch, OBLRange )
{
  for( integer components = 1; components <= GEOS_MAX_OBL_COMPONENTS; ++components )
  {
    for( integer phases = 1; phases <= 3; ++phases )
    {
      for( bool energy : { false, true } )
      {
        reactiveCompositionalMultiphaseOBLKernels::internal::kernelLaunchSelectorEnergySwitch(
          phases, components, energy, [&] ( auto NP, auto NC, auto ENERGY )
        {
          static_assert( NC() <= GEOS_MAX_OBL_COMPONENTS, "Disabled component counts must not be instantiated" );
          EXPECT_EQ( NC(), components );
          EXPECT_EQ( NP(), phases );
          EXPECT_EQ( ENERGY(), energy );
        } );
      }
    }
  }
}

TEST( ComponentDispatch, HybridRange )
{
  for( integer components = 1; components <= GEOS_MAX_COMPONENTS; ++components )
  {
    for( integer phases : { 2, 3 } )
    {
      for( integer faces = 4; faces <= 13; ++faces )
      {
        integer actualFaces = 0;
        integer actualComponents = 0;
        integer actualPhases = 0;
        compositionalMultiphaseHybridFVMKernels::kernelLaunchSelector< HybridDispatchProbe, void >(
          faces, components, phases, actualFaces, actualComponents, actualPhases );
        EXPECT_EQ( actualFaces, faces );
        EXPECT_EQ( actualComponents, components );
        EXPECT_EQ( actualPhases, phases );
      }
    }
  }
}

} // namespace

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
