/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2020-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// using some utility classes from the following unit test
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"

#include "common/DataTypes.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/wavePropagation/WaveSolverBase.hpp"
#include "physicsSolvers/wavePropagation/WaveSolverBaseFields.hpp"
#include "physicsSolvers/wavePropagation/AcousticFirstOrderWaveEquationSEM.hpp"

#include <gtest/gtest.h>

using namespace geosx;
using namespace geosx::dataRepository;
using namespace geosx::testing;

CommandLineOptions g_commandLineOptions;

// This unit test checks the interpolation done to extract seismic traces from a wavefield.
// It computes a seismogram at a receiver co-located with the source and compares it to the surrounding receivers.
char const * xmlInput =
  "<?xml version=\"1.0\" ?>\n"
  "<Problem>\n"
  "  <Solvers>\n"
  "    <AcousticFirstOrderSEM\n"
  "      name=\"acousticFirstOrderSolver\"\n"
  "      cflFactor=\"0.25\"\n"
  "      discretization=\"FE1\"\n"
  "      targetRegions=\"{ Region }\"\n"
  "      sourceCoordinates=\"{ { 30, 30, 30 } }\"\n"
  "      timeSourceFrequency=\"2\"\n"
  "      receiverCoordinates=\"{ { 0.1, 0.1, 0.1 }, { 0.1, 0.1, 99.9 }, { 0.1, 99.9, 0.1 }, { 0.1, 99.9, 99.9 },\n"
  "                              { 99.9, 0.1, 0.1 }, { 99.9, 0.1, 99.9 }, { 99.9, 99.9, 0.1 }, { 99.9, 99.9, 99.9 },\n"
  "                              { 50, 50, 50 } }\"\n"
  "      outputSeismoTrace=\"0\"\n"
  "      dtSeismoTrace=\"0.05\"\n"
  "      rickerOrder=\"1\"/>\n"
  "  </Solvers>\n"
  "  <Mesh>\n"
  "    <InternalMesh\n"
  "      name=\"mesh\"\n"
  "      elementTypes=\"{ C3D8 }\"\n"
  "      xCoords=\"{ 0, 100 }\"\n"
  "      yCoords=\"{ 0, 100 }\"\n"
  "      zCoords=\"{ 0, 100 }\"\n"
  "      nx=\"{ 1 }\"\n"
  "      ny=\"{ 1 }\"\n"
  "      nz=\"{ 1 }\"\n"
  "      cellBlockNames=\"{ cb }\"/>\n"
  "  </Mesh>\n"
  "  <Events\n"
  "    maxTime=\"1\">\n"
  "    <PeriodicEvent\n"
  "      name=\"solverApplications\"\n"
  "      forceDt=\"0.05\"\n"
  "      targetExactStartStop=\"0\"\n"
  "      targetExactTimestep=\"0\"\n"
  "      target=\"/Solvers/acousticFirstOrderSolver\"/>\n"
  "    <PeriodicEvent\n"
  "      name=\"waveFieldUxCollection\"\n"
  "      timeFrequency=\"0.05\"\n"
  "      targetExactTimestep=\"0\"\n"
  "      target=\"/Tasks/waveFieldUxCollection\" />\n"
  "    <PeriodicEvent\n"
  "      name=\"waveFieldUyCollection\"\n"
  "      timeFrequency=\"0.05\"\n"
  "      targetExactTimestep=\"0\"\n"
  "      target=\"/Tasks/waveFieldUyCollection\" />\n"
  "    <PeriodicEvent\n"
  "      name=\"waveFieldUzCollection\"\n"
  "      timeFrequency=\"0.05\"\n"
  "      targetExactTimestep=\"0\"\n"
  "      target=\"/Tasks/waveFieldUzCollection\" />\n"
  "    <PeriodicEvent\n"
  "      name=\"waveFieldPressureCollection\"\n"
  "      timeFrequency=\"0.05\"\n"
  "      targetExactTimestep=\"0\"\n"
  "      target=\"/Tasks/waveFieldPressureCollection\" />\n"
  "  </Events>\n"
  "  <NumericalMethods>\n"
  "    <FiniteElements>\n"
  "      <FiniteElementSpace\n"
  "        name=\"FE1\"\n"
  "        order=\"1\"\n"
  "        formulation=\"SEM\" />\n"
  "    </FiniteElements>\n"
  "  </NumericalMethods>\n"
  "  <ElementRegions>\n"
  "    <CellElementRegion\n"
  "      name=\"Region\"\n"
  "      cellBlocks=\"{ cb }\"\n"
  "      materialList=\"{ nullModel }\"/>\n"
  "  </ElementRegions>\n"
  "  <Constitutive>\n"
  "    <NullModel\n"
  "      name=\"nullModel\"/>\n"
  "  </Constitutive>\n"
  "  <FieldSpecifications>\n"
  "    <FieldSpecification\n"
  "      name=\"cellVelocity\"\n"
  "      initialCondition=\"1\"\n"
  "      objectPath=\"ElementRegions/Region/cb\"\n"
  "      fieldName=\"mediumVelocity\"\n"
  "      scale=\"1500\"\n"
  "      setNames=\"{ all }\"/>\n"
  "    <FieldSpecification\n"
  "      name=\"cellDensity\"\n"
  "      initialCondition=\"1\"\n"
  "      objectPath=\"ElementRegions/Region/cb\"\n"
  "      fieldName=\"mediumDensity\"\n"
  "      scale=\"1\"\n"
  "      setNames=\"{ all }\"/>\n"
  "  </FieldSpecifications>\n"
  "  <Tasks>\n"
  "    <PackCollection\n"
  "      name=\"waveFieldPressureCollection\"\n"
  "      objectPath=\"nodeManager\"\n"
  "      fieldName=\"pressure_np1\"/>\n"
  "    <PackCollection\n"
  "      name=\"waveFieldUxCollection\"\n"
  "      objectPath=\"mesh/FE1/ElementRegions/Region/cb\"\n"
  "      fieldName=\"velocity_x\"/>\n"
  "    <PackCollection\n"
  "      name=\"waveFieldUyCollection\"\n"
  "      objectPath=\"mesh/FE1/ElementRegions/Region/cb\"\n"
  "      fieldName=\"velocity_y\"/>\n"
  "    <PackCollection\n"
  "      name=\"waveFieldUzCollection\"\n"
  "      objectPath=\"mesh/FE1/ElementRegions/Region/cb\"\n"
  "      fieldName=\"velocity_z\"/>\n"

  "  </Tasks>\n"
  "</Problem>\n";

class AcousticFirstOrderWaveEquationSEMTest : public ::testing::Test
{
public:

  AcousticFirstOrderWaveEquationSEMTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {
    setupProblemFromXML( state.getProblemManager(), xmlInput );
  }

  static real64 constexpr time = 0.0;
  static real64 constexpr dt = 5e-2;
  static real64 constexpr eps = std::numeric_limits< real64 >::epsilon();

  GeosxState state;
  AcousticFirstOrderWaveEquationSEM * propagator;
};

real64 constexpr AcousticFirstOrderWaveEquationSEMTest::time;
real64 constexpr AcousticFirstOrderWaveEquationSEMTest::dt;
real64 constexpr AcousticFirstOrderWaveEquationSEMTest::eps;

TEST_F( AcousticFirstOrderWaveEquationSEMTest, SeismoTrace )
{

  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  propagator = &state.getProblemManager().getPhysicsSolverManager().getGroup< AcousticFirstOrderWaveEquationSEM >( "acousticFirstOrderSolver" );
  real64 time_n = time;
  // run for 1s (20 steps)
  for( int i=0; i<20; i++ )
  {
    propagator->solverStep( time_n, dt, i, domain );
    time_n += dt;
  }
  // cleanup (triggers calculation of the remaining seismograms data points)
  propagator->cleanup( 1.0, 20, 0, 0, domain );

  // retrieve seismo
  arrayView2d< real32 > const pReceivers = propagator->getReference< array2d< real32 > >( AcousticFirstOrderWaveEquationSEM::viewKeyStruct::pressureNp1AtReceiversString() ).toView();
  arrayView2d< real32 > const uxReceivers = propagator->getReference< array2d< real32 > >( AcousticFirstOrderWaveEquationSEM::viewKeyStruct::uxNp1AtReceiversString() ).toView();
  arrayView2d< real32 > const uyReceivers = propagator->getReference< array2d< real32 > >( AcousticFirstOrderWaveEquationSEM::viewKeyStruct::uyNp1AtReceiversString() ).toView();
  arrayView2d< real32 > const uzReceivers = propagator->getReference< array2d< real32 > >( AcousticFirstOrderWaveEquationSEM::viewKeyStruct::uzNp1AtReceiversString() ).toView();

  // move it to CPU, if needed
  pReceivers.move( LvArray::MemorySpace::host, false );
  uxReceivers.move( LvArray::MemorySpace::host, false );
  uyReceivers.move( LvArray::MemorySpace::host, false );
  uzReceivers.move( LvArray::MemorySpace::host, false );



  // check number of seismos and trace length
  ASSERT_EQ( pReceivers.size( 1 ), 9 );
  ASSERT_EQ( pReceivers.size( 0 ), 21 );
  ASSERT_EQ( uxReceivers.size( 1 ), 9 );
  ASSERT_EQ( uxReceivers.size( 0 ), 21 );
  ASSERT_EQ( uyReceivers.size( 1 ), 9 );
  ASSERT_EQ( uyReceivers.size( 0 ), 21 );
  ASSERT_EQ( uzReceivers.size( 1 ), 9 );
  ASSERT_EQ( uzReceivers.size( 0 ), 21 );

  // check seismo content. The pressure and velocity values cannot be directly checked as the problem is too small.
  // Since the basis is linear, check that the seismograms are nonzero (for t>0) and the seismogram at the center is equal
  // to the average of the others.
  for( int i=0; i<21; i++ )
  {
    if( i > 0 )
    {
      ASSERT_TRUE( std::abs( pReceivers[i][8] ) > 0 );
      ASSERT_TRUE( std::abs( uxReceivers[i][8] ) > 0 );
      ASSERT_TRUE( std::abs( uyReceivers[i][8] ) > 0 );
      ASSERT_TRUE( std::abs( uzReceivers[i][8] ) > 0 );
    }
    double avgP = 0;
    double avgUx = 0;
    double avgUy = 0;
    double avgUz = 0;
    for( int r=0; r<8; r++ )
    {
      avgP += pReceivers[i][r];
      avgUx += uxReceivers[i][r];
      avgUy += uyReceivers[i][r];
      avgUz += uzReceivers[i][r];
    }
    avgP /=8.0;
    avgUx /=8.0;
    avgUy /=8.0;
    avgUz /=8.0;
    ASSERT_TRUE( std::abs( pReceivers[i][8] - avgP ) < 0.00001 );
    ASSERT_TRUE( std::abs( uxReceivers[i][8] - avgUx ) < 0.00001 );
    ASSERT_TRUE( std::abs( uyReceivers[i][8] - avgUy ) < 0.00001 );
    ASSERT_TRUE( std::abs( uzReceivers[i][8] - avgUz ) < 0.00001 );
  }
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geosx::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
