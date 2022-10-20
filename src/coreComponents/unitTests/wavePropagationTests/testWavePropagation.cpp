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
#include "physicsSolvers/wavePropagation/AcousticWaveEquationSEM.hpp"

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
  "    <AcousticSEM\n"
  "      name=\"acousticSolver\"\n"
  "      cflFactor=\"0.25\"\n"
  "      discretization=\"FE1\"\n"
  "      targetRegions=\"{ Region }\"\n"
  "      sourceCoordinates=\"{ { 50, 50, 50 } }\"\n"
  "      timeSourceFrequency=\"2\"\n"
  "      receiverCoordinates=\"{ { 0.1, 0.1, 0.1 }, { 0.1, 0.1, 99.9 }, { 0.1, 99.9, 0.1 }, { 0.1, 99.9, 99.9 },\n"
  "                              { 99.9, 0.1, 0.1 }, { 99.9, 0.1, 99.9 }, { 99.9, 99.9, 0.1 }, { 99.9, 99.9, 99.9 },\n"
  "                              { 50, 50, 50 } }\"\n"
  "      outputSeismoTrace=\"0\"\n"
  "      dtSeismoTrace=\"0.1\"/>\n"
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
  "      forceDt=\"0.1\"\n"
  "      targetExactStartStop=\"0\"\n"
  "      targetExactTimestep=\"0\"\n"
  "      target=\"/Solvers/acousticSolver\"/>\n"
  "    <PeriodicEvent\n"
  "      name=\"waveFieldNp1Collection\"\n"
  "      timeFrequency=\"0.1\"\n"
  "      targetExactTimestep=\"0\"\n"
  "      target=\"/Tasks/waveFieldNp1Collection\" />\n"
  "    <PeriodicEvent\n"
  "      name=\"waveFieldNCollection\"\n"
  "      timeFrequency=\"0.1\"\n"
  "      targetExactTimestep=\"0\"\n"
  "      target=\"/Tasks/waveFieldNCollection\" />\n"
  "    <PeriodicEvent\n"
  "      name=\"waveFieldNm1Collection\"\n"
  "      timeFrequency=\"0.1\"\n"
  "      targetExactTimestep=\"0\"\n"
  "      target=\"/Tasks/waveFieldNm1Collection\" />\n"
  "  </Events>\n"
  "  <NumericalMethods>\n"
  "    <FiniteElements>\n"
  "      <FiniteElementSpace\n"
  "        name=\"FE1\"\n"
  "        order=\"1\"/>\n"
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
  "      name=\"initialPressureN\"\n"
  "      initialCondition=\"1\"\n"
  "      setNames=\"{ all }\"\n"
  "      objectPath=\"nodeManager\"\n"
  "      fieldName=\"pressure_n\"\n"
  "      scale=\"0.0\"/>\n"
  "    <FieldSpecification\n"
  "      name=\"initialPressureNm1\"\n"
  "      initialCondition=\"1\"\n"
  "      setNames=\"{ all }\"\n"
  "      objectPath=\"nodeManager\"\n"
  "      fieldName=\"pressure_nm1\"\n"
  "      scale=\"0.0\"/>\n"
  "    <FieldSpecification\n"
  "      name=\"cellVelocity\"\n"
  "      initialCondition=\"1\"\n"
  "      objectPath=\"ElementRegions/Region/cb\"\n"
  "      fieldName=\"mediumVelocity\"\n"
  "      scale=\"1500\"\n"
  "      setNames=\"{ all }\"/>\n"
  "    <FieldSpecification\n"
  "      name=\"zposFreeSurface\"\n"
  "      objectPath=\"faceManager\"\n"
  "      fieldName=\"FreeSurface\"\n"
  "      scale=\"0.0\"\n"
  "      setNames=\"{ zpos }\"/>\n"
  "  </FieldSpecifications>\n"
  "  <Tasks>\n"
  "    <PackCollection\n"
  "      name=\"waveFieldNp1Collection\"\n"
  "      objectPath=\"nodeManager\"\n"
  "      fieldName=\"pressure_np1\"/>\n"
  "    <PackCollection\n"
  "      name=\"waveFieldNCollection\"\n"
  "      objectPath=\"nodeManager\"\n"
  "      fieldName=\"pressure_n\"/>\n"
  "    <PackCollection\n"
  "      name=\"waveFieldNm1Collection\"\n"
  "      objectPath=\"nodeManager\"\n"
  "      fieldName=\"pressure_nm1\"/>\n"
  "  </Tasks>\n"
  "</Problem>\n";

class AcousticWaveEquationSEMTest : public ::testing::Test
{
public:

  AcousticWaveEquationSEMTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {
    setupProblemFromXML( state.getProblemManager(), xmlInput );
  }

  static real64 constexpr time = 0.0;
  static real64 constexpr dt = 1e-1;
  static real64 constexpr eps = std::numeric_limits< real64 >::epsilon();

  GeosxState state;
  AcousticWaveEquationSEM * propagator;
};

real64 constexpr AcousticWaveEquationSEMTest::time;
real64 constexpr AcousticWaveEquationSEMTest::dt;
real64 constexpr AcousticWaveEquationSEMTest::eps;

TEST_F( AcousticWaveEquationSEMTest, SeismoTrace )
{

  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  propagator = &state.getProblemManager().getPhysicsSolverManager().getGroup< AcousticWaveEquationSEM >( "acousticSolver" );
  real64 time_n = time;
  // run for 1s (10 steps)
  for( int i=0; i<10; i++ )
  {
    propagator->explicitStepForward( time_n, dt, i, domain, false );
    time_n += dt;
  }
  // cleanup (triggers calculation of the remaining seismograms data points)
  propagator->cleanup( 1.0, 10, 0, 0, domain );

  // retrieve seismo
  arrayView2d< real32 > const pReceivers = propagator->getReference< array2d< real32 > >( AcousticWaveEquationSEM::viewKeyStruct::pressureNp1AtReceiversString() ).toView();

  // move it to CPU, if needed
  pReceivers.move( LvArray::MemorySpace::host, false );

  // check number of seismos and trace length
  ASSERT_EQ( pReceivers.size( 1 ), 9 );
  ASSERT_EQ( pReceivers.size( 0 ), 11 );

  // check seismo content. The pressure values cannot be directly checked as the problem is too small.
  // Since the basis is linear, check that the seismograms are nonzero (for t>0) and the seismogram at the center is equal
  // to the average of the others.
  for( int i=0; i<11; i++ )
  {
    if( i > 0 )
    {
      ASSERT_TRUE( std::abs( pReceivers[i][8] ) > 0 );
    }
    double avg = 0;
    for( int r=0; r<8; r++ )
    {
      avg += pReceivers[i][r];
    }
    avg /=8.0;
    ASSERT_TRUE( std::abs( pReceivers[i][8] - avg ) < 0.00001 );
  }
  // run adjoint solver
  for( int i=0; i<10; i++ )
  {
    propagator->explicitStepBackward( time_n, dt, i, domain, false );
    time_n += dt;
  }
  // check again the seismo content.
  for( int i=0; i<11; i++ )
  {
    if( i > 0 )
    {
      ASSERT_TRUE( std::abs( pReceivers[i][8] ) > 0 );
    }
    double avg = 0;
    for( int r=0; r<8; r++ )
    {
      avg += pReceivers[i][r];
    }
    avg /=8.0;
    ASSERT_TRUE( std::abs( pReceivers[i][8] - avg ) < 0.00001 );
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
