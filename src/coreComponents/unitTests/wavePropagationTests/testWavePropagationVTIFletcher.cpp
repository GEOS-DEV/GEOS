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
#include "physicsSolvers/wavePropagation/AcousticVTIFletcherWaveEquationSEM.hpp"

#include <gtest/gtest.h>

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

// This unit test checks the interpolation done to extract seismic traces from a wavefield.
// It computes a seismogram at a receiver co-located with the source and compares it to the surrounding receivers.
char const * xmlInput =
  R"xml(
  <Problem>
    <Solvers>
      <AcousticVTIFletcherSEM
        name="acousticVTIFletcherSolver"
        cflFactor="0.25"
        discretization="FE1"
        targetRegions="{ Region }"
        sourceCoordinates="{ { 50, 50, 50 } }"
        timeSourceFrequency="2"
        receiverCoordinates="{ { 0.1, 0.1, 0.1 }, { 0.1, 0.1, 99.9 }, { 0.1, 99.9, 0.1 }, { 0.1, 99.9, 99.9 },
                                { 99.9, 0.1, 0.1 }, { 99.9, 0.1, 99.9 }, { 99.9, 99.9, 0.1 }, { 99.9, 99.9, 99.9 },
                                { 50, 50, 50 } }"
        outputSeismoTrace="0"
        dtSeismoTrace="0.1"/>
    </Solvers>
    <Mesh>
      <InternalMesh
        name="mesh"
        elementTypes="{ C3D8 }"
        xCoords="{ 0, 100 }"
        yCoords="{ 0, 100 }"
        zCoords="{ 0, 100 }"
        nx="{ 1 }"
        ny="{ 1 }"
        nz="{ 1 }"
        cellBlockNames="{ cb }"/>
    </Mesh>
    <Events
      maxTime="1">
      <PeriodicEvent
        name="solverApplications"
        forceDt="0.1"
        targetExactStartStop="0"
        targetExactTimestep="0"
        target="/Solvers/acousticVTIFletcherSolver"/>
      <PeriodicEvent
        name="waveFieldNp1Collection"
        timeFrequency="0.1"
        targetExactTimestep="0"
        target="/Tasks/waveFieldNp1Collection" />
      <PeriodicEvent
        name="waveFieldNCollection"
        timeFrequency="0.1"
        targetExactTimestep="0"
        target="/Tasks/waveFieldNCollection" />
      <PeriodicEvent
        name="waveFieldNm1Collection"
        timeFrequency="0.1"
        targetExactTimestep="0"
        target="/Tasks/waveFieldNm1Collection" />
    </Events>
    <NumericalMethods>
      <FiniteElements>
        <FiniteElementSpace
          name="FE1"
          order="1"
          formulation="SEM"/>
      </FiniteElements>
    </NumericalMethods>
    <ElementRegions>
      <CellElementRegion
        name="Region"
        cellBlocks="{ cb }"
        materialList="{ nullModel }"/>
    </ElementRegions>
    <Constitutive>
      <NullModel
        name="nullModel"/>
    </Constitutive>
    <FieldSpecifications>
      <FieldSpecification
        name="initialPressure_p_N"
        initialCondition="1"
        setNames="{ all }"
        objectPath="nodeManager"
        fieldName="pressure_p_n"
        scale="0.0"/>
      <FieldSpecification
        name="initialPressure_p_Nm1"
        initialCondition="1"
        setNames="{ all }"
        objectPath="nodeManager"
        fieldName="pressure_p_nm1"
        scale="0.0"/>
      <FieldSpecification
        name="initialPressure_q_N"
        initialCondition="1"
        setNames="{ all }"
        objectPath="nodeManager"
        fieldName="pressure_q_n"
        scale="0.0"/>
      <FieldSpecification
        name="initialPressure_q_Nm1"
        initialCondition="1"
        setNames="{ all }"
        objectPath="nodeManager"
        fieldName="pressure_q_nm1"
        scale="0.0"/>
      <FieldSpecification
        name="cellVelocity"
        initialCondition="1"
        objectPath="ElementRegions/Region/cb"
        fieldName="acousticVelocity"
        scale="1500"
        setNames="{ all }"/>
      <FieldSpecification
        name="cellDensity"
        initialCondition="1"
        objectPath="ElementRegions/Region/cb"
        fieldName="acousticDensity"
        scale="1"
        setNames="{ all }"/>
      <FieldSpecification
        name="cellEpsilon"
        initialCondition="1"
        objectPath="ElementRegions/Region/cb"
        fieldName="acousticEpsilon"
        scale="0.24"
        setNames="{ all }"/>
      <FieldSpecification
        name="cellDelta"
        initialCondition="1"
        objectPath="ElementRegions/Region/cb"
        fieldName="acousticDelta"
        scale="0.1"
        setNames="{ all }"/>
      <FieldSpecification
        name="cellSigma"
        initialCondition="1"
        objectPath="ElementRegions/Region/cb"
        fieldName="acousticSigma"
        scale="0.75"
        setNames="{ all }"/>
      <FieldSpecification
        name="zposFreeSurface"
        objectPath="faceManager"
        fieldName="acousticFreeSurface"
        scale="0.0"
        setNames="{ zpos }"/>
    </FieldSpecifications>
    <Tasks>
      <PackCollection
        name="waveFieldNp1Collection"
        objectPath="nodeManager"
        fieldName="pressure_p_np1"/>
      <PackCollection
        name="waveFieldNCollection"
        objectPath="nodeManager"
        fieldName="pressure_p_n"/>
      <PackCollection
        name="waveFieldNm1Collection"
        objectPath="nodeManager"
        fieldName="pressure_p_nm1"/>
    </Tasks>
  </Problem>
  )xml";

class AcousticVTIFletcherWaveEquationSEMTest : public ::testing::Test
{
public:

  AcousticVTIFletcherWaveEquationSEMTest():
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
  AcousticVTIFletcherWaveEquationSEM * propagator;
};

real64 constexpr AcousticVTIFletcherWaveEquationSEMTest::time;
real64 constexpr AcousticVTIFletcherWaveEquationSEMTest::dt;
real64 constexpr AcousticVTIFletcherWaveEquationSEMTest::eps;

TEST_F( AcousticVTIFletcherWaveEquationSEMTest, SeismoTrace )
{

  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  propagator = &state.getProblemManager().getPhysicsSolverManager().getGroup< AcousticVTIFletcherWaveEquationSEM >( "acousticVTIFletcherSolver" );
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
  arrayView2d< real32 > const pReceivers = propagator->getReference< array2d< real32 > >( AcousticVTIFletcherWaveEquationSEM::viewKeyStruct::pressureNp1AtReceiversString() ).toView();

  // move it to CPU, if needed
  pReceivers.move( hostMemorySpace, false );

  // check number of seismos and trace length
  ASSERT_EQ( pReceivers.size( 1 ), 10 );
  ASSERT_EQ( pReceivers.size( 0 ), 11 );

  // check seismo content. The pressure values cannot be directly checked as the problem is too small.
  // Since the basis is linear, check that the seismograms are nonzero (for t>0) and the seismogram at the center is equal
  // to the average of the others.
  for( int i = 0; i < 11; i++ )
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
    avg /= 8.0;
    ASSERT_TRUE( std::abs( pReceivers[i][8] - avg ) < 0.00001 );
  }
  // (not supported yet) run adjoint solver
/*
  for( int i = 0; i < 10; i++ )
  {
    propagator->explicitStepBackward( time_n, dt, i, domain, false );
    time_n += dt;
  }
  // check again the seismo content.
  for( int i = 0; i < 11; i++ )
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
    avg /= 8.0;
    ASSERT_TRUE( std::abs( pReceivers[i][8] - avg ) < 0.00001 );
  }*/  
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
