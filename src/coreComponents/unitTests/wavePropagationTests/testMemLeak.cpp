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

// using some utility classes from the following unit test
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"

#include <iostream>
#include <cstdio>
#include <sys/resource.h>
#include "common/DataTypes.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mesh/ElementType.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/wavePropagation/shared/WaveSolverBase.hpp"
#include "physicsSolvers/wavePropagation/sem/acoustic/secondOrderEqn/isotropic/AcousticWaveEquationSEM.hpp"
#include "events/EventManager.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"

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
      <AcousticSEM
        name="acousticSolver"
        discretization="FE1"
        targetRegions="{ Region }"
        sourceCoordinates="{ { 50, 50, 50 } }"
        timeSourceFrequency="2"
        receiverCoordinates="{ { 0.1, 0.1, 0.1 }, { 0.1, 0.1, 99.9 }, { 0.1, 99.9, 0.1 }, { 0.1, 99.9, 99.9 },
                                { 99.9, 0.1, 0.1 }, { 99.9, 0.1, 99.9 }, { 99.9, 99.9, 0.1 }, { 99.9, 99.9, 99.9 },
                                { 50, 50, 50 } }"
        outputSeismoTrace="0"
        timestepStabilityLimit="1"
        cflFactor="0.95"
        dtSeismoTrace="0.1"/>
    </Solvers>
    <Mesh>
      <InternalMesh
        name="mesh"
        elementTypes="{ C3D8 }"
        xCoords="{ 0, 100 }"
        yCoords="{ 0, 100 }"
        zCoords="{ 0, 100 }"
        nx="{ 10 }"
        ny="{ 10 }"
        nz="{ 10 }"
        cellBlockNames="{ cb }"/>
    </Mesh>
    <Events
      maxTime="1">
      <PeriodicEvent
        name="solverApplications"
        forceDt="0.1"
        targetExactStartStop="0"
        targetExactTimestep="0"
        target="/Solvers/acousticSolver"/>
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
        cellBlocks="{ * }"
        materialList="{ nullModel }"/>
    </ElementRegions>
    <Constitutive>
      <NullModel
        name="nullModel"/>
    </Constitutive>
    <FieldSpecifications>
      <FieldSpecification
        name="initialPressureN"
        initialCondition="1"
        setNames="{ all }"
        objectPath="nodeManager"
        fieldName="pressure_n"
        scale="0.0"/>
      <FieldSpecification
        name="initialPressureNm1"
        initialCondition="1"
        setNames="{ all }"
        objectPath="nodeManager"
        fieldName="pressure_nm1"
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
        name="zposFreeSurface"
        objectPath="faceManager"
        fieldName="FreeSurface"
        scale="0.0"
        setNames="{ zpos }"/>
    </FieldSpecifications>
    <Tasks>
      <PackCollection
        name="waveFieldNp1Collection"
        objectPath="nodeManager"
        fieldName="pressure_np1"/>
      <PackCollection
        name="waveFieldNCollection"
        objectPath="nodeManager"
        fieldName="pressure_n"/>
      <PackCollection
        name="waveFieldNm1Collection"
        objectPath="nodeManager"
        fieldName="pressure_nm1"/>
    </Tasks>
  </Problem>
  )xml";

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

  GeosxState state;
  AcousticWaveEquationSEM * propagator;
};

TEST_F( AcousticWaveEquationSEMTest, SeismoTrace )
{

  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  propagator = &state.getProblemManager().getPhysicsSolverManager().getGroup< AcousticWaveEquationSEM >( "acousticSolver" );

  int rnk = MpiWrapper::commRank( MPI_COMM_GEOS );

  printf( "                             usage0    usage1    usage2    usage3    usage4    usage5    usage6 \n");
  for(int i = 0; i < 1000000; i++ )
  {
    struct rusage usage0;
    struct rusage usage1;
    struct rusage usage2;
    struct rusage usage3;
    struct rusage usage4;
    struct rusage usage5;
    struct rusage usage6;

    getrusage(RUSAGE_SELF, &usage0);

    FieldIdentifiers fieldsToBeSync;
    getrusage(RUSAGE_SELF, &usage1);
    fieldsToBeSync.addFields( FieldLocation::Node, { fields::acousticfields::Pressure_np1::key() } );
    getrusage(RUSAGE_SELF, &usage2);
    CommunicationTools & syncFields = CommunicationTools::getInstance();

    getrusage(RUSAGE_SELF, &usage3);

    MeshLevel * meshLevel;
    propagator->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                              MeshLevel & mesh,
                                                                              string_array const & )
    {
      meshLevel = &mesh;
    } );
    MeshLevel & mesh = *meshLevel;
  
    {

      getrusage(RUSAGE_SELF, &usage4);
      syncFields.synchronizeFields( fieldsToBeSync,
                                  mesh,
                                  domain.getNeighbors(),
                                  true );

      getrusage(RUSAGE_SELF, &usage5);
    }// );

    getrusage(RUSAGE_SELF, &usage6);

    if( usage0.ru_maxrss != usage1.ru_maxrss ||
        usage1.ru_maxrss != usage2.ru_maxrss ||
        usage2.ru_maxrss != usage3.ru_maxrss ||
        usage3.ru_maxrss != usage4.ru_maxrss ||
        usage4.ru_maxrss != usage5.ru_maxrss ||
        usage5.ru_maxrss != usage6.ru_maxrss )
    {
      printf( "Cycle %6d rank %6d Memory usage: %ld  %ld  %ld  %ld  %ld  %ld  %ld \n", 
              i,
              rnk, 
              usage0.ru_maxrss,
              usage1.ru_maxrss,
              usage2.ru_maxrss,
              usage3.ru_maxrss,
              usage4.ru_maxrss,
              usage5.ru_maxrss,
              usage6.ru_maxrss );
    }
  }
  ASSERT_TRUE( true );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
