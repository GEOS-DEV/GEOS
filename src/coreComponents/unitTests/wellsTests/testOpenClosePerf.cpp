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

#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"

#include "common/DataTypes.hpp"
#include "mainInterface/initialization.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/multiphysics/CompositionalMultiphaseReservoirAndWells.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWell.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/CompositionalMultiphaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::constitutive;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

char const * PreXmlInput =
  R"xml(
  <Problem>
    <Solvers gravityVector="{ 0.0, 0.0, -9.81 }">
      <CompositionalMultiphaseReservoir name="reservoirSystem"
                 flowSolverName="compositionalMultiphaseFlow"
                 wellSolverName="compositionalMultiphaseWell"
                 logLevel="1"
                 targetRegions="{Region1,wellRegion1}">
        <NonlinearSolverParameters newtonMaxIter="40"/>
        <LinearSolverParameters solverType="direct"
                                logLevel="2"/>
      </CompositionalMultiphaseReservoir>
      <CompositionalMultiphaseFVM name="compositionalMultiphaseFlow"
                                  logLevel="1"
                                  discretization="fluidTPFA"
                                  targetRegions="{Region1}"
                                  temperature="297.15"
                                  useMass="0">
      </CompositionalMultiphaseFVM>
      <CompositionalMultiphaseWell name="compositionalMultiphaseWell"
                                   logLevel="1"
                                   targetRegions="{wellRegion1}"
                                   useMass="0">
          <WellControls name="wellControls1"
                        type="producer"
                        referenceElevation="1.25"
                        control="BHP"
                        targetBHP="2e6"
                        targetPhaseRate="1"
                        targetPhaseName="oil"/>
      </CompositionalMultiphaseWell>
    </Solvers>
    <Mesh>
      <InternalMesh name="mesh1"
                    elementTypes="{C3D8}"
                    xCoords="{ 0, 100 }"
                    yCoords="{ 0, 100 }"
                    zCoords="{ 0, 5 }"
                    nx="{ 1 }"
                    ny="{ 1 }"
                    nz="{ 5 }"
                    cellBlockNames="{cb1}">)xml";

char const * PostXmlInput =
  R"xml(
      </InternalMesh>
    </Mesh>
    <NumericalMethods>
      <FiniteVolume>
        <TwoPointFluxApproximation name="fluidTPFA"/>
      </FiniteVolume>
    </NumericalMethods>
    <ElementRegions>
      <CellElementRegion name="Region1"
                         cellBlocks="{cb1}"
                         materialList="{fluid1, rock, relperm}"/>
      <WellElementRegion name="wellRegion1"
                         materialList="{fluid1, relperm}"/>
    </ElementRegions>
    <Constitutive>
      <CompositionalMultiphaseFluid name="fluid1"
                                    phaseNames="{oil, gas}"
                                    equationsOfState="{PR, PR}"
                                    componentNames="{N2, C10, C20, H2O}"
                                    componentCriticalPressure="{34e5, 25.3e5, 14.6e5, 220.5e5}"
                                    componentCriticalTemperature="{126.2, 622.0, 782.0, 647.0}"
                                    componentAcentricFactor="{0.04, 0.443, 0.816, 0.344}"
                                    componentMolarWeight="{28e-3, 134e-3, 275e-3, 18e-3}"
                                    componentVolumeShift="{0, 0, 0, 0}"
                                    componentBinaryCoeff="{ {0, 0, 0, 0},
                                                            {0, 0, 0, 0},
                                                            {0, 0, 0, 0},
                                                            {0, 0, 0, 0} }"/>
      <CompressibleSolidConstantPermeability name="rock"
          solidModelName="nullSolid"
          porosityModelName="rockPorosity"
          permeabilityModelName="rockPerm"/>
     <NullModel name="nullSolid"/>
     <PressurePorosity name="rockPorosity"
                       defaultReferencePorosity="0.05"
                       referencePressure = "0.0"
                       compressibility="1.0e-9"/>
    <ConstantPermeability name="rockPerm"
                          permeabilityComponents="{2.0e-16, 2.0e-16, 2.0e-16}"/>
      <BrooksCoreyRelativePermeability name="relperm"
                                       phaseNames="{oil, gas}"
                                       phaseMinVolumeFraction="{0.1, 0.15}"
                                       phaseRelPermExponent="{2.0, 2.0}"
                                       phaseRelPermMaxValue="{0.8, 0.9}"/>
    </Constitutive>
    <FieldSpecifications>
      <FieldSpecification name="initialPressure"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/Region1/cb1"
                 fieldName="pressure"
                 scale="5e6"/>
      <FieldSpecification name="initialComposition_N2"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/Region1/cb1"
                 fieldName="globalCompFraction"
                 component="0"
                 scale="0.099"/>
      <FieldSpecification name="initialComposition_C10"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/Region1/cb1"
                 fieldName="globalCompFraction"
                 component="1"
                 scale="0.3"/>
      <FieldSpecification name="initialComposition_C20"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/Region1/cb1"
                 fieldName="globalCompFraction"
                 component="2"
                 scale="0.6"/>
      <FieldSpecification name="initialComposition_H20"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/Region1/cb1"
                 fieldName="globalCompFraction"
                 component="3"
                 scale="0.001"/>
    </FieldSpecifications>
  </Problem>
  )xml";



template< typename LAMBDA >
void testPlugTopDownPerfCheck( CompositionalMultiphaseReservoirAndWells<> & solver,
                               DomainPartition & domain,

                               LAMBDA && perfFunction )
{
  CompositionalMultiphaseWell & wellSolver = *solver.wellSolver();

  typedef std::map< real64, std::vector< int > > map_type;
  map_type refVal;
  refVal[29800.0] = { 1, 1, 1, 1, 1};
  refVal[32400.0] = { 1, 1, 1, 1, 0};
  refVal[32500.0] = { 0, 0, 0, 0, 0};
  refVal[33000.0] = { 1, 1, 1, 1, 1};

  map_type::iterator it = refVal.begin();
  while( it != refVal.end() )
  {
    perfFunction( it->first );

    const std::vector< int > & refStatus = it->second;
    wellSolver.forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                             MeshLevel & mesh,
                                                                             string_array const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                          [&]( localIndex const,
                                                                               WellElementSubRegion & subRegion )
      {
        array1d< integer > & currentStatus = subRegion.getWellElementStatus();
        EXPECT_EQ( refStatus.size(), currentStatus.size() );
        const integer * curStat = currentStatus.data();
        for( integer i = 0; i < currentStatus.size(); ++i )
        {
          EXPECT_EQ( curStat[ i ], refStatus[ i ] );
        }
      } );
    } );
    it++;
  }
}

template< typename LAMBDA >
void testPlugBottomUpPerfCheck( CompositionalMultiphaseReservoirAndWells<> & solver,
                                DomainPartition & domain,

                                LAMBDA && perfFunction )
{
  CompositionalMultiphaseWell & wellSolver = *solver.wellSolver();

  typedef std::map< real64, std::vector< int > > map_type;
  map_type refVal;
  refVal[4800.0] =  { 1, 1, 1, 1, 1};
  refVal[14800.0] = { 1, 1, 1, 1, 0};
  refVal[19600.0] = { 1, 1, 1, 0, 0};
  refVal[32400.0] = { 1, 1, 0, 0, 0};
  refVal[39800.0] = { 0, 0, 0, 0, 0};

  map_type::iterator it = refVal.begin();
  while( it != refVal.end() )
  {
    perfFunction( it->first );

    const std::vector< int > & refStatus = it->second;
    wellSolver.forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                             MeshLevel & mesh,
                                                                             string_array const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                          [&]( localIndex const,
                                                                               WellElementSubRegion & subRegion )
      {
        array1d< integer > & currentStatus = subRegion.getWellElementStatus();
        EXPECT_EQ( refStatus.size(), currentStatus.size() );
        const integer * curStat = currentStatus.data();
        for( integer i = 0; i < currentStatus.size(); ++i )
        {
          EXPECT_EQ( curStat[ i ], refStatus[ i ] );
        }
      } );
    } );
    it++;
  }
}

template< typename LAMBDA >
void testOpenTopDownPerfCheck( CompositionalMultiphaseReservoirAndWells<> & solver,
                               DomainPartition & domain,

                               LAMBDA && perfFunction )
{
  CompositionalMultiphaseWell & wellSolver = *solver.wellSolver();

  typedef std::map< real64, std::vector< int > > map_type;
  map_type refVal;
  refVal[4800.0] = { 0, 0, 0, 0, 0};
  refVal[14800.0] = { 1, 0, 0, 0, 0};
  refVal[19800.0] = { 1, 1, 0, 0, 0};
  refVal[29800.0] = { 1, 1, 1, 0, 0};
  refVal[39800.0] = { 1, 1, 1, 1, 1};
  map_type::iterator it = refVal.begin();
  while( it != refVal.end() )
  {
    perfFunction( it->first );

    const std::vector< int > & refStatus = it->second;
    wellSolver.forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                             MeshLevel & mesh,
                                                                             string_array const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                          [&]( localIndex const,
                                                                               WellElementSubRegion & subRegion )
      {
        array1d< integer > & currentStatus = subRegion.getWellElementStatus();
        EXPECT_EQ( refStatus.size(), currentStatus.size() );
        const integer * curStat = currentStatus.data();
        for( integer i = 0; i < currentStatus.size(); ++i )
        {
          EXPECT_EQ( curStat[ i ], refStatus[ i ] );
        }
      } );
    } );
    it++;
  }
}

template< typename LAMBDA >
void testOpenBottomUpPerfCheck( CompositionalMultiphaseReservoirAndWells<> & solver,
                                DomainPartition & domain,

                                LAMBDA && perfFunction )
{
  CompositionalMultiphaseWell & wellSolver = *solver.wellSolver();

  typedef std::map< real64, std::vector< int > > map_type;
  map_type refVal;
  refVal[4800.0] = { 0, 0, 0, 0, 0};
  refVal[33000.0] = { 1, 1, 1, 1, 1};

  map_type::iterator it = refVal.begin();
  while( it != refVal.end() )
  {
    perfFunction( it->first );

    const std::vector< int > & refStatus = it->second;
    wellSolver.forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                             MeshLevel & mesh,
                                                                             string_array const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                          [&]( localIndex const,
                                                                               WellElementSubRegion & subRegion )
      {
        array1d< integer > & currentStatus = subRegion.getWellElementStatus();
        EXPECT_EQ( refStatus.size(), currentStatus.size() );
        const integer * curStat = currentStatus.data();
        for( integer i = 0; i < currentStatus.size(); ++i )
        {
          EXPECT_EQ( curStat[ i ], refStatus[ i ] );
        }
      } );
    } );
    it++;
  }
}

class CompositionalMultiphaseReservoirSolverTest : public ::testing::Test
{
public:

  CompositionalMultiphaseReservoirSolverTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {}
  char const * plugTopDown  =
    R"(
    <InternalWell name="well_producer1"
    name="well_producer1"
    wellRegionName="wellRegion1"
    wellControlsName="wellControls1"
  logLevel="1"
    polylineNodeCoords="{ { 0.5, 0.5, 4.9 },
                          { 0.5, 0.5, 0.1 } }"
    polylineSegmentConn="{ { 0, 1 } }"
    radius="0.1"
    numElementsPerSegment="5">
    <Perforation
      name="producer1_perf0"
      perfStatusTable="{{0, 0.5e4},{1.0,0.0}}"
      distanceFromHead="0.5"/>
    <Perforation
      name="producer1_perf1"
      perfStatusTable="{{0, 1.50e4},{1.0,0.0}}"
      distanceFromHead="1.50"/>
    <Perforation
      name="producer1_perf2"
      perfStatusTable="{{0, 2.00e4},{1.0,0.0}}"
      distanceFromHead="2.50"/>
    <Perforation
      name="producer1_perf3"
      perfStatusTable="{{0, 3.25e4},{1.0,0.0}}"
      distanceFromHead="3.50"/>
    <Perforation
      name="producer1_perf4"
      perfStatusTable="{{0, 3.0e4,3.3e4},{1.0,0.0,1.0}}"
      distanceFromHead="4.00"/>
    </InternalWell>)";

  char const * plugBottomUp  =
    R"(
      <InternalWell
        name="well_producer1"
        wellRegionName="wellRegion1"
        wellControlsName="wellControls1"
      logLevel="1"
        polylineNodeCoords="{ { 0.5, 0.5, 4.9 },
                              { 0.5, 0.5, 0.1 } }"
        polylineSegmentConn="{ { 0, 1 } }"
        radius="0.1"
        numElementsPerSegment="5">
        <Perforation
          name="producer1_perf0"
          perfStatusTable="{{0, 3.0e4},{1.0,0.0}}"
          distanceFromHead="0.5"/>
        <Perforation
          name="producer1_perf1"
          perfStatusTable="{{0, 3.25e4},{1.0,0.0}}"
          distanceFromHead="1.50"/>
        <Perforation
          name="producer1_perf2"
          perfStatusTable="{{0, 2.00e4},{1.0,0.0}}"
          distanceFromHead="2.50"/>
        <Perforation
          name="producer1_perf3"
          perfStatusTable="{{0, 1.50e4},{1.0,0.0}}"
          distanceFromHead="3.50"/>
        <Perforation
          name="producer1_perf4"
          perfStatusTable="{{0, 0.5e4},{1.0,0.0}}"
          distanceFromHead="4.00"/>
      </InternalWell>)";

  char const * openTopDown  =
    R"(
      <InternalWell
        name="well_producer1"
        wellRegionName="wellRegion1"
        wellControlsName="wellControls1"
      logLevel="1"
        polylineNodeCoords="{ { 0.5, 0.5, 4.9 },
                              { 0.5, 0.5, 0.1 } }"
        polylineSegmentConn="{ { 0, 1 } }"
        radius="0.1"
        numElementsPerSegment="5">
        <Perforation
          name="producer1_perf0"
          perfStatusTable="{{0, 0.5e4},{0.0,1.0}}"
          distanceFromHead="0.5"/>
        <Perforation
          name="producer1_perf1"
          perfStatusTable="{{0, 1.50e4},{0.0,1.0}}"
          distanceFromHead="1.50"/>
        <Perforation
          name="producer1_perf2"
          perfStatusTable="{{0, 2.00e4},{0.0,1.0}}"
          distanceFromHead="2.50"/>
        <Perforation
          name="producer1_perf3"
          perfStatusTable="{{0, 3.25e4},{0.0,1.0}}"
          distanceFromHead="3.50"/>
        <Perforation
          name="producer1_perf4"
          perfStatusTable="{{0, 3.0e4},{0.0,1.0}}"
          distanceFromHead="4.00"/>
      </InternalWell>)";

  char const * openBottomUp  =
    R"(
      <InternalWell
        name="well_producer1"
        wellRegionName="wellRegion1"
        wellControlsName="wellControls1"
      logLevel="1"
        polylineNodeCoords="{ { 0.5, 0.5, 4.9 },
                              { 0.5, 0.5, 0.1 } }"
        polylineSegmentConn="{ { 0, 1 } }"
        radius="0.1"
        numElementsPerSegment="5">
        <Perforation
          name="producer1_perf0"
          perfStatusTable="{{0, 3.0e4},{0.0,1.0}}"
          distanceFromHead="0.5"/>
        <Perforation
          name="producer1_perf1"
          perfStatusTable="{{0, 3.25e4},{0.0,1.0}}"
          distanceFromHead="1.50"/>
        <Perforation
          name="producer1_perf2"
          perfStatusTable="{{0, 2.00e4},{0.0,1.0}}"
          distanceFromHead="2.50"/>
        <Perforation
          name="producer1_perf3"
          perfStatusTable="{{0, 1.50e4},{0.0,1.0}}"
          distanceFromHead="3.50"/>
        <Perforation
          name="producer1_perf4"
          perfStatusTable="{{0, 0.5e4},{0.0,1.0}}"
          distanceFromHead="4.00"/>
      </InternalWell>)";

  void SetupPerfs( char const * internalWells )
  {

    string const xmlInput = std::string( PreXmlInput ) + internalWells + PostXmlInput;


    setupProblemFromXML( state.getProblemManager(), xmlInput.c_str() );
    solver = &state.getProblemManager().getPhysicsSolverManager().getGroup< CompositionalMultiphaseReservoirAndWells<> >( "reservoirSystem" );

    DomainPartition & domain = state.getProblemManager().getDomainPartition();

    solver->setupSystem( domain,
                         solver->getDofManager(),
                         solver->getLocalMatrix(),
                         solver->getSystemRhs(),
                         solver->getSystemSolution() );

  }

  GeosxState state;
  CompositionalMultiphaseReservoirAndWells<> * solver;
};


TEST_F( CompositionalMultiphaseReservoirSolverTest, plugTopDownPerfCheck )
{

  SetupPerfs( plugTopDown );
  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testPlugTopDownPerfCheck( *solver, domain,
                            [&] ( real64 time )
  {
    WellSolverBase * wellSolverBase = solver->wellSolver();
    wellSolverBase->setPerforationStatus( time, domain );
  } );
}

TEST_F( CompositionalMultiphaseReservoirSolverTest, plugBottomUpPerfCheck )
{

  SetupPerfs( plugBottomUp );
  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testPlugBottomUpPerfCheck( *solver, domain,
                             [&] ( real64 time )
  {
    WellSolverBase * wellSolverBase = solver->wellSolver();
    wellSolverBase->setPerforationStatus( time, domain );
  } );
}
TEST_F( CompositionalMultiphaseReservoirSolverTest, openTopDownPerfCheck )
{

  SetupPerfs( openTopDown );
  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testOpenTopDownPerfCheck( *solver, domain,
                            [&] ( real64 time )
  {
    WellSolverBase * wellSolverBase = solver->wellSolver();
    wellSolverBase->setPerforationStatus( time, domain );
  } );
}

TEST_F( CompositionalMultiphaseReservoirSolverTest, openBottomUpPerfCheck )
{

  SetupPerfs( openBottomUp );
  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testOpenBottomUpPerfCheck( *solver, domain,
                             [&] ( real64 time )
  {
    WellSolverBase * wellSolverBase = solver->wellSolver();
    wellSolverBase->setPerforationStatus( time, domain );
  } );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
