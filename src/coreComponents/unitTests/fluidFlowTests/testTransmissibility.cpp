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


#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "codingUtilities/UnitTestUtilities.hpp"
#include "unitTests/fluidFlowTests/testSingleFlowUtils.hpp"
#include "physicsSolvers/fluidFlow/StencilDataCollection.hpp"
#include "mainInterface/ProblemManager.hpp"


using namespace geos;
using namespace geos::dataRepository;
using namespace geos::testing;


CommandLineOptions g_commandLineOptions;


// TODO remove everything that is not useful (wells? solver application?)
/// Provide every common xml input for the transmissibility tests
constexpr string_view xmlInputCommon =
  R"xml(
<Problem>
  <Solvers gravityVector="{ 0.0, 0.0, -9.81 }">
    <SinglePhaseReservoir name="reservoirSystem"
                flowSolverName="singlePhaseFlow"
                wellSolverName="singlePhaseWell"
                logLevel="1"
                targetRegions="{Region1,wellRegion1,wellRegion2}">
      <NonlinearSolverParameters newtonMaxIter="40"/>
      <LinearSolverParameters solverType="direct"
                              logLevel="1"/>
    </SinglePhaseReservoir>
    <SinglePhaseFVM name="singlePhaseFlow"
                    logLevel="1"
                    discretization="singlePhaseTPFA"
                    targetRegions="{Region1}">
    </SinglePhaseFVM>
    <SinglePhaseWell name="singlePhaseWell"
                      logLevel="1"
                      targetRegions="{wellRegion1,wellRegion2}">
        <WellControls name="wellControls1"
                      type="producer"
                      referenceElevation="2"
                      control="BHP"
                      targetBHP="5e5"
                      targetTotalRate="1e-3"/>
        <WellControls name="wellControls2"
                      type="injector"
                      referenceElevation="2"
                      control="totalVolRate"
                      targetBHP="2e7"
                      targetTotalRate="1e-4"/>
    </SinglePhaseWell>
  </Solvers>
  <NumericalMethods>
    <FiniteVolume>
      <TwoPointFluxApproximation name="singlePhaseTPFA"/>
    </FiniteVolume>
  </NumericalMethods>
  <ElementRegions>
    <CellElementRegion name="Region1"
                        cellBlocks="{cb1}"
                        materialList="{water, rock}"/>
    <WellElementRegion name="wellRegion1"
                        materialList="{water}"/>
    <WellElementRegion name="wellRegion2"
                        materialList="{water}"/>
  </ElementRegions>
  <Constitutive>
    <CompressibleSinglePhaseFluid name="water"
                                  defaultDensity="1000"
                                  defaultViscosity="0.001"
                                  referencePressure="0.0"
                                  referenceDensity="1000"
                                  compressibility="5e-10"
                                  referenceViscosity="0.001"
                                  viscosibility="0.0"/>
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
  </Constitutive>
  <FieldSpecifications>
    <FieldSpecification name="initialPressure"
                        initialCondition="1"
                        setNames="{all}"
                        objectPath="ElementRegions/Region1/cb1"
                        fieldName="pressure"
                        scale="5e6"/>
  </FieldSpecifications>

  <Events maxTime="3.0">
    <PeriodicEvent
      name="solverApplications"
      forceDt="10"
      target="/Solvers/reservoirSystem" />
    <SoloEvent
      name="cellToCellDataCollectionEvent"
      beginTime="1.0"
      target="/Tasks/cellToCellDataCollection" />
  </Events>

  <Tasks>
    <CellToCellDataCollection
      name="cellToCellDataCollection"
      flowSolverName="singlePhaseFlow"
      logLevel="1" />
  </Tasks>

)xml";

/// Provide the ending of the xml input for the transmissibility tests
constexpr string_view xmlInputEnd =
  R"xml(
</Problem>
)xml";

constexpr string_view stencilDataCollectionPath = "/Tasks/cellToCellDataCollection";


using Vector3 = std::array< real64, 3 >;
using TransmissibilityMap = std::map< std::pair< globalIndex, globalIndex >, real64 >;

enum class Axis : integer { X = 0, Y = 1, Z = 2 };

struct TestParams
{
  std::array< integer, 3 > cellCount;
  Vector3 cellDistance;
};

constexpr real64 g_transmissibilityTolerance = 1.0e-11;
constexpr Vector3 g_testPermeability = { 2.0e-16, 2.0e-16, 2.0e-16 };
constexpr real64 g_testNetToGross = 1.0;


void testStencilOutputStructured( string_view xmlInput, TestParams const & params );


TEST( TransmissibilityTest, stencilOutputVerificationIso )
{
  static constexpr string_view meshInput =
    R"xml(
  <Mesh>
    <InternalMesh name="mesh1"
                  elementTypes="{C3D8}"
                  xCoords="{0, 30}"
                  yCoords="{0, 30}"
                  zCoords="{0, 30}"
                  nx="{3}"
                  ny="{3}"
                  nz="{3}"
                  cellBlockNames="{cb1}">
      <InternalWell name="well_producer1"
                    wellRegionName="wellRegion1"
                    wellControlsName="wellControls1"
                    polylineNodeCoords="{ { 4.5, 0, 2   },
                                          { 4.5, 0, 0.5 } }"
                    polylineSegmentConn="{ { 0, 1 } }"
                    radius="0.1"
                    numElementsPerSegment="1">
          <Perforation name="producer1_perf1"
                        distanceFromHead="1.45"/>
      </InternalWell>
      <InternalWell name="well_injector1"
                    wellRegionName="wellRegion2"
                    wellControlsName="wellControls2"
                    polylineNodeCoords="{ { 0.5, 0, 2   },
                                          { 0.5, 0, 0.5 } }"
                    polylineSegmentConn="{ { 0, 1 } }"
                    radius="0.1"
                    numElementsPerSegment="1">
          <Perforation name="injector1_perf1"
                        distanceFromHead="1.45"/>
      </InternalWell>
    </InternalMesh>
  </Mesh>
)xml";
  std::ostringstream xmlInput;
  xmlInput << xmlInputCommon << meshInput << xmlInputEnd;

  TestParams const params = {
    { 3, 3, 3 }, // cellCount
    { 10.0, 10.0, 10.0 }, // cellDistance
  };

  testStencilOutputStructured( xmlInput.str(), params );
}

// TEST_F( TransmissibilityTest, StencilOutputVerificationAniso )
// {
//   ...
//
//    testStencilOutput( { 2.0, 3.0, 4.0 } );
// }



/**
 * @return The theorical half transmissiblity (from A to B or B to A)
 * @param params
 * @param axis The axis in which we want to compute the transmissibility
 */
real64 computeTransmissiblityStructured( TestParams const & params, Axis axis )
{
  real64 const faceArea = axis == Axis::X ? params.cellDistance[1]*params.cellDistance[2]:
                          axis == Axis::Y ? params.cellDistance[0]*params.cellDistance[2]:
                          params.cellDistance[0]*params.cellDistance[1];
  real64 const halfDistance = 0.5 * params.cellDistance[integer( axis )];

  real64 const transmissibility = g_testPermeability[integer( axis )] * g_testNetToGross * ( faceArea / halfDistance );
  // return half transmissibility
  return g_testNetToGross * transmissibility * 0.5;
}

void verifyConnectionsByDim( TransmissibilityMap transmissibilities, Axis axis,
                             TestParams const & params )
{
  integer const cellBIdOffset = axis == Axis::X ? 1:
                                axis == Axis::Y ? params.cellCount[0]:
                                params.cellCount[0]*params.cellCount[1];
  integer const endX = axis == Axis::X ? params.cellCount[0] - 1 : params.cellCount[0];
  integer const endY = axis == Axis::Y ? params.cellCount[1] - 1 : params.cellCount[1];
  integer const endZ = axis == Axis::Z ? params.cellCount[2] - 1 : params.cellCount[2];
  real64 const expectedAxisTransmissibility = computeTransmissiblityStructured( params, axis );

  for( integer z = 0; z < endZ; ++z )
  {
    for( integer y = 0; y < endY; ++y )
    {
      for( integer x = 0; x < endX; ++x )
      {
        integer const cellAId = x + params.cellCount[0] * (y + z * params.cellCount[1]);
        integer const cellBId = cellAId + cellBIdOffset;
        real64 const transmissibilityOutput = transmissibilities[std::make_pair( cellAId, cellBId )];
        real64 const transmissibilityTolerance = g_transmissibilityTolerance * std::abs( transmissibilityOutput );
        EXPECT_NEAR( transmissibilityOutput, expectedAxisTransmissibility,
                     transmissibilityTolerance )<< GEOS_FMT( "Transmissibility data from {} does not match with expectation.",
                                                             StencilDataCollection::catalogName() );
      }
    }
  }
}

TransmissibilityMap getTransmissibilityMap( StencilDataCollection const & stencilData,
                                            TestParams const & params )
{
  using VK = StencilDataCollection::viewKeyStruct;
  auto const & cellAGlobalId = stencilData.getReference< array1d< globalIndex > >( VK::cellAGlobalIdString() );
  auto const & cellBGlobalId = stencilData.getReference< array1d< globalIndex > >( VK::cellBGlobalIdString() );
  auto const & transmissibilityAB = stencilData.getReference< array1d< real64 > >( VK::transmissibilityABString() );
  auto const & transmissibilityBA = stencilData.getReference< array1d< real64 > >( VK::transmissibilityBAString() );

  // verify data size
  integer const nx = params.cellCount[0];
  integer const ny = params.cellCount[1];
  integer const nz = params.cellCount[2];
  integer const arraySize = transmissibilityAB.size();
  EXPECT_EQ( arraySize, 3*nx*ny*nz - nx*ny - ny*nz - nx*nz );
  // verify that array size is always equal for all buffers
  EXPECT_EQ( cellAGlobalId.size(), arraySize );
  EXPECT_EQ( cellBGlobalId.size(), arraySize );
  EXPECT_EQ( transmissibilityBA.size(), arraySize );

  TransmissibilityMap transmissibilities;
  for( int i = 0; i < arraySize; ++i )
  {
    transmissibilities[std::make_pair( cellAGlobalId[i], cellBGlobalId[i] )] = transmissibilityAB[i];
  }
  return transmissibilities;
}


void testStencilOutputStructured( string_view xmlInput, TestParams const & params )
{
  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  ProblemManager & problem = state.getProblemManager();
  setupProblemFromXML( problem, xmlInput.data() );

  EXPECT_FALSE( problem.runSimulation() ) << "Simulation exited early.";

  { // checking results
    StencilDataCollection const & stencilData =
      problem.getGroupByPath< StencilDataCollection >( string( stencilDataCollectionPath ) );
    TransmissibilityMap const transmissibilities = getTransmissibilityMap( stencilData, params );

    // let's check each couple of cell in the x direction
    verifyConnectionsByDim( transmissibilities, Axis::X, params );
    verifyConnectionsByDim( transmissibilities, Axis::Y, params );
    verifyConnectionsByDim( transmissibilities, Axis::Z, params );
  }
}



int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
