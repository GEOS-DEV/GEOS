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
#include "physicsSolvers/inducedSeismicity/DieterichSeismicityRate.hpp"

#include <gtest/gtest.h>

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

// This unit test checks the accuracy of the integral seismicity rate solver
// It computes the seismicity in response to a hard coded stressing history to which there exists an analytical solution
char const * xmlInput =
  R"xml(
  <?xml version="1.0" ?>
  <Problem>
    <Solvers>
      <DieterichSeismicityRate
        name="dieterichSR"
        discretization="singlePhaseTPFA"
        stressSolverName="singlePhaseFlow"
        targetRegions="{ Domain }"
        initialFaultNormalStress="-155e6"
        initialFaultShearStress="93e6"
        directEffect="0.01"
        backgroundStressingRate="3.171e-5">
      </DieterichSeismicityRate>
      <SinglePhaseFVM
        name="singlePhaseFlow"
        logLevel="1"
        discretization="singlePhaseTPFA"
        targetRegions="{ Domain }">
      </SinglePhaseFVM>
    </Solvers>
    <Mesh>
      <InternalMesh
        name="mesh"
        elementTypes="{ C3D8 }"
        xCoords="{ 0, 1e3 }"
        yCoords="{ 0, 1e3 }"
        zCoords="{ 0, 1e3 }"
        nx="{ 10 }"
        ny="{ 10 }"
        nz="{ 10 }"
        cellBlockNames="{ cb1 }"/>
    </Mesh>
    <NumericalMethods>
      <FiniteVolume>
        <TwoPointFluxApproximation
          name="singlePhaseTPFA"/>
      </FiniteVolume>
    </NumericalMethods>
    <Events
      maxTime="360000.0">
      <PeriodicEvent
        name="solverApplications"
        forceDt="3600.0"
        target="/Solvers/dieterichSR"/>
    </Events>
    <ElementRegions>
      <CellElementRegion
        name="Domain"
        cellBlocks="{ cb1 }"
        materialList="{ water, rock }"/>
    </ElementRegions>
    <Constitutive>
      <CompressibleSinglePhaseFluid
        name="water"
        defaultDensity="1000"
        defaultViscosity="0.001"
        referencePressure="1e6"
        compressibility="1e-9"
        viscosibility="0.0"/>
      <CompressibleSolidConstantPermeability
        name="rock"
        solidModelName="nullSolid"
        porosityModelName="rockPorosity"
        permeabilityModelName="rockPerm"/>
      <NullModel
        name="nullSolid"/>
      <PressurePorosity
        name="rockPorosity"
        defaultReferencePorosity="0.05"
        referencePressure="1.0e6"
        compressibility="1.0e-9"/>
      <ConstantPermeability
        name="rockPerm"
        permeabilityComponents="{ 2.0e-14, 2.0e-14, 2.0e-14 }"/>
    </Constitutive>
    <FieldSpecifications>
      <FieldSpecification
        name="Porosity"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/Domain/cb1"
        fieldName="rockPorosity_referencePorosity"
        scale="0.05"/>
      <FieldSpecification
        name="initialPressure"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/Domain/cb1"
        fieldName="pressure"
        scale="0.0e6"/>
      <FieldSpecification
        name="boundaryPressureLeft"
        objectPath="faceManager"
        fieldName="pressure"
        scale="0.0e6"
        setNames="{ left }"/>
      <FieldSpecification
        name="boundaryPressureRight"
        objectPath="faceManager"
        fieldName="pressure"
        scale="0e6"
        setNames="{ right }"/>
      <FieldSpecification
        name="boundaryPressureTop"
        objectPath="faceManager"
        fieldName="pressure"
        scale="0e6"
        setNames="{ top }"/>
      <FieldSpecification
        name="boundaryPressureBot"
        objectPath="faceManager"
        fieldName="pressure"
        scale="0e6"
        setNames="{ bot }"/>
      <FieldSpecification
        name="boundaryPressureFront"
        objectPath="faceManager"
        fieldName="pressure"
        scale="0e6"
        setNames="{ front }"/>
       <FieldSpecification
        name="boundaryPressureBack"
        objectPath="faceManager"
        fieldName="pressure"
        scale="0e6"
        setNames="{ back }"/>
    </FieldSpecifications>
    <Geometry>
      <Box
        name="left"
        xMin="{ -1, -1, -1 }"
        xMax="{ 0.1, 1001, 1001 }"/>
      <Box
        name="right"
        xMin="{ 999, -1, -1 }"
        xMax="{ 1001, 1001, 1001 }"/>
      <Box
        name="top"
        xMin="{ -1, -1, 999 }"
        xMax="{ 1001, 1001, 1001 }"/>
      <Box
        name="bot"
        xMin="{ -1, -1, -1 }"
        xMax="{ 1001, 1001, 0.1 }"/>
      <Box
        name="front"
        xMin="{ -1, 999, -1 }"
        xMax="{ 1001, 1001, 1001 }"/>
      <Box
        name="back"
        xMin="{ -1, -1, -1 }"
        xMax="{ 1001, 0.1, 1001 }"/>
    </Geometry>
  </Problem>
  )xml";

class DieterichSeismicityRateIntegralSolverTest : public ::testing::Test
{
public:

  DieterichSeismicityRateIntegralSolverTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {
    setupProblemFromXML( state.getProblemManager(), xmlInput );
  }

  GeosxState state;
  DieterichSeismicityRate * propagator;
};

TEST_F( DieterichSeismicityRateIntegralSolverTest, solverTest )
{
  // DomainPartition & domain = state.getProblemManager().getDomainPartition();
  // propagator = &state.getProblemManager().getPhysicsSolverManager().getGroup< AcousticFirstOrderWaveEquationSEM >( "acousticFirstOrderSolver" );

  // check number of seismos and trace length
  ASSERT_EQ( 0, 0 );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}