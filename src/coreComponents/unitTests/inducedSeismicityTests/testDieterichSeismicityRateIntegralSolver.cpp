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

static constexpr real64 relTol = 1.0e-2;

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
        initialFaultNormalTraction="-100e6"
        initialFaultShearTraction="60e6"
        directEffect="0.01"
        backgroundStressingRate="3.171e-5">
      </DieterichSeismicityRate>
      <SinglePhaseFVM
        name="singlePhaseFlow"
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

  static real64 constexpr time = 0.0;
  static real64 constexpr dt = 3600.0;
  static real64 constexpr eps = std::numeric_limits< real64 >::epsilon();

  // Parameters of model and analytical solution
  static real64 constexpr cTau = 1e-3;
  static real64 constexpr aSigma = 1e6;
  static real64 constexpr t_a = 1e6/3.171e-5;
  static real64 constexpr iniSig = -100e6;
  static real64 constexpr iniTau = 60e6;

  GeosxState state;
  DieterichSeismicityRate * propagator;
};

TEST_F( DieterichSeismicityRateIntegralSolverTest, solverTest )
{
  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  propagator = &state.getProblemManager().getPhysicsSolverManager().getGroup< DieterichSeismicityRate >( "dieterichSR" );
  real64 time_n = time;

  // run for 100 hours (100 steps)
  for( int i=0; i<100; i++ )
  {

    // Loop through CellElementSubRegions of domain to pass the subRegion as input to seismciity rate solver
    domain.forMeshBodies( [&] ( MeshBody & meshBody )
    {
      meshBody.forMeshLevels( [&] ( MeshLevel & mesh )
      {
        ElementRegionManager & elemManager = mesh.getElemManager();
        for( localIndex er = 0; er < elemManager.numRegions(); ++er )
        {
          ElementRegionBase & elemRegion = elemManager.getRegion( er );
          elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const, CellElementSubRegion & subRegion )
          {

            // Initialize shear and normal stresses acting on fault
            // (field variables in seismicity rate solver) as specified in XML file
            if( i==0 )
            {
              arrayView1d< real64 > const tempSigIni = subRegion.getField< inducedSeismicity::initialProjectedNormalTraction >();
              tempSigIni.setValues< parallelHostPolicy >( iniSig );
              arrayView1d< real64 > const tempSig_n = subRegion.getField< inducedSeismicity::projectedNormalTraction_n >();
              tempSig_n.setValues< parallelHostPolicy >( iniSig );
              arrayView1d< real64 > const tempSig = subRegion.getField< inducedSeismicity::projectedNormalTraction >();
              tempSig.setValues< parallelHostPolicy >( iniSig );

              arrayView1d< real64 > const tempTauIni = subRegion.getField< inducedSeismicity::initialProjectedShearTraction >();
              tempTauIni.setValues< parallelHostPolicy >( iniTau );
              arrayView1d< real64 > const tempTau_n = subRegion.getField< inducedSeismicity::projectedShearTraction_n >();
              tempTau_n.setValues< parallelHostPolicy >( iniTau );
              arrayView1d< real64 > const tempTau = subRegion.getField< inducedSeismicity::projectedShearTraction >();
              tempTau.setValues< parallelHostPolicy >( iniTau );
            }

            // Retrieve shear stress fields from seismicity rate solver, to be hard coded
            arrayView1d< real64 const > const tau_i = subRegion.getField< geos::fields::inducedSeismicity::initialProjectedShearTraction >();
            arrayView1d< real64 > const tau = subRegion.getField< geos::fields::inducedSeismicity::projectedShearTraction >();
            arrayView1d< real64 > const tau_n = subRegion.getField< geos::fields::inducedSeismicity::projectedShearTraction_n >();

            // Shear stress history will follow a log history,
            // to which one can compute an analytical solution for the seismicity rate
            real64 curTau = aSigma*std::log( cTau*(time_n+dt)+1 ) + tau_i[0];
            real64 curTau_n = aSigma*std::log( cTau*(time_n)+1 ) + tau_i[0];
            tau.setValues< parallelHostPolicy >( curTau );
            tau_n.setValues< parallelHostPolicy >( curTau_n );

            // Call seismicity rate solver with the modified shear stress history
            propagator->integralSolverStep( time_n, dt, subRegion );

            // retrive seismicity rate from solver
            arrayView1d< real64 const > const R = subRegion.getField< geos::fields::inducedSeismicity::seismicityRate >();

            // compute analytical solution to the hard coded shear stress history
            real64 K = std::exp( curTau/aSigma - tau_i[0]/aSigma );
            real64 denom = (cTau*(time_n+dt - t_a) + 1)*std::exp((time_n+dt)/t_a ) + cTau*t_a;
            real64 Ranly = K/denom;

            // Check relative error against analytical solution
            checkRelativeError( R[0], Ranly, relTol );
          } );
        }
      } );
    } );
    time_n += dt;
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
