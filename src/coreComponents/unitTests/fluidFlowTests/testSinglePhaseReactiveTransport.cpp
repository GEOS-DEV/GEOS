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

#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseReactiveTransportFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseReactiveTransport.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "unitTests/fluidFlowTests/testSinglePhaseReactiveTransportUtils.hpp"

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::constitutive;
using namespace geos::constitutive::reactivefluid;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

// Sphinx start after input XML

char const * xmlInputSimple =
  R"xml(
  <Problem>
    <Solvers>
      <SinglePhaseReactiveTransport name="SinglePhaseReactiveFlow"
                                    logLevel="1"
                                    discretization="fluidTPFA"
                                    targetRegions="{ region }">
        <NonlinearSolverParameters newtonTol="1.0e-6"
                                   newtonMaxIter="100" />
        <LinearSolverParameters solverType="gmres"
                                krylovTol="1.0e-10" />
      </SinglePhaseReactiveTransport>
    </Solvers>
    <Mesh>
      <InternalMesh name="mesh"
                    elementTypes="{ C3D8 }"
                    xCoords="{ 0, 3 }"
                    yCoords="{ 0, 1 }"
                    zCoords="{ 0, 1 }"
                    nx="{ 3 }"
                    ny="{ 1 }"
                    nz="{ 1 }"
                    cellBlockNames="{ cb }" />
    </Mesh>
    <Geometry>
      <Box name="source"
           xMin="{ -0.01, -0.01, -0.01 }"
           xMax="{ 1.01, 1.01, 1.01 }"/>
      <Box name="sink"
           xMin="{ 1.99, -0.01, -0.01 }"
           xMax="{ 3.01, 1.01, 1.01 }"/>
    </Geometry>
    <NumericalMethods>
      <FiniteVolume>
        <TwoPointFluxApproximation name="fluidTPFA" />
      </FiniteVolume>
    </NumericalMethods>
    <ElementRegions>
      <CellElementRegion name="region"
                         cellBlocks="{ * }"
                         materialList="{ water, rock }" />
    </ElementRegions>
    <Constitutive>
      <CompressibleSolidConstantPermeability name="rock"
                                             solidModelName="nullSolid"
                                             porosityModelName="rockPorosity"
                                             permeabilityModelName="rockPerm" />
      <NullModel name="nullSolid" />
      <PressurePorosity name="rockPorosity"
                        defaultReferencePorosity="0.1"
                        referencePressure="0.0"
                        compressibility="0.0" />
      <ConstantPermeability name="rockPerm"
                            permeabilityComponents="{ 2.0e-10, 2.0e-10, 2.0e-10 }" />
      <ReactiveCompressibleSinglePhaseFluid name="water"
                                            defaultDensity="1000"
                                            defaultViscosity="0.001"
                                            referencePressure="0.0"
                                            compressibility="0.0"
                                            viscosibility="0.0"
                                            chemicalSystemType="simple"/>
    </Constitutive>
    <FieldSpecifications>
      <FieldSpecification name="Porosity"
                          initialCondition="1"
                          setNames="{ all }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="rockPorosity_referencePorosity"
                          scale="0.1"/>
      <FieldSpecification name="initialPressure"
                          initialCondition="1"
                          setNames="{ all }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="pressure"
                          scale="1e6" />
      <FieldSpecification name="sinkPressure"
                          setNames="{ sink }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="pressure"
                          scale="1e6" />
      <FieldSpecification name="sourcePressure"
                          setNames="{ source }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="pressure"
                          scale="1.2e6" />
      <FieldSpecification
        name="initialAggregate_0"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="0"
        scale="0.5"/>

      <FieldSpecification
        name="initialAggregate_1"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="1"
        scale="1"/>

      <FieldSpecification
        name="initialAggregate_2"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="2"
        scale="1e-16"/>

      <FieldSpecification
        name="sourceLog_0"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="0"
        scale="-1.386294"/>

      <FieldSpecification
        name="sourceLog_1"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="1"
        scale="0.0"/>
      
      <FieldSpecification
        name="sourceLog_2"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="2"
        scale="-36.841361"/>

      <FieldSpecification
        name="sinkLog_0"
        setNames="{ sink }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="0"
        scale="-1.386294"/>

      <FieldSpecification
        name="sinkLog_1"
        setNames="{ sink }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="1"
        scale="0.0"/>
      
      <FieldSpecification
        name="sinkLog_2"
        setNames="{ sink }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="2"
        scale="-36.841361"/>
    </FieldSpecifications>
  </Problem>
  )xml";

char const * xmlInputCarbonate =
  R"xml(
  <Problem>
    <Solvers>
      <SinglePhaseReactiveTransport name="SinglePhaseReactiveFlow"
                                    logLevel="1"
                                    discretization="fluidTPFA"
                                    targetRegions="{ region }">
        <NonlinearSolverParameters newtonTol="1.0e-6"
                                   newtonMaxIter="100" />
        <LinearSolverParameters directParallel="0"/>
      </SinglePhaseReactiveTransport>
    </Solvers>
    <Mesh>
      <InternalMesh name="mesh"
                    elementTypes="{ C3D8 }"
                    xCoords="{ 0, 3 }"
                    yCoords="{ 0, 1 }"
                    zCoords="{ 0, 1 }"
                    nx="{ 3 }"
                    ny="{ 1 }"
                    nz="{ 1 }"
                    cellBlockNames="{ cb }" />
    </Mesh>
    <Geometry>
      <Box name="source"
           xMin="{ -0.01, -0.01, -0.01 }"
           xMax="{ 1.01, 1.01, 1.01 }"/>
      <Box name="sink"
           xMin="{ 1.99, -0.01, -0.01 }"
           xMax="{ 3.01, 1.01, 1.01 }"/>
    </Geometry>
    <NumericalMethods>
      <FiniteVolume>
        <TwoPointFluxApproximation name="fluidTPFA" />
      </FiniteVolume>
    </NumericalMethods>
    <ElementRegions>
      <CellElementRegion name="region"
                         cellBlocks="{ * }"
                         materialList="{ water, rock }" />
    </ElementRegions>
    <Constitutive>
      <CompressibleSolidConstantPermeability name="rock"
                                             solidModelName="nullSolid"
                                             porosityModelName="rockPorosity"
                                             permeabilityModelName="rockPerm" />
      <NullModel name="nullSolid" />
      <PressurePorosity name="rockPorosity"
                        defaultReferencePorosity="0.1"
                        referencePressure="0.0"
                        compressibility="0.0" />
      <ConstantPermeability name="rockPerm"
                            permeabilityComponents="{ 2.0e-13, 2.0e-13, 2.0e-13 }" />
      <ReactiveCompressibleSinglePhaseFluid name="water"
                                            defaultDensity="1000"
                                            defaultViscosity="0.001"
                                            referencePressure="0.0"
                                            compressibility="5e-10"
                                            viscosibility="0.0"
                                            chemicalSystemType="carbonate"/>
    </Constitutive>
    <FieldSpecifications>
      <FieldSpecification name="Porosity"
                          initialCondition="1"
                          setNames="{ all }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="rockPorosity_referencePorosity"
                          scale="0.1"/>
      <FieldSpecification name="initialPressure"
                          initialCondition="1"
                          setNames="{ all }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="pressure"
                          scale="1e6" />
      <FieldSpecification name="sinkPressure"
                          setNames="{ sink }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="pressure"
                          scale="1e6" />
      <FieldSpecification name="sourcePressure"
                          setNames="{ source }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="pressure"
                          scale="1.2e6" />

      <FieldSpecification
        name="initialAggregate_CaCO3"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="0"
        scale="3.76e-3"/>

      <FieldSpecification
        name="initialAggregate_H"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="1"
        scale="3.76e-1"/>

      <FieldSpecification
        name="initialAggregate_HCO3"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="2"
        scale="3.76e-1"/>

      <FieldSpecification
        name="initialAggregate_Ca"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="3"
        scale="3.87e-2"/>

      <FieldSpecification
        name="initialAggregate_SO4"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="4"
        scale="3.21e-2"/>

      <FieldSpecification
        name="initialAggregate_Cl"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="5"
        scale="1.89"/>

      <FieldSpecification
        name="initialAggregate_Mg"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="6"
        scale="1.65e-2"/>
      
      <FieldSpecification
        name="initialAggregate_Na"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="7"
        scale="1.09"/>

      <FieldSpecification
        name="sourceLog_CaCO3"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="0"
        scale="-5.5833363216"/>

      <FieldSpecification
        name="sourceLog_H"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="1"
        scale="-7.7294305998"/>

      <FieldSpecification
        name="sourceLog_HCO3"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="2"
        scale="-7.8958055517"/>
      
      <FieldSpecification
        name="sourceLog_Ca"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="3"
        scale="-4.2187815037"/>

      <FieldSpecification
        name="sourceLog_SO4"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="4"
        scale="-5.9949216104"/>

      <FieldSpecification
        name="sourceLog_Cl"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="5"
        scale="0.61982840898"/>

      <FieldSpecification
        name="sourceLog_Mg"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="6"
        scale="-4.6170530778"/>

      <FieldSpecification
        name="sourceLog_Na"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="7"
        scale="0.069813174359"/>

      <FieldSpecification
        name="sinkLog_CaCO3"
        setNames="{ sink }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="0"
        scale="-5.5833363216"/>

      <FieldSpecification
        name="sinkLog_H"
        setNames="{ sink }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="1"
        scale="-7.7294305998"/>

      <FieldSpecification
        name="sinkLog_HCO3"
        setNames="{ sink }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="2"
        scale="-7.8958055517"/>
      
      <FieldSpecification
        name="sinkLog_Ca"
        setNames="{ sink }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="3"
        scale="-4.2187815037"/>

      <FieldSpecification
        name="sinkLog_SO4"
        setNames="{ sink }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="4"
        scale="-5.9949216104"/>

      <FieldSpecification
        name="sinkLog_Cl"
        setNames="{ sink }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="5"
        scale="0.61982840898"/>

      <FieldSpecification
        name="sinkLog_Mg"
        setNames="{ sink }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="6"
        scale="-4.6170530778"/>

      <FieldSpecification
        name="sinkLog_Na"
        setNames="{ sink }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="7"
        scale="0.069813174359"/>
    </FieldSpecifications>
  </Problem>
  )xml";
// Sphinx end before input XML

template< typename LAMBDA >
void testNumericalJacobian( SinglePhaseReactiveTransport & solver,
                            DomainPartition & domain,
                            real64 const perturbParameter,
                            real64 const relTol,
                            LAMBDA assembleFunction )
{
  CRSMatrix< real64, globalIndex > const & jacobian = solver.getLocalMatrix();
  array1d< real64 > residual( jacobian.numRows() );

  // assemble the analytical residual
  solver.resetStateToBeginningOfStep( domain );

  residual.zero();
  jacobian.zero();

  assembleFunction( jacobian.toViewConstSizes(), residual.toView() );
  residual.move( hostMemorySpace, false );

  // copy the analytical residual
  array1d< real64 > residualOrig( residual );

  // create the numerical jacobian
  jacobian.move( hostMemorySpace );
  CRSMatrix< real64, globalIndex > jacobianFD( jacobian );
  jacobianFD.zero();

  // fill jacobian FD
  fillCellCenteredNumericalJacobian( solver,
                                     domain,
                                     false,
                                     perturbParameter,
                                     residual.toView(),
                                     residualOrig.toView(),
                                     jacobian.toView(),
                                     jacobianFD.toView(),
                                     assembleFunction );

  // assemble the analytical jacobian
  solver.resetStateToBeginningOfStep( domain );

  residual.zero();
  jacobian.zero();
  assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

  compareLocalMatrices( jacobian.toViewConst(), jacobianFD.toViewConst(), relTol, 1e-6 );
}

class SinglePhaseReactiveTransportTest : public ::testing::Test
{
public:

SinglePhaseReactiveTransportTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {

    setupProblemFromXML( state.getProblemManager(), xmlInputCarbonate );
    solver = &state.getProblemManager().getPhysicsSolverManager().getGroup< SinglePhaseReactiveTransport >( "SinglePhaseReactiveFlow" );

    DomainPartition & domain = state.getProblemManager().getDomainPartition();

    solver->setupSystem( domain,
                         solver->getDofManager(),
                         solver->getLocalMatrix(),
                         solver->getSystemRhs(),
                         solver->getSystemSolution() );

    solver->implicitStepSetup( time, dt, domain );
  }

  static real64 constexpr time = 0.0;
  static real64 constexpr dt = 1000;
  static real64 constexpr eps = std::numeric_limits< real64 >::epsilon();

  GeosxState state;
  SinglePhaseReactiveTransport * solver;
};

real64 constexpr SinglePhaseReactiveTransportTest::time;
real64 constexpr SinglePhaseReactiveTransportTest::dt;
real64 constexpr SinglePhaseReactiveTransportTest::eps;

TEST_F( SinglePhaseReactiveTransportTest, jacobianNumericalCheck_flux )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-10; // 1% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    // The first input parameter denotes t_n, which is unused. Just input something here.
    solver->assembleFluxTerms( dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}

TEST_F( SinglePhaseReactiveTransportTest, jacobianNumericalCheck_accumulationBalance )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-6; // 1% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->assembleAccumulationTermsInMassBalanceAndSpeciesAmountEqs( dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
