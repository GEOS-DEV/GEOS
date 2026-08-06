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
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/PhysicsSolverBase.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseReactiveTransportFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseReactiveTransport.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "integrationTests/fluidFlowTests/testSinglePhaseReactiveTransportUtils.hpp"

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::constitutive;
using namespace geos::constitutive::reactivefluid;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

// Sphinx start after input XML

char const * xmlInputCarbonate =
  R"xml(
  <Problem>
    <Solvers>
      <SinglePhaseReactiveTransport name="SinglePhaseReactiveFlow"
                                    logLevel="1"
                                    discretization="fluidTPFA"
                                    targetRegions="{ region }">
        <NonlinearSolverParameters newtonTol="1.0e-4"
                                   newtonMaxIter="100"
                                   maxAllowedResidualNorm="1e100" />
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
    </Geometry>
    <NumericalMethods>
      <FiniteVolume>
        <TwoPointFluxApproximation name="fluidTPFA" />
      </FiniteVolume>
    </NumericalMethods>
    <ElementRegions>
      <CellElementRegion name="region"
                         cellBlocks="{ * }"
                         materialList="{ water, rock, diffusion }" />
    </ElementRegions>
    <Constitutive>
      <ReactiveSolidConstantPermeability name="rock"
                                         solidModelName="nullSolid"
                                         porosityModelName="rockPorosity"
                                         permeabilityModelName="rockPerm"
                                         surfaceAreaModelName="rockSurfaceArea" />
      <NullModel name="nullSolid" />
      <ReactivePorosity name="rockPorosity"
                        defaultReferencePorosity="0.1"
                        defaultInitialVolumeFractions="{ 0.1 }"
                        molarWeights="{ 0.100087 }"
                        mineralDensities="{ 2710 }"/>
      <ConstantPermeability name="rockPerm"
                            permeabilityComponents="{ 3.0e-16, 3.0e-16, 3.0e-16 }" />
      <ConstantSurfaceArea name="rockSurfaceArea"
                           defaultInitialSurfaceArea="{ 2.4573e+06 }"/>
      <ReactiveCompressibleSinglePhaseFluid name="water"
                                            defaultDensity="1000"
                                            defaultViscosity="0.001"
                                            referencePressure="0.0"
                                            compressibility="5e-10"
                                            viscosibility="0.0"
                                            chemicalSystemType="carbonate"/>
      <ConstantDiffusion name="diffusion"
                         phaseNames="{ water }"
                         defaultPhaseDiffusivityMultipliers="{ 1 }"
                         diffusivityComponents="{ 1e-2, 1e-2, 1e-2 }"/>
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
                          setNames="{ xpos }"
                          objectPath="faceManager"
                          fieldName="pressure"
                          scale="1.0e6"/>
      <FieldSpecification name="sourcePressure"
                          setNames="{ source }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="pressure"
                          scale="2.0e6" />

      <FieldSpecification name="surfaceArea"
                          initialCondition="1"
                          setNames="{ all }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="rockSurfaceArea_initialSurfaceArea"
                          component="0"
                          scale="2.4573e+06"/>

      <FieldSpecification name="surfaceAreaSource"
                          initialCondition="1"
                          setNames="{ source }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="rockSurfaceArea_initialSurfaceArea"
                          component="0"
                          scale="0.0"/>

      <FieldSpecification
        name="initialAggregate_H_all"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="0"
        scale="1.585e-7"/>

      <FieldSpecification
        name="initialAggregate_HCO3_all"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="1"
        scale="8.293e-3"/>

      <FieldSpecification
        name="initialAggregate_Ca_all"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="2"
        scale="2.171e-3"/>

      <FieldSpecification
        name="initialAggregate_SO4_all"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="3"
        scale="1.666e-6"/>

      <FieldSpecification
        name="initialAggregate_Cl_all"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="4"
        scale="2.821e-5"/>

      <FieldSpecification
        name="initialAggregate_Mg_all"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="5"
        scale="1.605e-3"/>
      
      <FieldSpecification
        name="initialAggregate_Na_all"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="6"
        scale="7.817e-4"/>
      
      <FieldSpecification
        name="initialAggregate_H_source"
        initialCondition="1"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="0"
        scale="1.585e-7"/>

      <FieldSpecification
        name="initialAggregate_HCO3_source"
        initialCondition="1"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="1"
        scale="7.317e-2"/>

      <FieldSpecification
        name="initialAggregate_Ca_source"
        initialCondition="1"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="2"
        scale="1.517e-2"/>

      <FieldSpecification
        name="initialAggregate_SO4_source"
        initialCondition="1"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="3"
        scale="1.666e-6"/>

      <FieldSpecification
        name="initialAggregate_Cl_source"
        initialCondition="1"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="4"
        scale="2.821e-5"/>

      <FieldSpecification
        name="initialAggregate_Mg_source"
        initialCondition="1"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="5"
        scale="1.605e-3"/>
      
      <FieldSpecification
        name="initialAggregate_Na_source"
        initialCondition="1"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_primarySpeciesAggregateConcentration"
        component="6"
        scale="7.817e-4"/>

      <FieldSpecification
        name="initialPrimary_H_all"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_initialPrimarySpeciesConcentration"
        component="0"
        scale="1.585e-7"/>

      <FieldSpecification
        name="initialPrimary_HCO3_all"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_initialPrimarySpeciesConcentration"
        component="1"
        scale="8.293e-3"/>

      <FieldSpecification
        name="initialPrimary_Ca_all"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_initialPrimarySpeciesConcentration"
        component="2"
        scale="2.171e-3"/>

      <FieldSpecification
        name="initialPrimary_SO4_all"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_initialPrimarySpeciesConcentration"
        component="3"
        scale="1.666e-6"/>

      <FieldSpecification
        name="initialPrimary_Cl_all"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_initialPrimarySpeciesConcentration"
        component="4"
        scale="2.821e-5"/>

      <FieldSpecification
        name="initialPrimary_Mg_all"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_initialPrimarySpeciesConcentration"
        component="5"
        scale="1.605e-3"/>
      
      <FieldSpecification
        name="initialPrimary_Na_all"
        initialCondition="1"
        setNames="{ all }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_initialPrimarySpeciesConcentration"
        component="6"
        scale="7.817e-4"/>

      <FieldSpecification
        name="initialPrimary_H_source"
        initialCondition="1"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_initialPrimarySpeciesConcentration"
        component="0"
        scale="1.585e-7"/>

      <FieldSpecification
        name="initialPrimary_HCO3_source"
        initialCondition="1"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_initialPrimarySpeciesConcentration"
        component="1"
        scale="7.317e-2"/>

      <FieldSpecification
        name="initialPrimary_Ca_source"
        initialCondition="1"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_initialPrimarySpeciesConcentration"
        component="2"
        scale="1.517e-2"/>

      <FieldSpecification
        name="initialPrimary_SO4_source"
        initialCondition="1"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_initialPrimarySpeciesConcentration"
        component="3"
        scale="1.666e-6"/>

      <FieldSpecification
        name="initialPrimary_Cl_source"
        initialCondition="1"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_initialPrimarySpeciesConcentration"
        component="4"
        scale="2.821e-5"/>

      <FieldSpecification
        name="initialPrimary_Mg_source"
        initialCondition="1"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_initialPrimarySpeciesConcentration"
        component="5"
        scale="1.605e-3"/>
      
      <FieldSpecification
        name="initialPrimary_Na_source"
        initialCondition="1"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="water_initialPrimarySpeciesConcentration"
        component="6"
        scale="7.817e-4"/>
      
      <FieldSpecification
        name="log_H_source"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="0"
        scale="-19.21448423"/>

      <FieldSpecification
        name="log_HCO3_source"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="1"
        scale="-2.73216432"/>

      <FieldSpecification
        name="log_Ca_source"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="2"
        scale="-4.77653747"/>

      <FieldSpecification
        name="log_SO4_source"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="3"
        scale="-14.28117504"/>

      <FieldSpecification
        name="log_Cl_source"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="4"
        scale="-10.47763277"/>

      <FieldSpecification
        name="log_Mg_source"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="5"
        scale="-6.43480025"/>
      
      <FieldSpecification
        name="log_Na_source"
        setNames="{ source }"
        objectPath="ElementRegions/region/cb"
        fieldName="logPrimarySpeciesConcentration"
        component="6"
        scale="-7.15404368"/>
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
  static real64 constexpr dt = 0.01;
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
