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
 * ------------------------------------------------------------------------------------------------------------
 */

#include <gtest/gtest.h>
#include <tuple>

#include "integrationTests/fluidFlowTests/testCompFlowUtils.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/DomainPartition.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"

#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVM.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

// Pressure and saturation L2 error tolerances
static constexpr real64 PRESSURE_L2_TOLERANCE = 1.0e-8;
static constexpr real64 SATURATION_L2_TOLERANCE = 1.0e-8;

// Single time step
static constexpr real64 TIME_STEP = 1.0e-2; // 0.01

static constexpr auto INNER_TPFA = "TPFA";
static constexpr auto INNER_quasiTPFA = "quasiTPFA";
static constexpr auto INNER_BDVLM = "beiraoDaVeigaLipnikovManzini";

// Generate XML for compositional multiphase FV-TPFA (simplified from the user-provided template)
static std::string generateXmlInputCompTPFA( std::string const & meshFile )
{
  std::ostringstream oss;
  oss <<
    R"xml(
<?xml version="1.0"?>
<Problem>
  <Solvers gravityVector="{ 0.0, 0.0, 0.0 }">
    <CompositionalMultiphaseFVM
      name="FlowSolverTPFA"
      discretization="TPFA"
      targetRegions="{ Domain }"
      logLevel="0"
      useMass="1"
      temperature="300">
      <NonlinearSolverParameters
        newtonTol="1e-9"
        lineSearchAction="None"
        maxTimeStepCuts="100"
        maxSubSteps="100"
        newtonMaxIter="10"/>
      <LinearSolverParameters
        directParallel="0"
        reuseFactorization="0"/>
    </CompositionalMultiphaseFVM>
  </Solvers>

  <Mesh>
    <VTKMesh
      name="mesh1"
      logLevel="0"
      file=")xml" << meshFile <<
    R"xml("/>
  </Mesh>

  <Geometry>
    <Box name="westBC" xMin="{ -0.001, 0.0, 0.0}" xMax="{ +0.001, 1.0, 1.0}"/>
    <Box name="eastBC" xMin="{ +0.999, 0.0, 0.0}" xMax="{ +1.001, 1.0, 1.0}"/>
  </Geometry>

  <NumericalMethods>
    <FiniteVolume>
      <TwoPointFluxApproximation
        name="TPFA"
        upwindingScheme="HU2PH"/>
    </FiniteVolume>
  </NumericalMethods>

  <ElementRegions>
    <CellElementRegion name="Domain" cellBlocks="{ * }" materialList="{ fluid, rock, relperm }"/>
  </ElementRegions>

  <Constitutive>
    <InvariantImmiscibleFluid name="fluid"
      componentNames="{ H2O, CH4 }"
      phaseNames="{ water, gas }"
      densities="{ 1.0, 1.0 }"
      componentMolarWeight="{ 1.0, 1.0 }"
      viscosities="{ 1.0, 1.0 }"/>

    <CompressibleSolidConstantPermeability
      name="rock"
      solidModelName="nullSolid"
      porosityModelName="rockPorosity"
      permeabilityModelName="rockPerm"/>

    <NullModel name="nullSolid"/>

    <PressurePorosity name="rockPorosity" defaultReferencePorosity="0.1" referencePressure="0" compressibility="0.0"/>

    <ConstantPermeability name="rockPerm" permeabilityComponents="{ 1.0, 1.0, 1.0 }"/>

    <TableRelativePermeability
      name="relperm"
      phaseNames="{ water, gas }"
      wettingNonWettingRelPermTableNames="{ waterRelativePermeabilityTable, gasRelativePermeabilityTable }"/>
  </Constitutive>

  <Functions>
    <!-- Linear relative permeability tables -->
    <TableFunction
      name="waterRelativePermeabilityTable"
      coordinates="{ 0.0, 0.5, 1.0 }"
      values="{ 0.0, 0.5, 1.0 }"/>
    <TableFunction
      name="gasRelativePermeabilityTable"
      coordinates="{ 0.0, 0.5, 1.0 }"
      values="{ 0.0, 0.5, 1.0 }"/>
  </Functions>

  <FieldSpecifications>
    <FieldSpecification name="Porosity" initialCondition="1" setNames="{ all }"
      objectPath="ElementRegions/Domain" fieldName="rockPorosity_referencePorosity" scale="0.1"/>

    <FieldSpecification name="initialPressure" initialCondition="1" setNames="{ all }"
      objectPath="ElementRegions/Domain" fieldName="pressure" scale="1.0"/>

    <!-- Uniform initial global component fractions to avoid external files -->
    <FieldSpecification name="initZH2O" initialCondition="1" setNames="{ all }"
      objectPath="ElementRegions/Domain" fieldName="globalCompFraction" component="0" scale="1.0"/>

    <FieldSpecification name="initZCH4" initialCondition="1" setNames="{ all }"
      objectPath="ElementRegions/Domain" fieldName="globalCompFraction" component="1" scale="0.0"/>

    <!-- Pressure and composition BCs on faces -->
    <FieldSpecification name="westPressure" objectPath="faceManager"
      fieldName="pressure" scale="2.0" setNames="{ westBC }"/>

    <FieldSpecification name="westZH2O" objectPath="faceManager"
      fieldName="globalCompFraction" component="0" scale="0.0" setNames="{ westBC }"/>

    <FieldSpecification name="westZCH4" objectPath="faceManager"
      fieldName="globalCompFraction" component="1" scale="1.0" setNames="{ westBC }"/>

    <FieldSpecification name="westT" objectPath="faceManager"
      fieldName="temperature" scale="300.0" setNames="{ westBC }"/>

    <FieldSpecification name="eastPressure" objectPath="faceManager"
      fieldName="pressure" scale="1.0" setNames="{ eastBC }"/>

    <FieldSpecification name="eastZH2O" objectPath="faceManager"
      fieldName="globalCompFraction" component="0" scale="0.0" setNames="{ eastBC }"/>

    <FieldSpecification name="eastZCH4" objectPath="faceManager"
      fieldName="globalCompFraction" component="1" scale="1.0" setNames="{ eastBC }"/>

    <FieldSpecification name="eastT" objectPath="faceManager"
      fieldName="temperature" scale="300.0" setNames="{ eastBC }"/>
  </FieldSpecifications>
</Problem>
)xml";
  return oss.str();
}

// Generate XML for compositional multiphase Hybrid MFD with configurable inner product (default TPFA)
static std::string generateXmlInputCompMFD( std::string const & innerProductType,
                                            std::string const & meshFile )
{
  std::ostringstream oss;
  oss <<
    R"xml(
<?xml version="1.0"?>
<Problem>
  <Solvers gravityVector="{ 0.0, 0.0, 0.0 }">
    <CompositionalMultiphaseHybridFVM
      name="FlowSolverMFD"
      discretization="compositionalHybridMimetic"
      targetRegions="{ Domain }"
      logLevel="0"
      useMass="1"
      temperature="300">
      <NonlinearSolverParameters
        newtonTol="1e-9"
        lineSearchAction="None"
        maxTimeStepCuts="100"
        maxSubSteps="100"
        newtonMaxIter="10"/>
      <LinearSolverParameters
        directParallel="0"
        reuseFactorization="0"/>
    </CompositionalMultiphaseHybridFVM>
  </Solvers>

  <Mesh>
    <VTKMesh
      name="mesh1"
      logLevel="0"
      file=")xml" << meshFile <<
    R"xml("/>
  </Mesh>

  <Geometry>
    <Box name="westBC" xMin="{ -0.001, 0.0, 0.0}" xMax="{ +0.001, 1.0, 1.0}"/>
    <Box name="eastBC" xMin="{ +0.999, 0.0, 0.0}" xMax="{ +1.001, 1.0, 1.0}"/>
  </Geometry>

  <NumericalMethods>
    <FiniteVolume>
      <HybridMimeticDiscretization name="compositionalHybridMimetic" innerProductType=")xml"
      << innerProductType <<
    R"xml("/>
    </FiniteVolume>
  </NumericalMethods>

  <ElementRegions>
    <CellElementRegion name="Domain" cellBlocks="{ * }" materialList="{ fluid, rock, relperm }"/>
  </ElementRegions>

  <Functions>
    <!-- Linear relative permeability tables -->
    <TableFunction
      name="waterRelativePermeabilityTable"
      coordinates="{ 0.0, 0.5, 1.0 }"
      values="{ 0.0, 0.5, 1.0 }"/>
    <TableFunction
      name="gasRelativePermeabilityTable"
      coordinates="{ 0.0, 0.5, 1.0 }"
      values="{ 0.0, 0.5, 1.0 }"/>
  </Functions>

  <Constitutive>
    <InvariantImmiscibleFluid name="fluid"
      componentNames="{ H2O, CH4 }"
      phaseNames="{ water, gas }"
      densities="{ 1.0, 1.0 }"
      componentMolarWeight="{ 1.0, 1.0 }"
      viscosities="{ 1.0, 1.0 }"/>

    <CompressibleSolidConstantPermeability
      name="rock"
      solidModelName="nullSolid"
      porosityModelName="rockPorosity"
      permeabilityModelName="rockPerm"/>

    <NullModel name="nullSolid"/>

    <PressurePorosity name="rockPorosity" defaultReferencePorosity="0.1" referencePressure="0" compressibility="0.0"/>

    <ConstantPermeability name="rockPerm" permeabilityComponents="{ 1.0, 1.0, 1.0 }"/>

    <TableRelativePermeability
      name="relperm"
      phaseNames="{ water, gas }"
      wettingNonWettingRelPermTableNames="{ waterRelativePermeabilityTable, gasRelativePermeabilityTable }"/>
  </Constitutive>

  <FieldSpecifications>
    <FieldSpecification name="Porosity" initialCondition="1" setNames="{ all }"
      objectPath="ElementRegions/Domain" fieldName="rockPorosity_referencePorosity" scale="0.1"/>

    <FieldSpecification name="initialPressure" initialCondition="1" setNames="{ all }"
      objectPath="ElementRegions/Domain" fieldName="pressure" scale="1.0"/>
    <FieldSpecification name="initialfacePressure" initialCondition="1" setNames="{ all }"
      objectPath="faceManager" fieldName="facePressure" scale="1.0"/>

    <!-- Uniform initial global component fractions to avoid external files -->
    <FieldSpecification name="initZH2O" initialCondition="1" setNames="{ all }"
      objectPath="ElementRegions/Domain" fieldName="globalCompFraction" component="0" scale="1.0"/>

    <FieldSpecification name="initZCH4" initialCondition="1" setNames="{ all }"
      objectPath="ElementRegions/Domain" fieldName="globalCompFraction" component="1" scale="0.0"/>

    <!-- Face-based initial component fractions for hybrid solver -->
    <FieldSpecification name="initFaceZH2O" initialCondition="1" setNames="{ all }"
      objectPath="faceManager" fieldName="faceGlobalCompFraction" component="0" scale="1.0"/>

    <FieldSpecification name="initFaceZCH4" initialCondition="1" setNames="{ all }"
      objectPath="faceManager" fieldName="faceGlobalCompFraction" component="1" scale="0.0"/>

    <!-- Face-based initial temperature for hybrid solver -->
    <FieldSpecification name="initFaceTemperature" initialCondition="1" setNames="{ all }"
      objectPath="faceManager" fieldName="faceTemperature" scale="300.0"/>

    <!-- Pressure BC regions for hybrid solver use bc* fields -->
    <FieldSpecification name="westPressure" objectPath="faceManager"
      fieldName="bcPressure" scale="2.0" setNames="{ westBC }"/>

    <FieldSpecification name="westZH2O" objectPath="faceManager"
      fieldName="bcGlobalCompFraction" component="0" scale="0.0" setNames="{ westBC }"/>

    <FieldSpecification name="westZCH4" objectPath="faceManager"
      fieldName="bcGlobalCompFraction" component="1" scale="1.0" setNames="{ westBC }"/>

    <FieldSpecification name="westT" objectPath="faceManager"
      fieldName="bcTemperature" scale="300.0" setNames="{ westBC }"/>

    <FieldSpecification name="eastPressure" objectPath="faceManager"
      fieldName="bcPressure" scale="1.0" setNames="{ eastBC }"/>

    <FieldSpecification name="eastZH2O" objectPath="faceManager"
      fieldName="bcGlobalCompFraction" component="0" scale="0.0" setNames="{ eastBC }"/>

    <FieldSpecification name="eastZCH4" objectPath="faceManager"
      fieldName="bcGlobalCompFraction" component="1" scale="1.0" setNames="{ eastBC }"/>

    <FieldSpecification name="eastT" objectPath="faceManager"
      fieldName="bcTemperature" scale="300.0" setNames="{ eastBC }"/>
  </FieldSpecifications>
</Problem>
)xml";
  return oss.str();
}

// Helper: copy arrayView1d<real64 const> to std::vector<real64>
static inline std::vector< real64 > arrayViewToVector( arrayView1d< real64 const > & arr )
{
  arr.move( hostMemorySpace, false );
  return std::vector< real64 >( arr.data(), arr.data() + arr.size() );
}

// Helper: copy a 2D view to flat vector (row-major by phase index), accepts any 2D view type
template< typename View2D >
static inline std::vector< real64 > arrayView2dToVector( View2D & arr )
{
  arr.move( hostMemorySpace, false );
  localIndex nCells = arr.size( 0 );
  localIndex nCols = arr.size( 1 );
  std::vector< real64 > out;
  out.reserve( static_cast< size_t >( nCells ) * static_cast< size_t >( nCols ) );
  for( localIndex i = 0; i < nCells; ++i )
  {
    for( localIndex j = 0; j < nCols; ++j )
    {
      out.push_back( static_cast< real64 >( arr( i, j ) ) );
    }
  }
  return out;
}

// Helper: compute normalized L2 error for linear pressure profile p(x) = 2*(1-x) + x = 2 - x
static inline real64 computeNormalizedPressureL2Error( CellElementSubRegion & subRegion )
{
  arrayView2d< real64 const > centers = subRegion.getElementCenter();
  arrayView1d< real64 const > volumes = subRegion.getElementVolume();
  arrayView1d< real64 const > const p_h = subRegion.getField< fields::flow::pressure >();

  RAJA::ReduceSum< parallelDeviceReduce, real64 > squaredPressureError_ReduceSum( 0.0 );
  RAJA::ReduceSum< parallelDeviceReduce, real64 > squaredExactPressure_ReduceSum( 0.0 );
  localIndex const n_cells = subRegion.size();

  forAll< geos::parallelDevicePolicy<> >( n_cells, [=] GEOS_HOST_DEVICE ( localIndex const i )
  {
    real64 const x = centers[i][0];
    real64 const cellVolume = volumes[i];
    real64 const pNumeric = p_h[i];
    real64 const pExact = 2.0 * (1.0 - x) + 1.0 * x;
    squaredPressureError_ReduceSum += (pNumeric - pExact) * (pNumeric - pExact) * cellVolume;
    squaredExactPressure_ReduceSum += pExact * pExact * cellVolume;
  } );

  // Gather global reductions
  real64 const squaredPressureError = squaredPressureError_ReduceSum.get();
  real64 const squaredExactPressure = squaredExactPressure_ReduceSum.get();

  // Compute pressure relative L2 error
  real64 const normalizedL2Error = std::sqrt( squaredPressureError ) / std::sqrt( squaredExactPressure );

  return normalizedL2Error;
}

// TPFA test: check linear pressure profile on regular polyhedral mesh
class CompositionalTPFAIntegrationTest : public ::testing::TestWithParam< const char * >
{
public:
  CompositionalTPFAIntegrationTest()
    : state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) ) {}

protected:
  void SetUp() override
  {
    testBinaryDir = TEST_BINARY_DIR;
    std::string meshFile = testBinaryDir + std::string( "/" ) + GetParam();
    std::string xmlInput = generateXmlInputCompTPFA( meshFile );
    setupProblemFromXML( state.getProblemManager(), xmlInput.c_str() );
  }

  GeosxState state;
  std::string testBinaryDir;
};

INSTANTIATE_TEST_SUITE_P(
  MeshFiles,
  CompositionalTPFAIntegrationTest,
  ::testing::Values(
    "polyhedral_voronoi_regular.vtu",
    "hex_pyr_tet_nested_mixed.vtu",
    "polyhedral_voronoi_lattice.vtu",
    "polyhedral_voronoi_complex.vtu"
    )
  );

TEST_P( CompositionalTPFAIntegrationTest, PressureFieldL2Error )
{
  ProblemManager & problemManager = state.getProblemManager();
  DomainPartition & domain = problemManager.getDomainPartition();

  auto & solver =
    dynamic_cast< CompositionalMultiphaseFVM & >(
      problemManager.getPhysicsSolverManager().getGroup< CompositionalMultiphaseFVM >( "FlowSolverTPFA" ) );

  // One implicit step with dt = 0.01
  solver.setupSystem( domain, solver.getDofManager(),
                      solver.getLocalMatrix(), solver.getSystemRhs(),
                      solver.getSystemSolution() );
  solver.implicitStepSetup( 0.0, TIME_STEP, domain );
  solver.solverStep( 0.0, TIME_STEP, 0, domain );
  solver.updateState( domain );
  solver.implicitStepComplete( 0.0, TIME_STEP, domain );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();
  CellElementSubRegion & subRegion = mesh.getElemManager().getRegion( 0 ).getSubRegion< CellElementSubRegion >( 0 );

  real64 l2Error = computeNormalizedPressureL2Error( subRegion );

  // Expect exact solution only on the k-orthogonal regular mesh; expect non-zero error on mixed mesh
  std::string meshFile = GetParam();
  if( meshFile == std::string( "polyhedral_voronoi_regular.vtu" ) )
  {
    EXPECT_NEAR( l2Error, 0.0, PRESSURE_L2_TOLERANCE );
  }
  else
  {
    EXPECT_GT( l2Error, PRESSURE_L2_TOLERANCE );
  }
}

// MFD-TPFA test: check linear pressure profile on polyhedral meshes
class CompositionalMFDTPFAIntegrationTest : public ::testing::TestWithParam< const char * >
{
public:
  CompositionalMFDTPFAIntegrationTest()
    : state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) ) {}

protected:
  void SetUp() override
  {
    testBinaryDir = TEST_BINARY_DIR;
    std::string meshFile = testBinaryDir + std::string( "/" ) + GetParam();
    std::string xmlInput = generateXmlInputCompMFD( INNER_TPFA, meshFile );
    setupProblemFromXML( state.getProblemManager(), xmlInput.c_str() );
  }

  GeosxState state;
  std::string testBinaryDir;
};

INSTANTIATE_TEST_SUITE_P(
  MeshFiles,
  CompositionalMFDTPFAIntegrationTest,
  ::testing::Values(
    "polyhedral_voronoi_regular.vtu",
    "hex_pyr_tet_nested_mixed.vtu",
    "polyhedral_voronoi_lattice.vtu",
    "polyhedral_voronoi_complex.vtu"
    )
  );

TEST_P( CompositionalMFDTPFAIntegrationTest, PressureFieldL2Error )
{
  ProblemManager & problemManager = state.getProblemManager();
  DomainPartition & domain = problemManager.getDomainPartition();

  auto & solver =
    dynamic_cast< CompositionalMultiphaseHybridFVM & >(
      problemManager.getPhysicsSolverManager().getGroup< CompositionalMultiphaseHybridFVM >( "FlowSolverMFD" ) );

  // One implicit step with dt = 0.01
  solver.setupSystem( domain, solver.getDofManager(),
                      solver.getLocalMatrix(), solver.getSystemRhs(),
                      solver.getSystemSolution() );
  solver.implicitStepSetup( 0.0, TIME_STEP, domain );
  solver.solverStep( 0.0, TIME_STEP, 0, domain );
  solver.updateState( domain );
  solver.implicitStepComplete( 0.0, TIME_STEP, domain );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();
  CellElementSubRegion & subRegion = mesh.getElemManager().getRegion( 0 ).getSubRegion< CellElementSubRegion >( 0 );

  real64 l2Error = computeNormalizedPressureL2Error( subRegion );

  // Expect exact solution only on the k-orthogonal regular mesh; expect non-zero error on mixed mesh
  std::string meshFile = GetParam();
  if( meshFile == std::string( "polyhedral_voronoi_regular.vtu" ) )
  {
    EXPECT_NEAR( l2Error, 0.0, PRESSURE_L2_TOLERANCE );
  }
  else
  {
    EXPECT_GT( l2Error, PRESSURE_L2_TOLERANCE );
  }
}

// Cross-check: TPFA vs MFD(TPFA): pressure and saturation profiles must match exactly
class CompositionalTPFAvsMFDTPFA : public ::testing::TestWithParam< const char * > {};

INSTANTIATE_TEST_SUITE_P(
  MeshFiles,
  CompositionalTPFAvsMFDTPFA,
  ::testing::Values(
    "polyhedral_voronoi_regular.vtu",
    "hex_pyr_tet_nested_mixed.vtu",
    "polyhedral_voronoi_lattice.vtu",
    "polyhedral_voronoi_complex.vtu"
    )
  );

TEST_P( CompositionalTPFAvsMFDTPFA, PressureAndSaturationComparison )
{
  std::string const testBinaryDir = TEST_BINARY_DIR;
  std::string const meshFile = testBinaryDir + std::string( "/" ) + GetParam();

  std::vector< real64 > p_tpfa, p_mfd;
  std::vector< real64 > sat_tpfa, sat_mfd;  // phase volume fractions flattened [cell,phase]
  localIndex n_cells_tpfa = 0, n_cells_mfd = 0;
  localIndex n_phases_tpfa = 2, n_phases_mfd = 2;  // From XML phaseNames

  // --- Run TPFA solver ---
  {
    GeosxState tpfaState( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
    std::string xml = generateXmlInputCompTPFA( meshFile );
    setupProblemFromXML( tpfaState.getProblemManager(), xml.c_str() );

    ProblemManager & pm = tpfaState.getProblemManager();
    DomainPartition & domain = pm.getDomainPartition();

    auto & solver =
      dynamic_cast< CompositionalMultiphaseFVM & >(
        pm.getPhysicsSolverManager().getGroup< CompositionalMultiphaseFVM >( "FlowSolverTPFA" ) );

    solver.setupSystem( domain, solver.getDofManager(),
                        solver.getLocalMatrix(), solver.getSystemRhs(),
                        solver.getSystemSolution() );
    solver.implicitStepSetup( 0.0, TIME_STEP, domain );
    solver.solverStep( 0.0, TIME_STEP, 0, domain );
    solver.updateState( domain );
    solver.implicitStepComplete( 0.0, TIME_STEP, domain );

    MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();
    CellElementSubRegion & subRegion = mesh.getElemManager().getRegion( 0 ).getSubRegion< CellElementSubRegion >( 0 );

    arrayView1d< real64 const > p_h = subRegion.getField< fields::flow::pressure >();
    p_tpfa = arrayViewToVector( p_h );

    auto sat = subRegion.getField< fields::flow::phaseVolumeFraction >();
    sat_tpfa = arrayView2dToVector( sat );
  }

  // --- Run MFD solver (innerProductType = TPFA) ---
  {
    GeosxState mfdState( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
    std::string xml = generateXmlInputCompMFD( INNER_TPFA, meshFile );
    setupProblemFromXML( mfdState.getProblemManager(), xml.c_str() );

    ProblemManager & pm = mfdState.getProblemManager();
    DomainPartition & domain = pm.getDomainPartition();

    auto & solver =
      dynamic_cast< CompositionalMultiphaseHybridFVM & >(
        pm.getPhysicsSolverManager().getGroup< CompositionalMultiphaseHybridFVM >( "FlowSolverMFD" ) );

    solver.setupSystem( domain, solver.getDofManager(),
                        solver.getLocalMatrix(), solver.getSystemRhs(),
                        solver.getSystemSolution() );
    solver.implicitStepSetup( 0.0, TIME_STEP, domain );
    solver.solverStep( 0.0, TIME_STEP, 0, domain );
    solver.updateState( domain );
    solver.implicitStepComplete( 0.0, TIME_STEP, domain );

    MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();
    CellElementSubRegion & subRegion = mesh.getElemManager().getRegion( 0 ).getSubRegion< CellElementSubRegion >( 0 );

    arrayView1d< real64 const > p_h = subRegion.getField< fields::flow::pressure >();
    p_mfd = arrayViewToVector( p_h );

    auto sat = subRegion.getField< fields::flow::phaseVolumeFraction >();
    sat_mfd = arrayView2dToVector( sat );
  }

  ASSERT_EQ( n_cells_tpfa, n_cells_mfd );
  ASSERT_EQ( n_phases_tpfa, n_phases_mfd );

  // Pressure profile equality
  for( localIndex i = 0; i < n_cells_tpfa; ++i )
  {
    real64 const diff = std::abs( p_tpfa[i] - p_mfd[i] );
    EXPECT_NEAR( diff, 0.0, PRESSURE_L2_TOLERANCE );
  }

  // Saturation profile (phase volume fraction) equality
  size_t const nTot = static_cast< size_t >( n_cells_tpfa ) * static_cast< size_t >( n_phases_tpfa );
  for( size_t k = 0; k < nTot; ++k )
  {
    real64 const diff = std::abs( sat_tpfa[k] - sat_mfd[k] );
    EXPECT_NEAR( diff, 0.0, SATURATION_L2_TOLERANCE );
  }
}

// MFD non-TPFA inner products: check linear pressure profile is exact on all meshes
class CompositionalMFDNonTPFAExactnessTest : public ::testing::TestWithParam< std::tuple< const char *, const char * > >
{
public:
  CompositionalMFDNonTPFAExactnessTest()
    : state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) ) {}

protected:
  void SetUp() override
  {
    testBinaryDir = TEST_BINARY_DIR;
    auto const params = GetParam();
    innerProduct = std::get< 1 >( params );
    std::string const meshRel = std::get< 0 >( params );
    std::string const meshFile = testBinaryDir + std::string( "/" ) + meshRel;
    std::string xmlInput = generateXmlInputCompMFD( innerProduct, meshFile );
    setupProblemFromXML( state.getProblemManager(), xmlInput.c_str() );
  }

  GeosxState state;
  std::string testBinaryDir;
  std::string innerProduct;
};

INSTANTIATE_TEST_SUITE_P(
  MeshAndInnerProducts,
  CompositionalMFDNonTPFAExactnessTest,
  ::testing::Combine(
    ::testing::Values(
      "polyhedral_voronoi_regular.vtu",
      "hex_pyr_tet_nested_mixed.vtu",
      "polyhedral_voronoi_lattice.vtu",
      "polyhedral_voronoi_complex.vtu"
      ),
    ::testing::Values(
      INNER_quasiTPFA,
      INNER_BDVLM
      )
    )
  );

TEST_P( CompositionalMFDNonTPFAExactnessTest, PressureFieldL2ErrorExact )
{
  ProblemManager & problemManager = state.getProblemManager();
  DomainPartition & domain = problemManager.getDomainPartition();

  auto & solver =
    dynamic_cast< CompositionalMultiphaseHybridFVM & >(
      problemManager.getPhysicsSolverManager().getGroup< CompositionalMultiphaseHybridFVM >( "FlowSolverMFD" ) );

  // One implicit step with dt = 0.01
  solver.setupSystem( domain, solver.getDofManager(),
                      solver.getLocalMatrix(), solver.getSystemRhs(),
                      solver.getSystemSolution() );
  solver.implicitStepSetup( 0.0, TIME_STEP, domain );
  solver.solverStep( 0.0, TIME_STEP, 0, domain );
  solver.updateState( domain );
  solver.implicitStepComplete( 0.0, TIME_STEP, domain );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();
  CellElementSubRegion & subRegion = mesh.getElemManager().getRegion( 0 ).getSubRegion< CellElementSubRegion >( 0 );

  real64 const normalizedL2Error = computeNormalizedPressureL2Error( subRegion );

  // Expect exact solution for these inner products on both meshes
  EXPECT_NEAR( normalizedL2Error, 0.0, PRESSURE_L2_TOLERANCE );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
