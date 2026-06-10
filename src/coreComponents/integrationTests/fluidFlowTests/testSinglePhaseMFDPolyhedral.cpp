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
#include "integrationTests/fluidFlowTests/testCompFlowUtils.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/DomainPartition.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVM.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseHybridFVM.hpp"


// This file implements integration tests for polyhedral discretizations of single-phase flow.
//
// Test summary:
// 1. Parameterized TPFA integration tests
// 2. Parameterized MFD integration tests with various inner products
// 3. Cross-check test ensuring that MFD with innerProductType="TPFA" reproduces
//    the same pressure field as the TPFA solver
//
// Tested Meshes:
// - polyhedral_voronoi_complex.vtu
// - polyhedral_voronoi_lattice.vtu
// - polyhedral_voronoi_regular.vtu
//
// Inner Products for MFD:
// - TPFA
// - QuasiTPFA
// - QuasiRT
// - Simple
// - BdVLM
//
// L2 error comparisons are performed to ensure:
// - TPFA produces machine-precision exact solutions on regular meshes
// - MFD reproduces machine-precision exact solutions on star-shaped meshes
// - MFD reproduces exactly TPFA results when innerProductType="TPFA"


//# Condition Number Estimates for Methods and Meshes
//
//| Method / Mesh | polyhedral_voronoi_complex.vtu | polyhedral_voronoi_lattice.vtu | polyhedral_voronoi_regular.vtu |
//|---|---:|---:|---:|
//| MFD(π) — TPFA | 3.0077238.e5 | 9.6678605.e2 | 2.5071172.e2 |
//| MFD(π) — QuasiTPFA | 6.6168582.e5 | 6.966971.e2 | 2.5071172.e2 |
//| MFD(π) — QuasiRT | 7.2227334.e5 | 1.2007323.e3 | 3.18867.e2 |
//| MFD(π) — Simple | 3.1486776.e6 | 8.7419486.e2 | 3.18867.e2 |
//| MFD(π) — BdVLM | 8.7411124.e6 | 6.1777283.e2 | 2.5071172.e2 |
//| MFD(p) — TPFA | 2.5071422.e3 | 1.6605592.e2 | 3.7053199.e1 |
//| MFD(p) — QuasiTPFA | 3.2336004.e3 | 7.2493851.e1 | 3.7053199.e1 |
//| MFD(p) — QuasiRT | 7.1975418.e3 | 1.7886776.e2 | 9.6124441.e1 |
//| MFD(p) — Simple | 6.3551279.e3 | 1.5516863.e2 | 9.6124441.e1 |
//| MFD(p) — BdVLM | 2.3821719.e3 | 6.7185947.e1 | 3.7053199.e1 |
//| FV — TPFA | 2.5071422.e3 | 1.6605592.e2 | 3.7053199.e1 |
//
//
//Note:
//
//MFD(π) Method
//For the MFD(π) method, the condition number estimate is computed on a system expressed only in terms of Lagrange multipliers. In the
// current construction of the method without applying upwinding in the elliptic components, the cell pressure block is diagonal. This
// structure allows the block to be eliminated exactly the cell pressure.
//MFD(p) Method
//In contrast, for the MFD(p) method the approach is reversed: the Lagrange multiplier is eliminated, and the cell-center pressure is kept.
// This elimination is performed directly in the current setting due to the smaller system size, this is not intended as a general
// recommendation to invert sparse blocks exactly. Instead, from a preconditioning perspective, one could approximate the inverse
// efficiently.
//Regarding mesh effects. all material parameters are set to unit values, resulting in each method constructing a discretization of the
// Laplacian. In this manner,
// the observed variations in the condition number estimates are attributable to mesh distortion effects only.

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

// Pressure L2 error tolerance
static constexpr real64 PRESSURE_L2_TOLERANCE = 1.0e-10;

// Maximum time step for events / solver steps (in seconds)
static constexpr real64 TIME_STEP = 1.0; // 1 day

static constexpr auto TPFA      = "TPFA";
static constexpr auto QuasiTPFA = "quasiTPFA";
static constexpr auto QuasiRT   = "quasiRT";
static constexpr auto Simple    = "simple";
static constexpr auto BdVLM     = "beiraoDaVeigaLipnikovManzini";


std::string generateXmlInputTPFA( std::string const & meshFile )
{
  std::ostringstream oss;
  oss << R"xml(
  <Problem>
    <Mesh>
      <VTKMesh
        name="mesh"  
        partitionRefinement="0"
        useGlobalIds="0"
        file=")xml" << meshFile <<
    R"xml("/>
    </Mesh>

    <Geometry>
      <Box name="westBC" xMin="{ -0.001, 0.0, 0.0}" xMax="{ +0.001, 1.0, 1.0}"/>
      <Box name="eastBC" xMin="{ +0.999, 0.0, 0.0}" xMax="{ +1.001, 1.0, 1.0}"/>
    </Geometry>

    <ElementRegions>
      <CellElementRegion name="Domain" cellBlocks="{ * }" materialList="{rock, fluid }"/>
    </ElementRegions>

    <Solvers gravityVector="{ 0.0, 0.0, 0.0}"> </Solvers>

    <Constitutive>
      <CompressibleSinglePhaseFluid name="fluid"
        defaultDensity="1.0" defaultViscosity="1.0"
        referenceDensity="1.0" referenceViscosity="1.0"
        referencePressure="0.0" compressibility="0.0" viscosibility="0.0"/>
      <CompressibleSolidConstantPermeability name="rock"
        solidModelName="nullSolid" porosityModelName="rockPorosity" permeabilityModelName="rockPerm"/>
      <NullModel name="nullSolid"/>
      <PressurePorosity name="rockPorosity"
        defaultReferencePorosity="0.1" referencePressure="0.0" compressibility="0.0"/>
      <ConstantPermeability name="rockPerm"
        permeabilityComponents="{ 1.0, 1.0, 1.0 }"/>
    </Constitutive>

    <FieldSpecifications>
      <FieldSpecification name="initialPressure" initialCondition="1"
        setNames="{ all }" objectPath="ElementRegions/Domain" fieldName="pressure" scale="1.0"/>  
      <FieldSpecification name="west_pressure"
        setNames="{ westBC }" objectPath="faceManager" fieldName="pressure" scale="2.0"/>
      <FieldSpecification name="east_pressure"
        setNames="{ eastBC }" objectPath="faceManager" fieldName="pressure" scale="1.0"/>      
    </FieldSpecifications>

    <NumericalMethods>
      <FiniteVolume>
        <TwoPointFluxApproximation name="singlePhaseTPFA"/>
      </FiniteVolume>
    </NumericalMethods>

    <Solvers>
      <SinglePhaseFVM name="SinglePhaseFlow" logLevel="1"
        discretization="singlePhaseTPFA" targetRegions="{ Domain }">
        <NonlinearSolverParameters newtonTol="1.0e-8" newtonMaxIter="2"/>
        <LinearSolverParameters
          directParallel="0"/>
      </SinglePhaseFVM>
    </Solvers>

  </Problem>
  )xml";
  return oss.str();
}

// Verifies that the standard TPFA solver produces consistent pressure fields
// on k-orthogonal meshes. L2 error is checked against the analytical linear pressure field.
class TPFAIntegrationTest : public ::testing::TestWithParam< const char * >
{
public:
  TPFAIntegrationTest()
    : state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) ) {}

protected:
  void SetUp() override
  {
    // Use the CMAKE-defined TEST_BINARY_DIR variable
    testBinaryDir = TEST_BINARY_DIR;

    std::string meshFile = testBinaryDir + "/" + GetParam();
    std::string xmlInput = generateXmlInputTPFA( meshFile );
    setupProblemFromXML( state.getProblemManager(), xmlInput.c_str());
  }

  GeosxState state;
  std::string testBinaryDir;
};

INSTANTIATE_TEST_SUITE_P(
  MeshFiles,
  TPFAIntegrationTest,
  ::testing::Values(
    "polyhedral_voronoi_complex.vtu",
    "polyhedral_voronoi_lattice.vtu",
    "polyhedral_voronoi_regular.vtu"
    )
  );

TEST_P( TPFAIntegrationTest, PressureFieldL2Error )
{
  ProblemManager & problemManager = state.getProblemManager();
  DomainPartition & domain = problemManager.getDomainPartition();

  // Retrieve the solver using the PhysicsSolverManager
  SinglePhaseFVM< SinglePhaseBase > & solver =
    dynamic_cast< SinglePhaseFVM< SinglePhaseBase > & >( problemManager.getPhysicsSolverManager().getGroup< SinglePhaseFVM< SinglePhaseBase > >( "SinglePhaseFlow" ) );

  // Run the simulation to compute the numerical pressure
  solver.execute( 0.0, TIME_STEP, 0, 0, 0, domain );

  // Access the mesh and subregion
  MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();
  CellElementSubRegion & subRegion = mesh.getElemManager().getRegion( 0 ).getSubRegion< CellElementSubRegion >( 0 );

  // Retrieve pressure field and cell centers
  arrayView2d< real64 const > centers = subRegion.getElementCenter();
  arrayView1d< real64 const > volumes = subRegion.getElementVolume();
  arrayView1d< real64 const > const p_h = subRegion.getField< fields::flow::pressure >();

  // Compute exact pressure and L2 error
  RAJA::ReduceSum< parallelDeviceReduce, real64 > l2Error_ReduceSum( 0.0 );
  RAJA::ReduceSum< parallelDeviceReduce, real64 > totalVolume_ReduceSum( 0.0 );
  localIndex const n_cells  = subRegion.size();
  forAll< geos::parallelDevicePolicy<> >( n_cells, [=] GEOS_HOST_DEVICE ( localIndex const i )
  {
    real64 x = centers[i][0];
    real64 volume = volumes[i];
    real64 pNumeric = p_h[i];
    real64 pExact = 2.0 * (1.0 - x) + 1.0 * x;
    l2Error_ReduceSum += (pNumeric - pExact) * (pNumeric - pExact) * volume;
    totalVolume_ReduceSum += volume;
  } );

  real64 const data[2] = { l2Error_ReduceSum.get(), totalVolume_ReduceSum.get() };
  real64 l2Error = std::sqrt( data[0] ) / data[1];

  std::string meshFile = GetParam();
  if( meshFile == "polyhedral_voronoi_regular.vtu" )
  {
    // Assert that the L2 error is within machine precision
    EXPECT_NEAR( l2Error, 0.0, PRESSURE_L2_TOLERANCE );
  }
  else
  {
    // Assert that the L2 error is not exact
    EXPECT_GT( l2Error, PRESSURE_L2_TOLERANCE );
  }

}

std::string generateXmlInputMFD( std::string const & innerProductType,
                                 std::string const & meshFile )
{
  std::ostringstream oss;
  oss << R"xml(
  <Problem>

  <Mesh>
    <VTKMesh
      name="mesh"
      logLevel="5"  
      partitionRefinement="0"
      useGlobalIds="0"
      file=")xml" << meshFile <<
    R"xml("/>
  </Mesh>

  <Geometry>
    <Box name="westBC" xMin="{ -0.001, 0.0, 0.0}" xMax="{ +0.001, 1.0, 1.0}"/>
    <Box name="eastBC" xMin="{ +0.999, 0.0, 0.0}" xMax="{ +1.001, 1.0, 1.0}"/>
  </Geometry>

  <ElementRegions>
    <CellElementRegion name="Domain" cellBlocks="{ * }" materialList="{rock, fluid }"/>
  </ElementRegions>

  <Solvers gravityVector="{ 0.0, 0.0, 0.0}"> </Solvers>

  <Constitutive>
    <CompressibleSinglePhaseFluid name="fluid"
      defaultDensity="1.0" defaultViscosity="1.0"
      referenceDensity="1.0" referenceViscosity="1.0"
      referencePressure="0.0" compressibility="0.0" viscosibility="0.0"/>

    <CompressibleSolidConstantPermeability name="rock"
      solidModelName="nullSolid" porosityModelName="rockPorosity" permeabilityModelName="rockPerm"/>

    <NullModel name="nullSolid"/>

    <PressurePorosity name="rockPorosity"
      defaultReferencePorosity="0.1" referencePressure="0.0" compressibility="0.0"/>

    <ConstantPermeability name="rockPerm"
      permeabilityComponents="{ 1.0, 1.0, 1.0 }"/>
  </Constitutive>

  <FieldSpecifications>
    <FieldSpecification name="initialPressure" initialCondition="1"
      setNames="{ all }" objectPath="ElementRegions/Domain" fieldName="pressure" scale="1.0"/>
    <FieldSpecification name="initialFacePressure" initialCondition="1"
      setNames="{ all }" objectPath="faceManager" fieldName="facePressure" scale="1.0"/>          
    <FieldSpecification name="west_pressure"
      setNames="{ westBC }" objectPath="faceManager" fieldName="bcPressure" scale="2.0"/>
    <FieldSpecification name="east_pressure"
      setNames="{ eastBC }" objectPath="faceManager" fieldName="bcPressure" scale="1.0"/>      
  </FieldSpecifications>

  <NumericalMethods>
    <FiniteVolume>
      <HybridMimeticDiscretization
        name="singlePhaseMFD"
        innerProductType=")xml"
      << innerProductType <<
    R"xml("/>
    </FiniteVolume>
  </NumericalMethods>

  <Solvers>
    <SinglePhaseHybridFVM name="SinglePhaseFlow" logLevel="1"
      discretization="singlePhaseMFD" targetRegions="{ Domain }">
      <NonlinearSolverParameters newtonTol="1.0e-8" newtonMaxIter="2"/>
      <LinearSolverParameters
        directParallel="0"/>
    </SinglePhaseHybridFVM>
  </Solvers>

  </Problem>
  )xml";

  return oss.str();
}

// Verifies MFD solver for various inner product types produces exact
// pressure fields for all test meshes. L2 error is compared with exact solution.
using MFDParams = std::tuple< const char *, const char * >;
class MFDIntegrationTest : public ::testing::TestWithParam< MFDParams >
{
public:
  MFDIntegrationTest()
    : state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) ) {}

protected:
  void SetUp() override
  {
    // Use the CMAKE-defined TEST_BINARY_DIR variable
    testBinaryDir = TEST_BINARY_DIR;

    auto [innerProduct, meshFile] = GetParam();
    std::string xmlInput = generateXmlInputMFD( innerProduct, testBinaryDir + "/" + meshFile );
    setupProblemFromXML( state.getProblemManager(), xmlInput.c_str() );
  }

  GeosxState state;
  std::string testBinaryDir;
};


INSTANTIATE_TEST_SUITE_P(
  InnerProductAndMeshes,
  MFDIntegrationTest,
  ::testing::Combine(
    ::testing::Values( TPFA, QuasiTPFA, QuasiRT, Simple, BdVLM ),
    ::testing::Values(
      "polyhedral_voronoi_complex.vtu",
      "polyhedral_voronoi_lattice.vtu",
      "polyhedral_voronoi_regular.vtu"
      )
    )
  );


TEST_P( MFDIntegrationTest, PressureFieldL2Error )
{
  ProblemManager & problemManager = state.getProblemManager();
  DomainPartition & domain = problemManager.getDomainPartition();

  // Retrieve the solver using the PhysicsSolverManager
  SinglePhaseHybridFVM & solver = dynamic_cast< SinglePhaseHybridFVM & >( problemManager.getPhysicsSolverManager().getGroup< SinglePhaseHybridFVM >( "SinglePhaseFlow" ) );

  // Run the simulation to compute the numerical pressure
  solver.execute( 0.0, TIME_STEP, 0, 0, 0, domain );

  // Access the mesh and subregion
  MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();
  CellElementSubRegion & subRegion = mesh.getElemManager().getRegion( 0 ).getSubRegion< CellElementSubRegion >( 0 );

  // Retrieve pressure field and cell centers
  arrayView2d< real64 const > centers = subRegion.getElementCenter();
  arrayView1d< real64 const > volumes = subRegion.getElementVolume();
  arrayView1d< real64 const > const p_h = subRegion.getField< fields::flow::pressure >();

  // Compute exact pressure and L2 error
  RAJA::ReduceSum< parallelDeviceReduce, real64 > l2Error_ReduceSum( 0.0 );
  RAJA::ReduceSum< parallelDeviceReduce, real64 > totalVolume_ReduceSum( 0.0 );
  localIndex const n_cells  = subRegion.size();
  forAll< geos::parallelDevicePolicy<> >( n_cells, [=] GEOS_HOST_DEVICE ( localIndex const i )
  {
    real64 x = centers[i][0];
    real64 volume = volumes[i];
    real64 pNumeric = p_h[i];
    real64 pExact = 2.0 * (1.0 - x) + 1.0 * x;
    l2Error_ReduceSum += (pNumeric - pExact) * (pNumeric - pExact) * volume;
    totalVolume_ReduceSum += volume;
  } );

  real64 const data[2] = { l2Error_ReduceSum.get(), totalVolume_ReduceSum.get() };
  real64 l2Error = std::sqrt( data[0] ) / data[1];

  auto [innerProduct, meshFile] = GetParam();
  if( innerProduct == TPFA and std::string( meshFile ) != "polyhedral_voronoi_regular.vtu" )
  {
    // Assert that the L2 error is not exact
    EXPECT_GT( l2Error, PRESSURE_L2_TOLERANCE );
  }
  else
  {
    // Assert that the L2 error is within machine precision
    EXPECT_NEAR( l2Error, 0.0, PRESSURE_L2_TOLERANCE );
  }
}

// Ensures that MFD with innerProductType="TPFA" reproduces exactly the
// same pressure field as the standard TPFA solver for each mesh.
// This test guarantees solver consistency between TPFA and MFD formulations.
//
class TPFAvsMFDTPFATest : public ::testing::TestWithParam< const char * >
{
protected:
  TPFAvsMFDTPFATest() = default;
};

// Helper function to copy arrayView1d<real64 const> to stdVector<real64>
static inline stdVector< real64 > arrayViewToVector( arrayView1d< real64 const > & arr, localIndex n )
{
  arr.move( hostMemorySpace, false );
  return stdVector< real64 >( arr.data(), arr.data() + n );
}


// Instantiate parameterized test for all mesh files
INSTANTIATE_TEST_SUITE_P(
  MeshFiles,
  TPFAvsMFDTPFATest,
  ::testing::Values(
    "polyhedral_voronoi_complex.vtu",
    "polyhedral_voronoi_lattice.vtu",
    "polyhedral_voronoi_regular.vtu"
    )
  );

TEST_P( TPFAvsMFDTPFATest, PressureFieldComparison )
{
  const char * meshFile = GetParam();
  // Use the CMAKE-defined TEST_BINARY_DIR variable
  std::string testBinaryDir = TEST_BINARY_DIR;

  stdVector< real64 > p_tpfa;
  stdVector< real64 > p_mfd;
  localIndex n_data_tpfa = 0;
  localIndex n_data_mfd = 0;

  // --- Run TPFA solver ---
  {
    GeosxState tpfaState( std::make_unique< CommandLineOptions >( g_commandLineOptions ));

    std::string xmlTPFA = generateXmlInputTPFA( testBinaryDir + "/" + meshFile );
    setupProblemFromXML( tpfaState.getProblemManager(), xmlTPFA.c_str());

    ProblemManager & pmTPFA = tpfaState.getProblemManager();
    DomainPartition & domainTPFA = pmTPFA.getDomainPartition();

    auto & solverTPFA =
      dynamic_cast< SinglePhaseFVM< SinglePhaseBase > & >(
        pmTPFA.getPhysicsSolverManager().getGroup< SinglePhaseFVM< SinglePhaseBase > >( "SinglePhaseFlow" ));

    solverTPFA.setupSystem( domainTPFA, solverTPFA.getDofManager(),
                            solverTPFA.getLocalMatrix(), solverTPFA.getSystemRhs(),
                            solverTPFA.getSystemSolution());
    solverTPFA.implicitStepSetup( 0.0, TIME_STEP, domainTPFA );
    solverTPFA.solverStep( 0.0, TIME_STEP, 0, domainTPFA );
    solverTPFA.implicitStepComplete( 0.0, TIME_STEP, domainTPFA );

    MeshLevel & meshTPFA = domainTPFA.getMeshBody( 0 ).getBaseDiscretization();
    CellElementSubRegion & subRegionTPFA =
      meshTPFA.getElemManager().getRegion( 0 ).getSubRegion< CellElementSubRegion >( 0 );

    arrayView1d< real64 const > p_h = subRegionTPFA.getField< fields::flow::pressure >();
    n_data_tpfa = subRegionTPFA.size();

    p_tpfa = arrayViewToVector( p_h, n_data_tpfa );

    // tpfaState destroyed here — CommunicationTools cleaned up
  }

  // --- Run MFD solver with innerProductType=TPFA ---
  {
    GeosxState mfdState( std::make_unique< CommandLineOptions >( g_commandLineOptions ));

    std::string xmlMFD = generateXmlInputMFD( TPFA, testBinaryDir + "/" + meshFile );
    setupProblemFromXML( mfdState.getProblemManager(), xmlMFD.c_str());

    ProblemManager & pmMFD = mfdState.getProblemManager();
    DomainPartition & domainMFD = pmMFD.getDomainPartition();

    auto & solverMFD =
      dynamic_cast< SinglePhaseHybridFVM & >(
        pmMFD.getPhysicsSolverManager().getGroup< SinglePhaseHybridFVM >( "SinglePhaseFlow" ));

    solverMFD.setupSystem( domainMFD, solverMFD.getDofManager(),
                           solverMFD.getLocalMatrix(), solverMFD.getSystemRhs(),
                           solverMFD.getSystemSolution());
    solverMFD.implicitStepSetup( 0.0, TIME_STEP, domainMFD );
    solverMFD.solverStep( 0.0, TIME_STEP, 0, domainMFD );
    solverMFD.implicitStepComplete( 0.0, TIME_STEP, domainMFD );

    MeshLevel & meshMFD = domainMFD.getMeshBody( 0 ).getBaseDiscretization();
    CellElementSubRegion & subRegionMFD =
      meshMFD.getElemManager().getRegion( 0 ).getSubRegion< CellElementSubRegion >( 0 );

    arrayView1d< real64 const > p_h = subRegionMFD.getField< fields::flow::pressure >();
    n_data_mfd = subRegionMFD.size();

    p_mfd = arrayViewToVector( p_h, n_data_mfd );

    // mfdState destroyed here
  }

  // --- Compare cellwise pressures ---
  ASSERT_EQ( n_data_tpfa, n_data_mfd );
  localIndex const n_cells = n_data_tpfa;
  for( localIndex i = 0; i < n_cells; ++i )
  {
    real64 p_num_tpfa = p_tpfa[i];
    real64 p_num_mfd  = p_mfd[i];
    real64 maxPressure = std::abs( p_num_tpfa - p_num_mfd );
    // Assert that the maxPressure is within machine precision
    EXPECT_NEAR( maxPressure, 0.0, PRESSURE_L2_TOLERANCE );
  }
}


int main( int argc, char * *argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
