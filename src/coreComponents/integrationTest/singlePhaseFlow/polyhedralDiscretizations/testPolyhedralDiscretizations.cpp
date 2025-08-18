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
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/DomainPartition.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVM.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseHybridFVM.hpp"

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

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
        logLevel="5"  
        partitionRefinement="0"
        useGlobalIds="0"
        file=")xml" << meshFile << R"xml("/>
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
        defaultDensity="1000" defaultViscosity="0.001"
        referencePressure="0.0" compressibility="0.0" viscosibility="0.0"/>
      <CompressibleSolidConstantPermeability name="rock"
        solidModelName="nullSolid" porosityModelName="rockPorosity" permeabilityModelName="rockPerm"/>
      <NullModel name="nullSolid"/>
      <PressurePorosity name="rockPorosity"
        defaultReferencePorosity="0.1" referencePressure="0.0" compressibility="0.0"/>
      <ConstantPermeability name="rockPerm"
        permeabilityComponents="{ 1.0e-13, 1.0e-13, 1.0e-13 }"/>
    </Constitutive>

    <FieldSpecifications>
      <FieldSpecification name="initialPressure" initialCondition="1"
        setNames="{ all }" objectPath="ElementRegions/Domain" fieldName="pressure" scale="1.0e7"/>    
      <FieldSpecification name="west_pressure"
        setNames="{ westBC }" objectPath="faceManager" fieldName="pressure" scale="2.0e7"/>
      <FieldSpecification name="east_pressure"
        setNames="{ eastBC }" objectPath="faceManager" fieldName="pressure" scale="1.0e7"/>      
    </FieldSpecifications>

    <NumericalMethods>
      <FiniteVolume>
        <TwoPointFluxApproximation name="singlePhaseTPFA"/>
      </FiniteVolume>
    </NumericalMethods>

    <Solvers>
      <SinglePhaseFVM name="SinglePhaseFlow" logLevel="1"
        discretization="singlePhaseTPFA" targetRegions="{ Domain }">
        <NonlinearSolverParameters newtonTol="1.0e-5" newtonMaxIter="2"/>
        <LinearSolverParameters directParallel="0"/>
      </SinglePhaseFVM>
    </Solvers>

    <Events minTime="0.0" maxTime="86400">
      <PeriodicEvent name="solverApplications"
        endTime="86400" maxEventDt="86400"
        target="/Solvers/SinglePhaseFlow"/>
    </Events>
  </Problem>
  )xml";
  return oss.str();
}

class TPFAIntegrationTest : public ::testing::TestWithParam<const char *>
{
public:
  TPFAIntegrationTest()
    : state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) ) {}

protected:
  void SetUp() override
  {
    std::string xmlInput = generateXmlInputTPFA( GetParam() );
    setupProblemFromXML( state.getProblemManager(), xmlInput.c_str() );
  }

  GeosxState state;
};

INSTANTIATE_TEST_SUITE_P(
  MeshFiles,
  TPFAIntegrationTest,
  ::testing::Values(
    "polyhedral_voronoi_complex.vtk",
    "polyhedral_voronoi_lattice.vtk",
    "polyhedral_voronoi_regular.vtk"
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
  solver.setupSystem( domain, solver.getDofManager(), solver.getLocalMatrix(), solver.getSystemRhs(), solver.getSystemSolution() );
  solver.implicitStepSetup( 0.0, 86400, domain );
  solver.solverStep( 0.0, 86400, 0, domain );
  solver.implicitStepComplete( 0.0, 86400, domain );

  // Access the mesh and subregion
  MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();
  CellElementSubRegion & subRegion = mesh.getElemManager().getRegion( 0 ).getSubRegion< CellElementSubRegion >( 0 );

  // Retrieve pressure field and cell centers
  arrayView2d< real64 const > centers = subRegion.getElementCenter();
  arrayView1d< real64 const > volumes = subRegion.getElementVolume();
  arrayView1d< real64 const > const p_h = subRegion.getField< fields::flow::pressure >();

  // Compute exact pressure and L2 error
  real64 l2Error = 0.0;
  real64 totalVolume = 0.0;
  for( localIndex i = 0; i < subRegion.size(); ++i )
  {
    real64 x = centers[i][0];
    real64 volume = volumes[i];
    real64 pNumeric = p_h[i] * 1.0e-6; // Convert pressure to MPa
    real64 pExact = 20.0 * (1.0 - x) + 10.0 * x;
    l2Error += (pNumeric - pExact) * (pNumeric - pExact) * volume;
    totalVolume += volume;
  }

  l2Error = std::sqrt( l2Error / totalVolume );
  
  std::string meshFile = GetParam();
  if (meshFile.compare("polyhedral_voronoi_regular.vtk") == 0)
  {
    // Assert that the L2 error is within machine precision
    EXPECT_NEAR( l2Error, 0.0, 1.0e-10 );
  }else{
    // Assert that the L2 error is not exact
    EXPECT_GT( l2Error, 1.0e-10 );
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
      file=")xml" << meshFile << R"xml("/>
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
      defaultDensity="1000" defaultViscosity="0.001"
      referencePressure="0.0" compressibility="0.0" viscosibility="0.0"/>

    <CompressibleSolidConstantPermeability name="rock"
      solidModelName="nullSolid" porosityModelName="rockPorosity" permeabilityModelName="rockPerm"/>

    <NullModel name="nullSolid"/>

    <PressurePorosity name="rockPorosity"
      defaultReferencePorosity="0.1" referencePressure="0.0" compressibility="0.0"/>

    <ConstantPermeability name="rockPerm"
      permeabilityComponents="{ 1.0e-13, 1.0e-13, 1.0e-13 }"/>
  </Constitutive>

  <FieldSpecifications>
    <FieldSpecification name="initialPressure" initialCondition="1"
      setNames="{ all }" objectPath="ElementRegions/Domain" fieldName="pressure" scale="1.0e7"/>    
    <FieldSpecification name="west_pressure"
      setNames="{ westBC }" objectPath="faceManager" fieldName="pressure" scale="2.0e7"/>
    <FieldSpecification name="east_pressure"
      setNames="{ eastBC }" objectPath="faceManager" fieldName="pressure" scale="1.0e7"/>      
  </FieldSpecifications>

  <NumericalMethods>
    <FiniteVolume>
      <HybridMimeticDiscretization
        name="singlePhaseMFD"
        innerProductType=")xml"
    << innerProductType << R"xml("/>
    </FiniteVolume>
  </NumericalMethods>

  <Solvers>
    <SinglePhaseHybridFVM name="SinglePhaseFlow" logLevel="1"
      discretization="singlePhaseMFD" targetRegions="{ Domain }">
      <NonlinearSolverParameters newtonTol="1.0e-5" newtonMaxIter="8"/>
      <LinearSolverParameters directParallel="0"/>
    </SinglePhaseHybridFVM>
  </Solvers>

  <Events minTime="0.0" maxTime="86400">
    <PeriodicEvent name="solverApplications"
      endTime="86400" maxEventDt="86400"
      target="/Solvers/SinglePhaseFlow"/>
  </Events>

  </Problem>
  )xml";

  return oss.str();
}

using MFDParams = std::tuple<const char *, const char *>;

class MFDIntegrationTest : public ::testing::TestWithParam<MFDParams>
{
public:
  MFDIntegrationTest()
    : state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) ) {}

protected:
  void SetUp() override
  {
    auto [innerProduct, meshFile] = GetParam();
    std::string xmlInput = generateXmlInputMFD(innerProduct, meshFile);
    setupProblemFromXML( state.getProblemManager(), xmlInput.c_str() );
  }

  GeosxState state;
};


INSTANTIATE_TEST_SUITE_P(
  InnerProductAndMeshes,
  MFDIntegrationTest,
  ::testing::Combine(
    ::testing::Values(TPFA, QuasiTPFA, QuasiRT, Simple, BdVLM),
    ::testing::Values(
      "polyhedral_voronoi_complex.vtk",
      "polyhedral_voronoi_lattice.vtk",
      "polyhedral_voronoi_regular.vtk"
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
  solver.setupSystem( domain, solver.getDofManager(), solver.getLocalMatrix(), solver.getSystemRhs(), solver.getSystemSolution() );
  solver.implicitStepSetup( 0.0, 86400, domain );
  solver.solverStep( 0.0, 86400, 0, domain );
  solver.implicitStepComplete( 0.0, 86400, domain );

  // Access the mesh and subregion
  MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();
  CellElementSubRegion & subRegion = mesh.getElemManager().getRegion( 0 ).getSubRegion< CellElementSubRegion >( 0 );

  // Retrieve pressure field and cell centers
  arrayView2d< real64 const > centers = subRegion.getElementCenter();
  arrayView1d< real64 const > volumes = subRegion.getElementVolume();
  arrayView1d< real64 const > const p_h = subRegion.getField< fields::flow::pressure >();

  // Compute exact pressure and L2 error
  real64 l2Error = 0.0;
  real64 totalVolume = 0.0;
  for( localIndex i = 0; i < subRegion.size(); ++i )
  {
    real64 x = centers[i][0];
    real64 volume = volumes[i];
    real64 pNumeric = p_h[i] * 1.0e-6; // Convert pressure to MPa
    real64 pExact = 20.0 * (1.0 - x) + 10.0 * x;
    l2Error += (pNumeric - pExact) * (pNumeric - pExact) * volume;
    totalVolume += volume;
  }

  l2Error = std::sqrt( l2Error / totalVolume );

  auto [innerProduct, meshFile] = GetParam();
  if (innerProduct == TPFA and std::string(meshFile).compare("polyhedral_voronoi_regular.vtk") != 0)
  {
    // Assert that the L2 error is not exact
    EXPECT_GT(l2Error, 1.0e-10 );
  }else{
    // Assert that the L2 error is within machine precision
    EXPECT_NEAR( l2Error, 0.0, 1.0e-10 );
  }
}

// cross-check test. Ensure that MFD with innerProductType="TPFA" produces exactly the same pressure field as the TPFA solver

class TPFAvsMFDTPFATest : public ::testing::TestWithParam<const char *>
{
public:
  TPFAvsMFDTPFATest()
    : tpfaState( std::make_unique< CommandLineOptions >( g_commandLineOptions ) ),
      mfdState( std::make_unique< CommandLineOptions >( g_commandLineOptions ) ) {}

protected:
  void SetUp() override
  {
    meshFile = GetParam();
  }

  std::string meshFile;
  GeosxState tpfaState;
  GeosxState mfdState;
};

INSTANTIATE_TEST_SUITE_P(
  MeshFiles,
  TPFAvsMFDTPFATest,
  ::testing::Values(
    "polyhedral_voronoi_complex.vtk",
    "polyhedral_voronoi_lattice.vtk",
    "polyhedral_voronoi_regular.vtk"
  )
);


TEST_P( TPFAvsMFDTPFATest, PressureFieldComparison )
{
  arrayView1d< real64> p_tpfa;
  arrayView1d< real64> p_mfd;
  geos::localIndex n_data_tpfa = 0;
  geos::localIndex n_data_mfd = 0;

  const char* meshFile = GetParam();

  // --- Run TPFA solver ---
  {
    GeosxState tpfaState( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
    std::string xmlTPFA = generateXmlInputTPFA( meshFile );
    setupProblemFromXML( tpfaState.getProblemManager(), xmlTPFA.c_str() );

    ProblemManager & pmTPFA = tpfaState.getProblemManager();
    DomainPartition & domainTPFA = pmTPFA.getDomainPartition();

    auto & solverTPFA =
      dynamic_cast< SinglePhaseFVM< SinglePhaseBase > & >(
        pmTPFA.getPhysicsSolverManager().getGroup< SinglePhaseFVM< SinglePhaseBase > >( "SinglePhaseFlow" ) );

    solverTPFA.setupSystem( domainTPFA, solverTPFA.getDofManager(),
                            solverTPFA.getLocalMatrix(), solverTPFA.getSystemRhs(), solverTPFA.getSystemSolution() );
    solverTPFA.implicitStepSetup( 0.0, 86400, domainTPFA );
    solverTPFA.solverStep( 0.0, 86400, 0, domainTPFA );
    solverTPFA.implicitStepComplete( 0.0, 86400, domainTPFA );

    MeshLevel & meshTPFA = domainTPFA.getMeshBody( 0 ).getBaseDiscretization();
    CellElementSubRegion & subRegionTPFA =
      meshTPFA.getElemManager().getRegion( 0 ).getSubRegion< CellElementSubRegion >( 0 );

    p_tpfa = std::move(subRegionTPFA.getField< fields::flow::pressure >());
    n_data_tpfa = subRegionTPFA.size();
  } // <--- tpfaState destroyed here, CommunicationTools cleaned up

  // --- Run MFD solver with innerProductType=TPFA ---
  {
    GeosxState mfdState( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
    std::string xmlMFD = generateXmlInputMFD( TPFA, meshFile );
    setupProblemFromXML( mfdState.getProblemManager(), xmlMFD.c_str() );

    ProblemManager & pmMFD = mfdState.getProblemManager();
    DomainPartition & domainMFD = pmMFD.getDomainPartition();

    auto & solverMFD =
      dynamic_cast< SinglePhaseHybridFVM & >(
        pmMFD.getPhysicsSolverManager().getGroup< SinglePhaseHybridFVM >( "SinglePhaseFlow" ) );

    solverMFD.setupSystem( domainMFD, solverMFD.getDofManager(),
                           solverMFD.getLocalMatrix(), solverMFD.getSystemRhs(), solverMFD.getSystemSolution() );
    solverMFD.implicitStepSetup( 0.0, 86400, domainMFD );
    solverMFD.solverStep( 0.0, 86400, 0, domainMFD );
    solverMFD.implicitStepComplete( 0.0, 86400, domainMFD );

    MeshLevel & meshMFD = domainMFD.getMeshBody( 0 ).getBaseDiscretization();
    CellElementSubRegion & subRegionMFD =
      meshMFD.getElemManager().getRegion( 0 ).getSubRegion< CellElementSubRegion >( 0 );

    p_mfd = std::move(subRegionMFD.getField< fields::flow::pressure >());
    n_data_mfd = subRegionMFD.size();
  }

  // --- Compare cellwise pressures ---
  ASSERT_EQ( n_data_tpfa, n_data_mfd );
  for( localIndex i = 0; i < n_data_tpfa; ++i )
  {
    real64 p_num_tpfa = p_tpfa[i];
    real64 p_num_mfd  = p_mfd[i];
    real64 p_diff     = (p_num_tpfa - p_num_mfd) * 1.0e-6; // Convert pressure to MPa
    EXPECT_NEAR( p_diff, 0.0, 1.0e-10 ) << "Mismatch at cell " << i;
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
