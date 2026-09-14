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

/**
 * @file testMixedVEMOutputConvention.cpp
 *
 * The mixed VEM cell stress must reach the output files under the same component names and in
 * the same order as the finite element solver, otherwise a viewer shows the shear components
 * under the wrong names.
 */

#include "constitutive/ConstitutiveManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshManager.hpp"

#include <gtest/gtest.h>

#include <vector>

using namespace geos;
using namespace geos::dataRepository;

namespace
{

CommandLineOptions g_commandLineOptions;

constexpr char const * const meshAndMaterial = R"xml(
  <Mesh>
    <InternalMesh
      name="mesh"
      elementTypes="{ C3D8 }"
      xCoords="{ 0, 1 }"
      yCoords="{ 0, 1 }"
      zCoords="{ 0, 1 }"
      nx="{ 2 }"
      ny="{ 2 }"
      nz="{ 2 }"
      cellBlockNames="{ cb1 }"/>
  </Mesh>
  <ElementRegions>
    <CellElementRegion name="Domain" cellBlocks="{ cb1 }" materialList="{ rock }"/>
  </ElementRegions>
  <Constitutive>
    <ElasticIsotropic name="rock" defaultDensity="2700" defaultBulkModulus="5.0e8" defaultShearModulus="3.0e8"/>
  </Constitutive>
  <Events maxTime="1.0">
    <PeriodicEvent name="solve" forceDt="1.0" target="/Solvers/solver"/>
  </Events>
)xml";

constexpr char const * const finiteElementSolver = R"xml(
  <Solvers gravityVector="{ 0.0, 0.0, 0.0 }">
    <SolidMechanicsLagrangianFEM name="solver" discretization="FE1" targetRegions="{ Domain }"/>
  </Solvers>
  <NumericalMethods>
    <FiniteElements>
      <FiniteElementSpace name="FE1" order="1"/>
    </FiniteElements>
  </NumericalMethods>
)xml";

constexpr char const * const mixedVEMSolver = R"xml(
  <Solvers gravityVector="{ 0.0, 0.0, 0.0 }">
    <SolidMechanicsMixedVEM name="solver" discretization="mixedVEM1" targetRegions="{ Domain }"/>
  </Solvers>
  <NumericalMethods>
    <MixedVEM>
      <MixedVEMDiscretization name="mixedVEM1" hybridization="0"/>
    </MixedVEM>
  </NumericalMethods>
)xml";

void setupProblemFromXML( ProblemManager & problemManager, string const & xmlInput )
{
  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult const xmlResult = xmlDocument.loadString( xmlInput.c_str() );
  ASSERT_TRUE( xmlResult ) << xmlResult.description();

  Group & commandLine = problemManager.getGroup< Group >( problemManager.groupKeys.commandLine );
  commandLine.registerWrapper< integer >( problemManager.viewKeys.xPartitionsOverride.key() ).
    setApplyDefaultValue( MpiWrapper::commSize( MPI_COMM_GEOS ) );

  xmlWrapper::xmlNode xmlProblemNode = xmlDocument.getChild( keys::ProblemManager );
  problemManager.processInputFileRecursive( xmlDocument, xmlProblemNode );

  DomainPartition & domain = problemManager.getDomainPartition();

  constitutive::ConstitutiveManager & constitutiveManager = domain.getConstitutiveManager();
  xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child( constitutiveManager.getName().c_str() );
  constitutiveManager.processInputFileRecursive( xmlDocument, topLevelNode );

  MeshManager & meshManager = problemManager.getGroup< MeshManager >( problemManager.groupKeys.meshManager );
  meshManager.generateMeshLevels( domain );

  ElementRegionManager & elementManager = domain.getMeshBody( 0 ).getBaseDiscretization().getElemManager();
  topLevelNode = xmlProblemNode.child( elementManager.getName().c_str() );
  elementManager.processInputFileRecursive( xmlDocument, topLevelNode );

  problemManager.problemSetup();
}

/// Component names of a cell field, or of the per quadrature point stress of the material.
struct StressLabels
{
  std::vector< string > cell;
  std::vector< string > material;
};

StressLabels readLabels( char const * const solverBlock, string const & cellField )
{
  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  ProblemManager & problemManager = state.getProblemManager();
  setupProblemFromXML( problemManager, string( "<Problem>" ) + solverBlock + meshAndMaterial + "</Problem>" );

  StressLabels labels;

  ElementRegionManager & elemManager =
    problemManager.getDomainPartition().getMeshBody( 0 ).getBaseDiscretization().getElemManager();

  elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
  {
    if( !labels.cell.empty() )
    {
      return;
    }

    Span< string const > const cell = subRegion.getWrapperBase( cellField ).getDimLabels( 1 );
    labels.cell.assign( cell.begin(), cell.end() );

    // the material stress is laid out as (cell, quadrature point, component)
    Span< string const > const material =
      subRegion.getConstitutiveModels().getGroup( "rock" ).getWrapperBase( "stress" ).getDimLabels( 2 );
    labels.material.assign( material.begin(), material.end() );
  } );

  return labels;
}

} // namespace

TEST( MixedVEMOutputConvention, cellStressLabelsMatchFiniteElement )
{
  StressLabels const fem = readLabels( finiteElementSolver, "averageStress" );
  StressLabels const vem = readLabels( mixedVEMSolver, "stress" );

  std::vector< string > const voigt = { "XX", "YY", "ZZ", "YZ", "XZ", "XY" };

  // the reference itself, so a change on the finite element side is caught here as well
  EXPECT_EQ( fem.cell, voigt );
  EXPECT_EQ( fem.material, voigt );

  EXPECT_EQ( vem.cell, fem.cell );
  EXPECT_EQ( vem.cell, vem.material );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
