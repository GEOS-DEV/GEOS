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


#include "integrationTests/fluidFlowTests/testCompFlowUtils.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "mainInterface/initialization.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include <gtest/gtest.h>

using namespace geos;
CommandLineOptions g_commandLineOptions;

char const * xmlInput =
  R"xml(
<Problem>
  <Solvers gravityVector="{ 0, 0, 0 }">
    <SolidMechanicsLagrangianFEM name="mechanics" discretization="FE1"
      targetRegions="{ Region }" timeIntegrationOption="QuasiStatic"/>
  </Solvers>
  <Mesh>
    <InternalMesh name="mesh" elementTypes="{ C3D8 }"
      xCoords="{ 0, 2 }" yCoords="{ 0, 1 }" zCoords="{ 0, 1 }"
      nx="{ 2 }" ny="{ 1 }" nz="{ 1 }" cellBlockNames="{ cells }"/>
  </Mesh>
  <NumericalMethods>
    <FiniteElements><FiniteElementSpace name="FE1" order="1"/></FiniteElements>
  </NumericalMethods>
  <ElementRegions>
    <CellElementRegion name="Region" cellBlocks="{ cells }" materialList="{ rock }"/>
  </ElementRegions>
  <Constitutive>
    <ElasticIsotropic name="rock" defaultDensity="2500"
      defaultBulkModulus="5e9" defaultShearModulus="3e9"/>
  </Constitutive>
  <FieldSpecifications>
    <FieldSpecification name="sxx" initialCondition="1" setNames="{ all }"
      objectPath="ElementRegions/Region" fieldName="rock_stress" component="0" scale="1e6"/>
    <FieldSpecification name="syy" initialCondition="1" setNames="{ all }"
      objectPath="ElementRegions/Region" fieldName="rock_stress" component="1" scale="2e6"/>
    <FieldSpecification name="szz" initialCondition="1" setNames="{ all }"
      objectPath="ElementRegions/Region" fieldName="rock_stress" component="2" scale="3e6"/>
    <FieldSpecification name="syz" initialCondition="1" setNames="{ all }"
      objectPath="ElementRegions/Region" fieldName="rock_stress" component="3" scale="4e6"/>
    <FieldSpecification name="sxz" initialCondition="1" setNames="{ all }"
      objectPath="ElementRegions/Region" fieldName="rock_stress" component="4" scale="5e6"/>
    <FieldSpecification name="sxy" initialCondition="1" setNames="{ all }"
      objectPath="ElementRegions/Region" fieldName="rock_stress" component="5" scale="6e6"/>
  </FieldSpecifications>
  <Events maxTime="1"/>
</Problem>
)xml";

TEST( InitialMechanicalOutput, PrescribedStressBeforeFirstStep )
{
  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  geos::testing::setupProblemFromXML( state.getProblemManager(), xmlInput );
  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();
  auto verify = [&]( real64 const expectedStrain[6] )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
    {
      auto const & stress = subRegion.getField< fields::solidMechanics::averageStress >();
      auto const & strain = subRegion.getField< fields::solidMechanics::averageStrain >();
      auto const & plastic = subRegion.getField< fields::solidMechanics::averagePlasticStrain >();
      stress.move( hostMemorySpace, false );
      strain.move( hostMemorySpace, false );
      plastic.move( hostMemorySpace, false );
      for( localIndex k = 0; k < subRegion.size(); ++k )
      {
        for( integer c = 0; c < 6; ++c )
        {
          EXPECT_NEAR( stress[k][c], (c + 1) * 1e6, 1e-7 );
          EXPECT_NEAR( strain[k][c], expectedStrain[c], 1e-14 );
          EXPECT_EQ( plastic[k][c], 0.0 );
        }
      }
    } );
  };
  real64 const zero[6] = {};
  verify( zero );

  // Repeat finalization with an affine initial displacement. It must update
  // derived strain without adding plastic strain or changing the prestress.
  auto & displacement = mesh.getNodeManager().getField< fields::solidMechanics::totalDisplacement >();
  auto & increment = mesh.getNodeManager().getField< fields::solidMechanics::incrementalDisplacement >();
  auto const & position = mesh.getNodeManager().referencePosition();
  displacement.move( hostMemorySpace, true );
  increment.move( hostMemorySpace, true );
  position.move( hostMemorySpace, false );
  for( localIndex k = 0; k < mesh.getNodeManager().size(); ++k )
  {
    displacement[k][0] = 0.001 * position[k][0] + 0.0004 * position[k][1];
    displacement[k][1] = 0.002 * position[k][1];
    displacement[k][2] = 0.003 * position[k][2];
    for( integer c = 0; c < 3; ++c )
    {
      increment[k][c] = displacement[k][c];
    }
  }
  state.getProblemManager().getPhysicsSolverManager().getGroup< SolidMechanicsLagrangianFEM >( "mechanics" ).finalizeInitialState( domain );
  real64 const affine[6] = { 0.001, 0.002, 0.003, 0, 0, 0.0002 };
  verify( affine );
}

int main( int argc, char * argv[] )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *basicSetup( argc, argv, false );
  int const result = RUN_ALL_TESTS();
  basicCleanup();
  return result;
}
