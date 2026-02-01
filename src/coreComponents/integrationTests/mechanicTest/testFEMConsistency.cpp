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
#include <fstream>
#include <tuple>
#include <limits>
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/initialization.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "mesh/SurfaceElementSubRegion.hpp"
#include "mesh/FaceElementSubRegion.hpp"
#include "physicsSolvers/solidMechanics/contact/ContactFields.hpp"
#include "common/format/Format.hpp"
#include "codingUtilities/Utilities.hpp"
#include "mesh/FaceManager.hpp"
#include "constitutive/solid/SolidFields.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"

using namespace geos;

constexpr real64 relative_tolerance = 1.0e-7;

CommandLineOptions g_commandLineOptions;

class ConsistencyTest : public ::testing::TestWithParam< std::tuple< std::string, real64, real64, real64 > >
{
protected:
  void SetUp() override
  {
    testBinaryDir = TEST_BINARY_DIR;
  }

  std::string generateXmlInput( std::string const & meshFile, std::string const & nodeSetNames, real64 s_xx, real64 s_yy, real64 s_zz )
  {
    std::ostringstream oss;
    oss << R"xml(<?xml version="1.0" ?>
<Problem>
  <Mesh>
    <VTKMesh name="mesh1" file=")xml" << meshFile << R"xml(" nodesetNames=")xml" << nodeSetNames <<
    R"xml("/>
  </Mesh>
  <Geometry>
    <Box name="xnegFace" xMin="{ -0.01, -0.01, -0.01 }" xMax="{  0.01,  1.01,  1.01 }"/>
    <Box name="xposFace" xMin="{  0.99, -0.01, -0.01 }" xMax="{  1.01,  1.01,  1.01 }"/>
    <Box name="ynegFace" xMin="{ -0.01, -0.01, -0.01 }" xMax="{  1.01,  0.01,  1.01 }"/>
    <Box name="yposFace" xMin="{ -0.01,  0.99, -0.01 }" xMax="{  1.01,  1.01,  1.01 }"/>
    <Box name="znegFace" xMin="{ -0.01, -0.01, -0.01 }" xMax="{  1.01,  1.01,  0.01 }"/>
    <Box name="zposFace" xMin="{ -0.01, -0.01,  0.99 }" xMax="{  1.01,  1.01,  1.01 }"/>
    <Box name="fracture" xMin="{ -0.01, -0.01, -0.01 }" xMax="{  1.01,  1.01,  1.01 }"/>
  </Geometry>
  <Solvers gravityVector="{0.0, 0.0, 0.0}">
    <SolidMechanicsAugmentedLagrangianContact 
      name="mechSolver"
      simultaneous="1"
      symmetric="1"
      iterPenaltyN="1.0e1"
      iterPenaltyT="1.0e-1"
      tolNormalTrac="1.e-08"
      tolTauLimit="1.e-08"
      tolJumpN="1.e-8"
      tolJumpT="1.e-8" 
      discretization="FE1" 
      targetRegions="{ Region, Fracture }" 
      logLevel="1">
      <NonlinearSolverParameters newtonTol="1.0e-5" newtonMaxIter="20" logLevel="1"/>
      <LinearSolverParameters solverType="gmres" preconditionerType="amg" krylovTol="1.0e-10" logLevel="1"/>
    </SolidMechanicsAugmentedLagrangianContact>
    <SurfaceGenerator name="SurfaceGen" targetRegions="{ Region, Fracture }" fractureRegion="Fracture" initialRockToughness="10.0e9"/>
  </Solvers>
  <NumericalMethods>
    <FiniteElements>
      <FiniteElementSpace name="FE1" order="1" />
    </FiniteElements>
  </NumericalMethods>
  <ElementRegions>
    <CellElementRegion name="Region" cellBlocks="{ * }" materialList="{ rock }"/>
    <SurfaceElementRegion name="Fracture" defaultAperture="1.0e-4" faceBlock="faceElementSubRegion" materialList="{ fractureContact }"/>
  </ElementRegions>
  <Constitutive>
    <ElasticIsotropic name="rock" defaultDensity="2500" defaultYoungModulus="10e9" defaultPoissonRatio="0.2"/>
    <Coulomb name="fractureContact" cohesion="1.0e20" frictionCoefficient="0.0"/>
  </Constitutive>
  <FieldSpecifications>
    <FieldSpecification name="separableFace" fieldName="isFaceSeparable" initialCondition="1" setNames=")xml" << nodeSetNames <<
    R"xml(" objectPath="faceManager" scale="1" />
    <FieldSpecification name="frac" initialCondition="1" setNames=")xml" << nodeSetNames <<
    R"xml(" objectPath="faceManager" fieldName="ruptureState" scale="1" />
    
    <FieldSpecification name="xneg_disp" component="0" setNames="{ xnegFace }" objectPath="nodeManager" fieldName="totalDisplacement" scale="0.0"/>
    <FieldSpecification name="yneg_disp" component="1" setNames="{ ynegFace }" objectPath="nodeManager" fieldName="totalDisplacement" scale="0.0"/>
    <FieldSpecification name="zneg_disp" component="2" setNames="{ znegFace }" objectPath="nodeManager" fieldName="totalDisplacement" scale="0.0"/>
    
    <Traction name="xpos_traction" setNames="{ xposFace }" objectPath="faceManager" tractionType="normal" scale=")xml" << s_xx <<
    R"xml("/>
    <Traction name="ypos_traction" setNames="{ yposFace }" objectPath="faceManager" tractionType="normal" scale=")xml" << s_yy <<
    R"xml("/>
    <Traction name="zpos_traction" setNames="{ zposFace }" objectPath="faceManager" tractionType="normal" scale=")xml" << s_zz <<
    R"xml("/>
  </FieldSpecifications>
  <Tasks>
    <SolidMechanicsAugmentedLagrangianContactInitialization name="ELASTICITY.PRE.INIT.STEP" solidSolverName="mechSolver" logLevel="1"/>
  </Tasks>  
  <Outputs>
  </Outputs>
  <Events minTime="-1.0e11" maxTime="1.0">
    <PeriodicEvent name="solverApplications" target="/Solvers/mechSolver" forceDt="1.0"/>
  </Events>
</Problem>
)xml";
    return oss.str();
  }

  std::string testBinaryDir;
};

TEST_P( ConsistencyTest, Run )
{
  auto const & params = GetParam();
  auto const & meshFileName = std::get< 0 >( params );
  real64 const s_xx = std::get< 1 >( params );
  real64 const s_yy = std::get< 2 >( params );
  real64 const s_zz = std::get< 3 >( params );
  std::string xmlPath = testBinaryDir + "/test_fem_consistency.xml";

  std::string nodeSetNames = "{ f1_node_set }";
  if( meshFileName.find( "_DFN_2.vtu" ) != std::string::npos )
  {
    nodeSetNames = "{ f2_node_set }";
  }
  else if( meshFileName.find( "_DFN_3.vtu" ) != std::string::npos )
  {
    nodeSetNames = "{ f3_node_set }";
  }
  else if( meshFileName.find( "_DFN_123.vtu" ) != std::string::npos )
  {
    nodeSetNames = "{ f1_node_set, f2_node_set, f3_node_set }";
  }
  else if( meshFileName.find( "_DFN_12.vtu" ) != std::string::npos )
  {
    nodeSetNames = "{ f1_node_set, f2_node_set }";
  }
  else if( meshFileName.find( "_DFN_13.vtu" ) != std::string::npos )
  {
    nodeSetNames = "{ f1_node_set, f3_node_set }";
  }

  std::string xmlContent = generateXmlInput( meshFileName, nodeSetNames, s_xx, s_yy, s_zz );

  {
    std::ofstream ofs( xmlPath );
    ofs << xmlContent;
  }

  auto options = std::make_unique< CommandLineOptions >( g_commandLineOptions );
  options->inputFileNames.push_back( xmlPath );
  options->problemName = "test_fem_consistency";

  // Scoped state to ensure cleanup
  {
    GeosxState state( std::move( options ) );
    ASSERT_TRUE( state.initializeDataRepository() );
    state.applyInitialConditions();
    state.run();

    // Verification using bulk stress and geometric normal
    ProblemManager & pm = state.getProblemManager();
    auto & fractureRegion = pm.getDomainPartition().getMeshBody( "mesh1" ).getBaseDiscretization().getElemManager().template getRegion< SurfaceElementRegion >( "Fracture" );
    auto & elemManager = pm.getDomainPartition().getMeshBody( "mesh1" ).getBaseDiscretization().getElemManager();
    auto & volumeRegion = elemManager.template getRegion< CellElementRegion >( "Region" );

    auto const & faceManager = pm.getDomainPartition().getMeshBody( "mesh1" ).getBaseDiscretization().getFaceManager();
    auto const & nodeManager = pm.getDomainPartition().getMeshBody( "mesh1" ).getBaseDiscretization().getNodeManager();

    auto const & faceToCell = faceManager.elementList();
    auto const & faceToRegion = faceManager.elementRegionList();
    auto const & faceToSubRegion = faceManager.elementSubRegionList();
    auto const & faceNodes = faceManager.nodeList();
    auto const & nodePos = nodeManager.referencePosition();

    fractureRegion.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      auto const & faces = subRegion.faceList();
      localIndex n_facets = faces.size();
      for( localIndex k=0; k < n_facets; ++k )
      {
        // Each fracture element corresponds to a set of faces (usually 2 for split)
        localIndex f = faces[0][k];

        // Compute Geometric Normal
        auto const & fnodes = faceNodes[f];
        if( fnodes.size() < 3 )
          continue;

        real64 nx = 0.0;
        real64 ny = 0.0;
        real64 nz = 0.0;

        if( fnodes.size() == 3 )
        {
          real64 p0[3] = { nodePos[fnodes[0]][0], nodePos[fnodes[0]][1], nodePos[fnodes[0]][2] };
          real64 p1[3] = { nodePos[fnodes[1]][0], nodePos[fnodes[1]][1], nodePos[fnodes[1]][2] };
          real64 p2[3] = { nodePos[fnodes[2]][0], nodePos[fnodes[2]][1], nodePos[fnodes[2]][2] };

          real64 v1[3] = { p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2] };
          real64 v2[3] = { p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2] };

          nx = v1[1]*v2[2] - v1[2]*v2[1];
          ny = v1[2]*v2[0] - v1[0]*v2[2];
          nz = v1[0]*v2[1] - v1[1]*v2[0];
        }
        else
        {
          real64 p0[3] = { nodePos[fnodes[0]][0], nodePos[fnodes[0]][1], nodePos[fnodes[0]][2] };
          real64 p1[3] = { nodePos[fnodes[1]][0], nodePos[fnodes[1]][1], nodePos[fnodes[1]][2] };
          real64 p2[3] = { nodePos[fnodes[2]][0], nodePos[fnodes[2]][1], nodePos[fnodes[2]][2] };
          real64 p3[3] = { nodePos[fnodes[3]][0], nodePos[fnodes[3]][1], nodePos[fnodes[3]][2] };

          real64 v1[3] = { p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2] };
          real64 v2[3] = { p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2] };

          nx = v1[1]*v2[2] - v1[2]*v2[1];
          ny = v1[2]*v2[0] - v1[0]*v2[2];
          nz = v1[0]*v2[1] - v1[1]*v2[0];
        }

        real64 norm = std::sqrt( nx*nx + ny*ny + nz*nz );
        nx /= norm; ny /= norm; nz /= norm;

        // Get neighbor cell
        localIndex neigh = 0;
        localIndex c = faceToCell[f][0];
        if( c == static_cast< localIndex >(-1) )
        {
          c = faceToCell[f][1];
          neigh = 1;
        }

        localIndex er = faceToRegion[f][neigh];
        localIndex esr = faceToSubRegion[f][neigh];

        ElementRegionBase & region = elemManager.getRegion( er );
        ElementSubRegionBase & cellSubRegion = region.getSubRegion( esr );

        // Get average stress from cell subregion
        auto const & avgStress = cellSubRegion.getField< fields::solidMechanics::averageStress >();

        // averageStress is array2d (element, component)
        // Components: XX, YY, ZZ, YZ, XZ, XY
        real64 sig_xx = avgStress( c, 0 );
        real64 sig_yy = avgStress( c, 1 );
        real64 sig_zz = avgStress( c, 2 );
        real64 sig_yz = avgStress( c, 3 );
        real64 sig_xz = avgStress( c, 4 );
        real64 sig_xy = avgStress( c, 5 );       // Warning: check Voigt order in SolidMechanicsLagrangianFEM

        // Compute t_sim = sigma * n
        real64 ts_x = sig_xx * nx + sig_xy * ny + sig_xz * nz;
        real64 ts_y = sig_xy * nx + sig_yy * ny + sig_yz * nz;         // Note symmetry
        real64 ts_z = sig_xz * nx + sig_yz * ny + sig_zz * nz;

        // Compute t_exact = S_input * n
        real64 te_x = s_xx * nx;
        real64 te_y = s_yy * ny;
        real64 te_z = s_zz * nz;

        // Compare (allow slightly larger relative_tolerance due to numerical precision of stress field)
        real64 const err_abs = std::sqrt( std::pow( ts_x - te_x, 2 ) + std::pow( ts_y - te_y, 2 ) + std::pow( ts_z - te_z, 2 ) );
        real64 const norm_te = std::sqrt( std::pow( te_x, 2 ) + std::pow( te_y, 2 ) + std::pow( te_z, 2 ) );
        real64 const err = ( norm_te > 1.0e-16 ) ? err_abs / norm_te : err_abs;
        EXPECT_LT( err, relative_tolerance ) << "Element " << k << " failed. t_sim=(" << ts_x << ", " << ts_y << ", " << ts_z << ") t_exact=(" << te_x << ", " << te_y << ", " << te_z << ")";
      }
    } );

    volumeRegion.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase & subRegion )
    {
      auto const & avgStress = subRegion.getField< fields::solidMechanics::averageStress >();
      localIndex n_cells = subRegion.size();
      for( localIndex k=0; k < n_cells; ++k )
      {
        real64 sig_xx = avgStress( k, 0 );
        real64 sig_yy = avgStress( k, 1 );
        real64 sig_zz = avgStress( k, 2 );
        real64 sig_yz = avgStress( k, 3 );
        real64 sig_xz = avgStress( k, 4 );
        real64 sig_xy = avgStress( k, 5 );

        EXPECT_NEAR( std::fabs( sig_xx - s_xx ) / s_xx, 0.0, relative_tolerance ) << "Volume Element " << k << " failed xx stress.";
        EXPECT_NEAR( std::fabs( sig_yy - s_yy ) / s_yy, 0.0, relative_tolerance ) << "Volume Element " << k << " failed yy stress.";
        EXPECT_NEAR( std::fabs( sig_zz - s_zz ) / s_zz, 0.0, relative_tolerance ) << "Volume Element " << k << " failed zz stress.";
        EXPECT_NEAR( std::fabs( sig_yz ), 0.0, relative_tolerance ) << "Volume Element " << k << " failed yz stress.";
        EXPECT_NEAR( std::fabs( sig_xz ), 0.0, relative_tolerance ) << "Volume Element " << k << " failed xz stress.";
        EXPECT_NEAR( std::fabs( sig_xy ), 0.0, relative_tolerance ) << "Volume Element " << k << " failed xy stress.";
      }
    } );
  }
}

INSTANTIATE_TEST_SUITE_P(
  FractureStressCases,
  ConsistencyTest,
  ::testing::Values(
    // Mesh, s_xx, s_yy, s_zz
    // Hex meshes
    std::make_tuple( "fractured_mesh_hex_DFN_1.vtu", -1.0e7, -0.5e7, -2.0e7 ), // Triaxial
    std::make_tuple( "fractured_mesh_hex_DFN_2.vtu", -1.0e7, -0.5e7, -2.0e7 ), // Triaxial
    std::make_tuple( "fractured_mesh_hex_DFN_3.vtu", -1.0e7, -0.5e7, -2.0e7 ), // Triaxial
    std::make_tuple( "fractured_mesh_hex_DFN_12.vtu", -1.0e7, -0.5e7, -2.0e7 ), // Triaxial
    std::make_tuple( "fractured_mesh_hex_DFN_13.vtu", -1.0e7, -0.5e7, -2.0e7 ), // Triaxial
    std::make_tuple( "fractured_mesh_hex_DFN_123.vtu", -1.0e7, -0.5e7, -2.0e7 ), // Triaxial

    // Tet meshes
    std::make_tuple( "fractured_mesh_tet_DFN_1.vtu", -1.0e7, -0.5e7, -2.0e7 ), // Triaxial
    std::make_tuple( "fractured_mesh_tet_DFN_2.vtu", -1.0e7, -0.5e7, -2.0e7 ), // Triaxial
    std::make_tuple( "fractured_mesh_tet_DFN_3.vtu", -1.0e7, -0.5e7, -2.0e7 ), // Triaxial
    std::make_tuple( "fractured_mesh_tet_DFN_12.vtu", -1.0e7, -0.5e7, -2.0e7 ), // Triaxial
    std::make_tuple( "fractured_mesh_tet_DFN_13.vtu", -1.0e7, -0.5e7, -2.0e7 ), // Triaxial
    std::make_tuple( "fractured_mesh_tet_DFN_123.vtu", -1.0e7, -0.5e7, -2.0e7 )  // Triaxial
    )
  );

int main( int argc, char * argv[] )
{
  MPI_Init( &argc, &argv );
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv, false );
  int result = RUN_ALL_TESTS();
  geos::basicCleanup();
  MPI_Finalize();
  return result;
}
