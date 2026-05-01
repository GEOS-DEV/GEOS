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

/**
 * @file testFEMConsistencySplitMeshFixture.hpp
 *
 * Shared test fixture and TEST_P body for the FEM consistency test using split (pre-fractured) meshes.
 * Both the serial and MPI test executables include this header and only
 * differ in the INSTANTIATE_TEST_SUITE_P parameters.
 *
 * The test verifies the consistency of the Finite Element Method (FEM)
 * implementation in GEOS for solid mechanics problems involving fractures
 * defined via faceBlocks in a VTM composite mesh.
 */

#ifndef GEOS_TEST_FEMCONSISTENCY_SPLITMESH_FIXTURE_HPP
#define GEOS_TEST_FEMCONSISTENCY_SPLITMESH_FIXTURE_HPP

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
#include "common/MpiWrapper.hpp"

using namespace geos;

constexpr real64 relative_tolerance_split = 1.0e-6; // Because MPa range is used

extern CommandLineOptions g_commandLineOptions;

/**
 * @brief This integration test verifies the consistency of the Finite Element Method (FEM) implementation
 *        in GEOS for solid mechanics problems involving pre-split fractures defined via faceBlocks in a VTM mesh.
 *
 * The test performs the following steps:
 * 1. Generates an XML input file for a given split mesh and stress boundary conditions.
 * 2. Runs the simulation using the SolidMechanicsAugmentedLagrangianContact solver.
 * 3. Verifies that the computed stress on fracture faces matches the expected traction derived from
 *    the applied stress field.
 * 4. Verifies that the stress in the volume elements matches the applied stress boundary conditions.
 *
 * The test is parameterized to run on various mesh types (hex, tet) with pre-defined fracture face blocks.
 *
 * Parametrized with std::tuple<std::string, real64, real64, real64, std::tuple<int, int, int>>:
 * - std::string: Mesh file name (e.g., "DFN_5_fractures_hex_binarized.vtm").
 * - real64: Applied stress component XX (s_xx).
 * - real64: Applied stress component YY (s_yy).
 * - real64: Applied stress component ZZ (s_zz).
 * - tuple<int, int, int>: Number of partitions in the x, y, z directions.
 */
class ConsistencySplitMeshTest : public ::testing::TestWithParam< std::tuple< std::string, real64, real64, real64, std::tuple< int, int, int > > >
{
protected:
  void SetUp() override
  {
    testBinaryDir = TEST_BINARY_DIR;
  }

  std::string generateXmlInput( std::string const & meshFile, real64 s_xx, real64 s_yy, real64 s_zz )
  {
    std::ostringstream oss;
    oss << R"xml(<?xml version="1.0" ?>
<Problem>
  <Mesh>
    <VTKMesh name="mesh1" checkEulerCharacteristic="1" useGlobalIds="1" file=")xml" << meshFile << R"xml(" faceBlocks="{ Fault }"/>
  </Mesh>

  <Geometry>
    <Box name="xnegFace" xMin="{ -0.01, -0.01, -0.01 }" xMax="{  0.01,  1.01,  1.01 }"/>
    <Box name="xposFace" xMin="{  0.99, -0.01, -0.01 }" xMax="{  1.01,  1.01,  1.01 }"/>
    <Box name="ynegFace" xMin="{ -0.01, -0.01, -0.01 }" xMax="{  1.01,  0.01,  1.01 }"/>
    <Box name="yposFace" xMin="{ -0.01,  0.99, -0.01 }" xMax="{  1.01,  1.01,  1.01 }"/>
    <Box name="znegFace" xMin="{ -0.01, -0.01, -0.01 }" xMax="{  1.01,  1.01,  0.01 }"/>
    <Box name="zposFace" xMin="{ -0.01, -0.01,  0.99 }" xMax="{  1.01,  1.01,  1.01 }"/>
  </Geometry>

  <Solvers gravityVector="{0.0, 0.0, 0.0}">
    <SolidMechanicsAugmentedLagrangianContact 
      name="mechSolver"
      simultaneous="1"
      symmetric="1"
      iterPenaltyN="1.0e1"
      iterPenaltyT="1.0e-1"
      tolNormalTrac="1.e-04"
      tolTauLimit="1.e-04"
      tolJumpN="1.e-4"
      tolJumpT="1.e-4" 
      discretization="FE1" 
      targetRegions="{ Region, Fracture }" 
      logLevel="1">
      <NonlinearSolverParameters newtonTol="1.0e-7" newtonMaxIter="20" logLevel="1"/>
      <LinearSolverParameters directParallel="0"/>
    </SolidMechanicsAugmentedLagrangianContact>
  </Solvers>

  <NumericalMethods>
    <FiniteElements>
      <FiniteElementSpace name="FE1" order="1" useHighOrderQuadratureRule="1" />
    </FiniteElements>
  </NumericalMethods>

  <ElementRegions>
    <CellElementRegion name="Region" cellBlocks="{ * }" materialList="{ rock }"/>
    <SurfaceElementRegion name="Fracture" 
                          defaultAperture="1.0e-4" 
                          faceBlock="Fault" 
                          materialList="{ fractureContact }"/>
  </ElementRegions>

  <Constitutive>
    <ElasticIsotropic name="rock" 
                      defaultDensity="2500" 
                      defaultYoungModulus="1.0e9" 
                      defaultPoissonRatio="0.25"/>
    <Coulomb name="fractureContact" 
             cohesion="1.0e10" 
             frictionCoefficient="0.5"/>
  </Constitutive>

  <FieldSpecifications>
    <FieldSpecification name="xneg_disp" component="0" setNames="{ xnegFace }" 
                        objectPath="nodeManager" fieldName="totalDisplacement" scale="0.0"/>
    <FieldSpecification name="yneg_disp" component="1" setNames="{ ynegFace }" 
                        objectPath="nodeManager" fieldName="totalDisplacement" scale="0.0"/>
    <FieldSpecification name="zneg_disp" component="2" setNames="{ znegFace }" 
                        objectPath="nodeManager" fieldName="totalDisplacement" scale="0.0"/>
    
    <Traction name="xpos_traction" setNames="{ xposFace }" objectPath="faceManager" 
              tractionType="normal" scale=")xml" << s_xx << R"xml("/>
    <Traction name="ypos_traction" setNames="{ yposFace }" objectPath="faceManager" 
              tractionType="normal" scale=")xml" << s_yy << R"xml("/>
    <Traction name="zpos_traction" setNames="{ zposFace }" objectPath="faceManager" 
              tractionType="normal" scale=")xml" << s_zz << R"xml("/>
  </FieldSpecifications>

  <Tasks>
    <SolidMechanicsAugmentedLagrangianContactInitialization name="ELASTICITY.POST.SPLIT.STEP" 
                                                            solidSolverName="mechSolver" logLevel="0" reset_displacements="0"/>
  </Tasks>  
  <Outputs>
    <VTK name="vtkOutputM" outputRegionType="cell" />
    <VTK name="vtkOutputF" outputRegionType="surface" />
  </Outputs>
  <Events minTime="0.0" maxTime="1.0">
    <SoloEvent name="ELASTICITY.POST.SPLIT.STEP" targetTime="0.0" beginTime="0.0" 
               target="/Tasks/ELASTICITY.POST.SPLIT.STEP"/>
    <PeriodicEvent name="outputsF" timeFrequency="1.0" target="/Outputs/vtkOutputF"/>
    <PeriodicEvent name="outputsM" timeFrequency="1.0" target="/Outputs/vtkOutputM"/>
    <PeriodicEvent name="solverApplications" target="/Solvers/mechSolver" forceDt="1.0"/>
  </Events>
</Problem>
)xml";
    return oss.str();
  }

  std::string testBinaryDir;
};

TEST_P( ConsistencySplitMeshTest, Run )
{
  auto const & params = GetParam();
  auto const & meshFileName = std::get< 0 >( params );
  real64 const s_xx = std::get< 1 >( params );
  real64 const s_yy = std::get< 2 >( params );
  real64 const s_zz = std::get< 3 >( params );
  std::tuple< int, int, int > const & partitions = std::get< 4 >( params );
  int const xPartitions = std::get< 0 >( partitions );
  int const yPartitions = std::get< 1 >( partitions );
  int const zPartitions = std::get< 2 >( partitions );
  std::string xmlPath = testBinaryDir + "/test_fem_consistency_splitmesh_" + meshFileName + ".xml";

  std::string xmlContent = generateXmlInput( meshFileName, s_xx, s_yy, s_zz );

  if( MpiWrapper::commRank() == 0 )
  {
    std::ofstream ofs( xmlPath );
    ofs << xmlContent;
  }
  MpiWrapper::barrier();

  std::unique_ptr< CommandLineOptions > options = std::make_unique< CommandLineOptions >( g_commandLineOptions );
  options->inputFileNames.push_back( xmlPath );
  options->problemName = "test_fem_consistency_splitmesh";

  options->xPartitionsOverride = xPartitions;
  options->yPartitionsOverride = yPartitions;
  options->zPartitionsOverride = zPartitions;
  options->overridePartitionNumbers = true;

  // Scoped state to ensure cleanup
  {
    GeosxState state( std::move( options ) );
    ASSERT_TRUE( state.initializeDataRepository() );
    state.applyInitialConditions();
    state.run();

    std::cout << "Finished running test. Now verifying results..." << std::endl;

    // Verification using bulk stress and geometric normal
    ProblemManager & pm = state.getProblemManager();
    SurfaceElementRegion & fractureRegion = pm.getDomainPartition().getMeshBody( "mesh1" ).getBaseDiscretization().getElemManager().template getRegion< SurfaceElementRegion >( "Fracture" );
    auto & elemManager = pm.getDomainPartition().getMeshBody( "mesh1" ).getBaseDiscretization().getElemManager();
    CellElementRegion & volumeRegion = elemManager.template getRegion< CellElementRegion >( "Region" );

    FaceManager const & faceManager = pm.getDomainPartition().getMeshBody( "mesh1" ).getBaseDiscretization().getFaceManager();
    auto const & nodeManager = pm.getDomainPartition().getMeshBody( "mesh1" ).getBaseDiscretization().getNodeManager();

    auto const & faceToCell = faceManager.elementList();
    auto const & faceToRegion = faceManager.elementRegionList();
    auto const & faceToSubRegion = faceManager.elementSubRegionList();
    auto const & faceNodes = faceManager.nodeList();
    auto const & nodePos = nodeManager.referencePosition();

    // Split meshes are expected to have conforming fracture geometry, so use standard tolerance
    real64 const shear_tolerance = relative_tolerance_split;
    real64 const normal_tolerance = relative_tolerance_split;

    fractureRegion.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      auto const & fractureTraction = subRegion.getField< fields::contact::traction >();
      auto const & rotationMatrix = subRegion.getField< fields::contact::rotationMatrix >();
      auto const & faces = subRegion.faceList();
      localIndex const numFractureElements = subRegion.size();

      // ===================================================================================
      // Test 1: Verify traction from volume element stress matches expected applied stress
      // This iterates over all faces adjacent to fracture elements (2 per fracture element)
      // ===================================================================================
      localIndex const numFaces = faces.size( 1 );  // Number of faces per fracture element (typically 2)
      for( localIndex k = 0; k < numFractureElements; ++k )
      {
        for( localIndex faceIdx = 0; faceIdx < numFaces; ++faceIdx )
        {
          localIndex f = faces( k, faceIdx );

          // Compute Geometric Normal using Dziuk method
          auto const & fnodes = faceNodes[f];
          if( fnodes.size() < 3 )
            continue;

          real64 nx = 0.0;
          real64 ny = 0.0;
          real64 nz = 0.0;

          if( fnodes.size() == 4 )
          {
            // For quadrilaterals: Use diagonal-based triangulation (more stable for warped quads)
            real64 p0[3] = { nodePos[fnodes[0]][0], nodePos[fnodes[0]][1], nodePos[fnodes[0]][2] };
            real64 p1[3] = { nodePos[fnodes[1]][0], nodePos[fnodes[1]][1], nodePos[fnodes[1]][2] };
            real64 p2[3] = { nodePos[fnodes[2]][0], nodePos[fnodes[2]][1], nodePos[fnodes[2]][2] };
            real64 p3[3] = { nodePos[fnodes[3]][0], nodePos[fnodes[3]][1], nodePos[fnodes[3]][2] };

            real64 v01[3] = { p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2] };
            real64 v02[3] = { p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2] };
            real64 cp1_x = v01[1]*v02[2] - v01[2]*v02[1];
            real64 cp1_y = v01[2]*v02[0] - v01[0]*v02[2];
            real64 cp1_z = v01[0]*v02[1] - v01[1]*v02[0];

            // Triangle 2: (V0, V2, V3)
            real64 v03[3] = { p3[0]-p0[0], p3[1]-p0[1], p3[2]-p0[2] };
            real64 cp2_x = v02[1]*v03[2] - v02[2]*v03[1];
            real64 cp2_y = v02[2]*v03[0] - v02[0]*v03[2];
            real64 cp2_z = v02[0]*v03[1] - v02[1]*v03[0];

            // Accumulate cross products
            nx = cp1_x + cp2_x;
            ny = cp1_y + cp2_y;
            nz = cp1_z + cp2_z;
          }
          else
          {
            // For non-quad polygons: Use centroid-based fan triangulation
            real64 centroid[3] = { 0.0, 0.0, 0.0 };
            for( localIndex i = 0; i < fnodes.size(); ++i )
            {
              centroid[0] += nodePos[fnodes[i]][0];
              centroid[1] += nodePos[fnodes[i]][1];
              centroid[2] += nodePos[fnodes[i]][2];
            }
            centroid[0] /= fnodes.size();
            centroid[1] /= fnodes.size();
            centroid[2] /= fnodes.size();

            for( localIndex i = 0; i < fnodes.size(); ++i )
            {
              localIndex const next = ( i + 1 ) % fnodes.size();

              real64 v1[3] = { nodePos[fnodes[i]][0] - centroid[0],
                               nodePos[fnodes[i]][1] - centroid[1],
                               nodePos[fnodes[i]][2] - centroid[2] };
              real64 v2[3] = { nodePos[fnodes[next]][0] - centroid[0],
                               nodePos[fnodes[next]][1] - centroid[1],
                               nodePos[fnodes[next]][2] - centroid[2] };

              real64 cp_x = v1[1]*v2[2] - v1[2]*v2[1];
              real64 cp_y = v1[2]*v2[0] - v1[0]*v2[2];
              real64 cp_z = v1[0]*v2[1] - v1[1]*v2[0];

              nx += cp_x;
              ny += cp_y;
              nz += cp_z;
            }
          }

          real64 norm = std::sqrt( nx*nx + ny*ny + nz*nz );
          nx /= norm; ny /= norm; nz /= norm;

          // Get neighbor cell - try to find a valid one from either side of the face
          localIndex c = -1;
          localIndex er = -1;
          localIndex esr = -1;

          localIndex c0 = faceToCell[f][0];
          localIndex er0 = faceToRegion[f][0];
          localIndex esr0 = faceToSubRegion[f][0];

          localIndex c1 = faceToCell[f][1];
          localIndex er1 = faceToRegion[f][1];
          localIndex esr1 = faceToSubRegion[f][1];

          // Use the first valid neighbor (one with valid cell, region, and subregion indices)
          if( c0 != static_cast< localIndex >(-1) &&
              er0 != static_cast< localIndex >(-1) &&
              esr0 != static_cast< localIndex >(-1) )
          {
            c = c0;
            er = er0;
            esr = esr0;
          }
          else if( c1 != static_cast< localIndex >(-1) &&
                   er1 != static_cast< localIndex >(-1) &&
                   esr1 != static_cast< localIndex >(-1) )
          {
            c = c1;
            er = er1;
            esr = esr1;
          }
          else
          {
            // No valid neighbor cell found on this rank, skip this face
            continue;
          }

          ElementRegionBase & region = elemManager.getRegion( er );
          ElementSubRegionBase & cellSubRegion = region.getSubRegion( esr );
          arrayView2d< real64 const, cells::RANK2_TENSOR_USD > const avgStress = cellSubRegion.getField< fields::solidMechanics::averageStress >();
          avgStress.move( LvArray::MemorySpace::host );

          real64 sig_xx = avgStress( c, 0 );
          real64 sig_yy = avgStress( c, 1 );
          real64 sig_zz = avgStress( c, 2 );
          real64 sig_yz = avgStress( c, 3 );
          real64 sig_xz = avgStress( c, 4 );
          real64 sig_xy = avgStress( c, 5 );

          // Compute t_sim = sigma * n
          real64 ts_x = sig_xx * nx + sig_xy * ny + sig_xz * nz;
          real64 ts_y = sig_xy * nx + sig_yy * ny + sig_yz * nz;
          real64 ts_z = sig_xz * nx + sig_yz * ny + sig_zz * nz;

          // Compute t_exact = S_input * n
          real64 te_x = s_xx * nx;
          real64 te_y = s_yy * ny;
          real64 te_z = s_zz * nz;

          real64 const err_abs = std::sqrt( std::pow( ts_x - te_x, 2 ) + std::pow( ts_y - te_y, 2 ) + std::pow( ts_z - te_z, 2 ) );
          real64 const norm_te = std::sqrt( std::pow( te_x, 2 ) + std::pow( te_y, 2 ) + std::pow( te_z, 2 ) );
          real64 const err = ( norm_te > 1.0e-16 ) ? err_abs / norm_te : err_abs;
          EXPECT_LT( err, normal_tolerance ) << "Fracture element " << k << " face " << faceIdx
                                             << " failed. t_sim=(" << ts_x << ", " << ts_y << ", " << ts_z
                                             << ") t_exact=(" << te_x << ", " << te_y << ", " << te_z << ")";
        }
      }

      // ===================================================================================
      // Test 2: Verify fracture traction field matches expected traction
      // This iterates only over fracture elements (not faces)
      // ===================================================================================
      for( localIndex k = 0; k < numFractureElements; ++k )
      {
        // Get the fracture element's normal from the rotation matrix
        // The rotation matrix columns are: [normal, tangent1, tangent2]
        real64 nx = rotationMatrix( k, 0, 0 );
        real64 ny = rotationMatrix( k, 1, 0 );
        real64 nz = rotationMatrix( k, 2, 0 );

        // Compute t_exact = S_input * n (in global coordinates)
        real64 te_x = s_xx * nx;
        real64 te_y = s_yy * ny;
        real64 te_z = s_zz * nz;

        // Project exact traction onto normal (should match tf_n)
        real64 te_n = te_x * nx + te_y * ny + te_z * nz;

        // Project exact traction onto tangential plane
        real64 te_tang_x = te_x - te_n * nx;
        real64 te_tang_y = te_y - te_n * ny;
        real64 te_tang_z = te_z - te_n * nz;
        real64 te_tang_mag = std::sqrt( te_tang_x*te_tang_x + te_tang_y*te_tang_y + te_tang_z*te_tang_z );

        // Get fracture traction from contact solver (in local coordinates: normal, tangent1, tangent2)
        real64 tf_n = fractureTraction( k, 0 );
        real64 tf_t1 = fractureTraction( k, 1 );
        real64 tf_t2 = fractureTraction( k, 2 );
        real64 tf_tang_mag = std::sqrt( tf_t1*tf_t1 + tf_t2*tf_t2 );

        real64 const norm_te = std::sqrt( te_x*te_x + te_y*te_y + te_z*te_z );
        real64 const err_n = std::fabs( tf_n - te_n );
        real64 const err_t = std::fabs( tf_tang_mag - te_tang_mag );

        real64 const rel_err_n = ( norm_te > 1.0e-16 ) ? err_n / norm_te : err_n;
        real64 const rel_err_t = ( norm_te > 1.0e-16 ) ? err_t / norm_te : err_t;

        EXPECT_LT( rel_err_n, normal_tolerance ) << "Fracture element " << k << " failed normal traction check. tf_n=" << tf_n << " te_n=" << te_n;
        EXPECT_LT( rel_err_t, shear_tolerance ) << "Fracture element " << k << " failed tangential traction check. tf_t=" << tf_tang_mag << " te_t=" << te_tang_mag;
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

        EXPECT_NEAR( std::fabs( (sig_xx - s_xx ) / s_xx ), 0.0, normal_tolerance ) << "Volume Element " << k << " failed xx stress.";
        EXPECT_NEAR( std::fabs( (sig_yy - s_yy ) / s_yy ), 0.0, normal_tolerance ) << "Volume Element " << k << " failed yy stress.";
        EXPECT_NEAR( std::fabs( (sig_zz - s_zz ) / s_zz ), 0.0, normal_tolerance ) << "Volume Element " << k << " failed zz stress.";
        EXPECT_NEAR( std::fabs( sig_yz ) / 1.0e6, 0.0, shear_tolerance ) << "Volume Element " << k << " failed yz stress.";
        EXPECT_NEAR( std::fabs( sig_xz ) / 1.0e6, 0.0, shear_tolerance ) << "Volume Element " << k << " failed xz stress.";
        EXPECT_NEAR( std::fabs( sig_xy ) / 1.0e6, 0.0, shear_tolerance ) << "Volume Element " << k << " failed xy stress.";
      }
    } );
  }
}

#endif /* GEOS_TEST_FEMCONSISTENCY_SPLITMESH_FIXTURE_HPP */
