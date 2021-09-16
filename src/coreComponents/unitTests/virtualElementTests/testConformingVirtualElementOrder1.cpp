/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "virtualElement/ConformingVirtualElementOrder1.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/MeshManager.hpp"

// TPL includes
#include "gtest/gtest.h"

using namespace geosx;
using namespace finiteElement;

CommandLineOptions g_commandLineOptions;

template< typename VEM >
GEOSX_HOST_DEVICE
static void checkIntegralMeanConsistency( FiniteElementBase const & feBase,
                                          typename VEM::StackVariables const & stack )
{
  static constexpr localIndex
    maxSupportPoints = FiniteElementBase::getMaxSupportPoints< VEM >();
  real64 basisFunctionsIntegralMean[maxSupportPoints];
  VEM::calcN( 0, stack, basisFunctionsIntegralMean );
  real64 sum = 0;
  for( localIndex iBasisFun = 0;
       iBasisFun < feBase.template numSupportPoints< VEM >( stack ); ++iBasisFun )
  {
    sum += basisFunctionsIntegralMean[iBasisFun];
  }
  EXPECT_TRUE( LvArray::math::abs( sum-1 ) < 1e-15 )
    << "Sum of basis functions integral mean is not 1, but " << sum << ". "
    << "The computed integral means are " << basisFunctionsIntegralMean;
}

template< typename VEM >
GEOSX_HOST_DEVICE
static void
checkIntegralMeanDerivativesConsistency( FiniteElementBase const & feBase,
                                         typename VEM::StackVariables const & stack )
{
  static constexpr localIndex
    maxSupportPoints = FiniteElementBase::getMaxSupportPoints< VEM >();
  real64 const dummy[VEM::numNodes][3] { { 0.0 } };
  localIndex const k = 0;
  for( localIndex q = 0; q < VEM::numQuadraturePoints; ++q )
  {
    real64 basisDerivativesIntegralMean[maxSupportPoints][3];
    feBase.template getGradN< VEM >( k, q, dummy, stack, basisDerivativesIntegralMean );
    real64 sumX = 0, sumY = 0, sumZ = 0;
    for( localIndex iBasisFun = 0;
         iBasisFun < feBase.template numSupportPoints< VEM >( stack );
         ++iBasisFun )
    {
      sumX += basisDerivativesIntegralMean[iBasisFun][0];
      sumY += basisDerivativesIntegralMean[iBasisFun][1];
      sumZ += basisDerivativesIntegralMean[iBasisFun][2];
    }
    EXPECT_TRUE( LvArray::math::abs( sumX ) < 1e-15 )
      << "Sum of the x-derivatives of basis functions integral mean is not 0, but " << sumX << ". "
      << "The computed integral means are " << basisDerivativesIntegralMean;
    EXPECT_TRUE( LvArray::math::abs( sumY ) < 1e-15 )
      << "Sum of the y-derivatives of basis functions integral mean is not 0, but " << sumY << ". "
      << "The computed integral means are " << basisDerivativesIntegralMean;
    EXPECT_TRUE( LvArray::math::abs( sumZ ) < 1e-15 )
      << "Sum of the z-derivatives of basis functions integral mean is not 0, but " << sumZ << ". "
      << "The computed integral means are " << basisDerivativesIntegralMean;
  }
}

template< typename VEM >
GEOSX_HOST_DEVICE
static void
checkStabilizationMatrixConsistency ( arrayView2d< real64 const,
                                                   nodes::REFERENCE_POSITION_USD > const & nodesCoords,
                                      localIndex const & cellIndex,
                                      CellElementSubRegion::NodeMapType const & cellToNodes,
                                      arrayView2d< real64 const > const & cellCenters,
                                      FiniteElementBase const & feBase,
                                      typename VEM::StackVariables const & stack )
{
  static constexpr localIndex
    maxSupportPoints = FiniteElementBase::getMaxSupportPoints< VEM >();
  localIndex const numCellPoints = cellToNodes[cellIndex].size();

  real64 cellDiameter = 0;
  for( localIndex numVertex = 0; numVertex < numCellPoints; ++numVertex )
  {
    for( localIndex numOthVertex = 0; numOthVertex < numVertex; ++numOthVertex )
    {
      array1d< real64 > vertDiff( 3 );
      LvArray::tensorOps::copy< 3 >( vertDiff, nodesCoords[cellToNodes( cellIndex, numVertex )] );
      LvArray::tensorOps::subtract< 3 >( vertDiff,
                                         nodesCoords[cellToNodes( cellIndex, numOthVertex )] );
      real64 const candidateDiameter = LvArray::tensorOps::l2NormSquared< 3 >( vertDiff );
      if( cellDiameter < candidateDiameter )
        cellDiameter = candidateDiameter;
    }
  }
  cellDiameter = LvArray::math::sqrt< real64 >( cellDiameter );
  real64 const invCellDiameter = 1.0/cellDiameter;

  array2d< real64 > monomialVemDofs( 3, numCellPoints );
  for( localIndex numVertex = 0; numVertex < numCellPoints; ++numVertex )
  {
    for( localIndex pos = 0; pos < 3; ++pos )
      monomialVemDofs( pos, numVertex ) = invCellDiameter*
                                          (nodesCoords( cellToNodes( cellIndex, numVertex ), pos )
                                           - cellCenters( cellIndex, pos ));
  }

  array1d< real64 > stabTimeMonomialDofs( numCellPoints );
  real64 stabilizationMatrix[maxSupportPoints][maxSupportPoints] { { 0.0 } };
  feBase.template addGradGradStabilizationMatrix< VEM >( stack, stabilizationMatrix );
  real64 stabTimeMonomialDofsNorm = 0;
  for( localIndex i = 0; i < numCellPoints; ++i )
  {
    stabTimeMonomialDofs( i ) = 0;
    stabTimeMonomialDofsNorm = 0;
    for( localIndex j = 0; j < numCellPoints; ++j )
    {
      stabTimeMonomialDofs( i ) += stabilizationMatrix[ i ][ j ];
    }
    stabTimeMonomialDofsNorm += stabTimeMonomialDofs( i ) * stabTimeMonomialDofs( i );
  }
  EXPECT_TRUE( LvArray::math::abs( stabTimeMonomialDofsNorm ) < 1e-15 )
    << "Product of stabilization matrix and monomial degrees of freedom is not zero for "
    << "monomial number 0. The computed product is " << stabTimeMonomialDofs;
  for( localIndex monomInd = 0; monomInd < 3; ++monomInd )
  {
    stabTimeMonomialDofsNorm = 0;
    for( localIndex i = 0; i < numCellPoints; ++i )
    {
      stabTimeMonomialDofs( i ) = 0;
      for( localIndex j = 0; j < numCellPoints; ++j )
      {
        stabTimeMonomialDofs( i ) += stabilizationMatrix[ i ][ j ] * monomialVemDofs( monomInd, j );
      }
      stabTimeMonomialDofsNorm += stabTimeMonomialDofs( i ) * stabTimeMonomialDofs( i );
    }
    EXPECT_TRUE( LvArray::math::abs( stabTimeMonomialDofsNorm ) < 1e-15 )
      << "Product of stabilization matrix and monomial degrees of freedom is not zero for "
      << "monomial number " << monomInd+1 << ". The computed product is " << stabTimeMonomialDofs;
  }
}

template< typename VEM >
GEOSX_HOST_DEVICE
static void checkSumOfQuadratureWeights( real64 const & cellVolume,
                                         typename VEM::StackVariables stack )
{
  static constexpr localIndex
    maxSupportPoints = FiniteElementBase::getMaxSupportPoints< VEM >();
  real64 sum = 0.0;
  real64 const dummy[maxSupportPoints][3] { { 0.0 } };
  for( localIndex q = 0; q < VEM::numQuadraturePoints; ++q )
  {
    real64 weight = VEM::transformedQuadratureWeight( q, dummy, stack );
    sum += weight;
  }
  EXPECT_TRUE( LvArray::math::abs( sum - cellVolume ) < 1e-15 )
    << "Sum of quadrature weights does not equal the cell volume. Sum is " << sum
    << ". Cell volume is " << cellVolume;
}

template< localIndex MAXCELLNODES, localIndex MAXFACENODES >
static void testCellsInMeshLevel( MeshLevel const & mesh )
{
  // Get managers.
  ElementRegionManager const & elementManager = mesh.getElemManager();
  CellElementRegion const & cellRegion =
    elementManager.getRegion< CellElementRegion >( 0 );
  CellElementSubRegion const & cellSubRegion =
    cellRegion.getSubRegion< CellElementSubRegion >( 0 );
  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager const & faceManager = mesh.getFaceManager();
  EdgeManager const & edgeManager = mesh.getEdgeManager();

  // Get geometric properties to be passed as inputs.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > nodesCoords =
    nodeManager.referencePosition();
  CellBlock::NodeMapType const & cellToNodeMap = cellSubRegion.nodeList();
  arrayView2d< real64 const > cellCenters = cellSubRegion.getElementCenter();
  arrayView1d< real64 const > cellVolumes = cellSubRegion.getElementVolume();
  // Loop over cells.
  localIndex const numCells = cellSubRegion.getElementVolume().size();

  using VEM = ConformingVirtualElementOrder1< MAXCELLNODES, MAXFACENODES >;
  typename VEM::Initialization initialization;
  FiniteElementBase::initialize< VEM >( nodeManager, edgeManager,
                                        faceManager, cellSubRegion,
                                        initialization );

  forAll< parallelDevicePolicy< > >( numCells, [=] GEOSX_HOST_DEVICE
                                       ( localIndex const cellIndex )
  {
    typename VEM::StackVariables stack;
    VEM virtualElement;
    virtualElement.template setup< VEM >( cellIndex, initialization, stack );

    checkIntegralMeanConsistency< VEM >( virtualElement, stack );
    checkIntegralMeanDerivativesConsistency< VEM >( virtualElement, stack );
    checkStabilizationMatrixConsistency< VEM >( nodesCoords, cellIndex, cellToNodeMap,
                                                cellCenters, virtualElement, stack );
    checkSumOfQuadratureWeights< VEM >( cellVolumes[cellIndex], stack );
  } );
}

TEST( ConformingVirtualElementOrder1, hexagons )
{
  string const inputStream=
    "<Problem>"
    "  <Mesh>"
    "    <InternalMesh"
    "      name=\"cube\""
    "      elementTypes=\"{C3D8}\""
    "      xCoords=\"{0.0, 1.0}\""
    "      yCoords=\"{0.0, 1.0}\""
    "      zCoords=\"{0.0, 1.0}\""
    "      nx=\"{1}\""
    "      ny=\"{1}\""
    "      nz=\"{1}\""
    "      cellBlockNames=\"{cb1}\""
    "    />"
    "  </Mesh>"
    "  <ElementRegions>"
    "    <CellElementRegion name=\"region1\" cellBlocks=\"{cb1}\""
    "                       materialList=\"{dummy_material}\" />"
    "  </ElementRegions>"
    "</Problem>";

  xmlWrapper::xmlDocument inputFile;
  xmlWrapper::xmlResult xmlResult = inputFile.load_buffer( inputStream.c_str(), inputStream.size());
  if( !xmlResult )
  {
    GEOSX_LOG_RANK_0( "XML parsed with errors!" );
    GEOSX_LOG_RANK_0( "Error description: " << xmlResult.description());
    GEOSX_LOG_RANK_0( "Error offset: " << xmlResult.offset );
  }
  xmlWrapper::xmlNode xmlProblemNode = inputFile.child( "Problem" );

  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );

  ProblemManager & problemManager = state.getProblemManager();
  problemManager.processInputFileRecursive( xmlProblemNode );

  // Open mesh levels
  DomainPartition & domain  = problemManager.getDomainPartition();
  MeshManager & meshManager = problemManager.getGroup< MeshManager >( problemManager.groupKeys
                                                                        .meshManager );
  meshManager.generateMeshLevels( domain );
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager & elementManager = mesh.getElemManager();
  xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child( elementManager.getName().c_str() );
  elementManager.processInputFileRecursive( topLevelNode );
  elementManager.postProcessInputRecursive();
  problemManager.problemSetup();

  // Test computed projectors for all cells in MeshLevel
  testCellsInMeshLevel< 10, 6 >( mesh );
}

TEST( ConformingVirtualElementOrder1, wedges )
{
  string const inputStream=
    "<Problem>"
    "  <Mesh>"
    "    <InternalMesh"
    "      name=\"wedges\""
    "      elementTypes=\"{C3D6}\""
    "      xCoords=\"{0.0, 1.0}\""
    "      yCoords=\"{0.0, 1.0}\""
    "      zCoords=\"{0.0, 1.0}\""
    "      nx=\"{1}\""
    "      ny=\"{1}\""
    "      nz=\"{1}\""
    "      cellBlockNames=\"{cb1}\""
    "    />"
    "  </Mesh>"
    "  <ElementRegions>"
    "    <CellElementRegion name=\"region1\" cellBlocks=\"{cb1}\""
    "                       materialList=\"{dummy_material}\" />"
    "  </ElementRegions>"
    "</Problem>";
  xmlWrapper::xmlDocument inputFile;
  xmlWrapper::xmlResult xmlResult = inputFile.load_buffer( inputStream.c_str(), inputStream.size());
  if( !xmlResult )
  {
    GEOSX_LOG_RANK_0( "XML parsed with errors!" );
    GEOSX_LOG_RANK_0( "Error description: " << xmlResult.description());
    GEOSX_LOG_RANK_0( "Error offset: " << xmlResult.offset );
  }
  xmlWrapper::xmlNode xmlProblemNode = inputFile.child( "Problem" );

  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );

  ProblemManager & problemManager = state.getProblemManager();
  problemManager.processInputFileRecursive( xmlProblemNode );

  // Open mesh levels
  DomainPartition & domain  = problemManager.getDomainPartition();
  MeshManager & meshManager = problemManager.getGroup< MeshManager >
                                ( problemManager.groupKeys.meshManager );
  meshManager.generateMeshLevels( domain );
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager & elementManager = mesh.getElemManager();
  xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child( elementManager.getName().c_str() );
  elementManager.processInputFileRecursive( topLevelNode );
  elementManager.postProcessInputRecursive();
  problemManager.problemSetup();

  // Test computed projectors for all cells in MeshLevel
  testCellsInMeshLevel< 8, 9 >( mesh );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geosx::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
