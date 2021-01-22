/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "managers/initialization.hpp"
#include "managers/ProblemManager.hpp"
#include "virtualElement/VirtualElementBase.hpp"
#include "virtualElement/ConformingVirtualElementOrder1.hpp"
#include "managers/DomainPartition.hpp"
#include "meshUtilities/MeshManager.hpp"

// TPL includes
#include "gtest/gtest.h"

using namespace geosx;
using namespace virtualElement;

static void checkIntegralMeanConsistency( localIndex numBasisFunctions,
                                          array1d< real64 > const & basisFunctionsIntegralMean )
{
  real64 sum = 0;
  for( localIndex iBasisFun = 0; iBasisFun <  numBasisFunctions; ++iBasisFun )
  {
    sum += basisFunctionsIntegralMean( iBasisFun );
  }
  EXPECT_TRUE( abs( sum-1 ) < 1e-15 )
    << "Sum of basis functions integral mean is not 1, but " << sum << ". "
    << "The computed integral means are " << basisFunctionsIntegralMean;
}

static void
checkIntegralMeanDerivativesConsistency( localIndex numBasisFunctions,
                                         array2d< real64 > const & basisDerivativesIntegralMean )
{
  real64 sumX = 0, sumY = 0, sumZ = 0;
  for( localIndex iBasisFun = 0; iBasisFun <  numBasisFunctions; ++iBasisFun )
  {
    sumX += basisDerivativesIntegralMean( 0, iBasisFun );
    sumY += basisDerivativesIntegralMean( 1, iBasisFun );
    sumZ += basisDerivativesIntegralMean( 2, iBasisFun );
  }
  EXPECT_TRUE( abs( sumX ) < 1e-15 )
    << "Sum of the x-derivatives of basis functions integral mean is not 0, but " << sumX << ". "
    << "The computed integral means are " << basisDerivativesIntegralMean[0];
  EXPECT_TRUE( abs( sumY ) < 1e-15 )
    << "Sum of the y-derivatives of basis functions integral mean is not 0, but " << sumY << ". "
    << "The computed integral means are " << basisDerivativesIntegralMean[1];
  EXPECT_TRUE( abs( sumZ ) < 1e-15 )
    << "Sum of the z-derivatives of basis functions integral mean is not 0, but " << sumZ << ". "
    << "The computed integral means are " << basisDerivativesIntegralMean[2];
}

static void
checkStabilizationMatrixConsistency
  ( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoords,
  localIndex const & cellIndex,
  CellElementSubRegion::NodeMapType const & cellToNodes,
  arrayView2d< real64 const > const & cellCenters,
  array2d< real64 > const & stabilizationMatrix )
{
  localIndex const numCellPoints = cellToNodes[cellIndex].size();

  real64 cellDiameter = 0;
  for( localIndex numVertex = 0; numVertex < numCellPoints; ++numVertex )
  {
    for( localIndex numOthVertex = 0; numOthVertex < numVertex; ++numOthVertex )
    {
      array1d< real64 > vertDiff( 3 );
      LvArray::tensorOps::copy< 3 >( vertDiff, nodesCoords[cellToNodes( cellIndex, numVertex )] );
      LvArray::tensorOps::subtract< 3 >( vertDiff, nodesCoords[cellToNodes( cellIndex, numOthVertex )] );
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
                                          (nodesCoords( cellToNodes( cellIndex, numVertex ), pos ) - cellCenters( cellIndex, pos ));
  }

  array1d< real64 > stabTimeMonomialDofs( numCellPoints );
  real64 stabTimeMonomialDofsNorm = 0;
  for( localIndex i = 0; i < numCellPoints; ++i )
  {
    stabTimeMonomialDofs( i ) = 0;
    stabTimeMonomialDofsNorm = 0;
    for( localIndex j = 0; j < numCellPoints; ++j )
    {
      stabTimeMonomialDofs( i ) += stabilizationMatrix( i, j );
    }
    stabTimeMonomialDofsNorm += stabTimeMonomialDofs( i )*stabTimeMonomialDofs( i );
  }
  EXPECT_TRUE( abs( stabTimeMonomialDofsNorm ) < 1e-15 )
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
        stabTimeMonomialDofs( i ) += stabilizationMatrix( i, j )*monomialVemDofs( monomInd, j );
      }
      stabTimeMonomialDofsNorm += stabTimeMonomialDofs( i )*stabTimeMonomialDofs( i );
    }
    EXPECT_TRUE( abs( stabTimeMonomialDofsNorm ) < 1e-15 )
      << "Product of stabilization matrix and monomial degrees of freedom is not zero for "
      << "monomial number " << monomInd+1 << ". The computed product is " << stabTimeMonomialDofs;
  }
}

TEST( VirtualElementBase, unitCube )
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

  ProblemManager * problemManager = new ProblemManager( "Problem", nullptr );
  problemManager->initializePythonInterpreter();
  problemManager->processInputFileRecursive( xmlProblemNode );

  // Open mesh levels
  DomainPartition * domain  = problemManager->getDomainPartition();
  MeshManager * meshManager = problemManager->getGroup< MeshManager >
                                ( problemManager->groupKeys.meshManager );
  meshManager->generateMeshLevels( domain );
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * elementManager = mesh.getElemManager();
  xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child( elementManager->getName().c_str() );
  elementManager->processInputFileRecursive( topLevelNode );
  elementManager->postProcessInputRecursive();
  problemManager->problemSetup();

  ConformingVirtualElementOrder1 vemElement;
  vemElement.computeProjectors( mesh, 0, 0, 0 );

  checkIntegralMeanConsistency( vemElement.getNumSupportPoints(),
                                vemElement.m_basisFunctionsIntegralMean );
  checkIntegralMeanDerivativesConsistency( vemElement.getNumSupportPoints(),
                                           vemElement.m_basisDerivativesIntegralMean );

  NodeManager const & nodeManager = *mesh.getNodeManager();
  CellElementRegion const & cellRegion =
    *elementManager->getRegion< CellElementRegion >( 0 );
  CellElementSubRegion const & cellSubRegion =
    *cellRegion.getSubRegion< CellElementSubRegion >( 0 );
  CellElementSubRegion::NodeMapType const & cellToNodes = cellSubRegion.nodeList();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > nodesCoords =
    nodeManager.referencePosition();
  arrayView2d< real64 const > cellCenters = cellSubRegion.getElementCenter();
  localIndex const cellIndex = 0;
  checkStabilizationMatrixConsistency( nodesCoords, cellIndex, cellToNodes, cellCenters,
                                       vemElement.m_stabilizationMatrix );

  delete problemManager;
}

TEST( VirtualElementBase, wedges )
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

  ProblemManager * problemManager = new ProblemManager( "Problem", nullptr );
  problemManager->initializePythonInterpreter();
  problemManager->processInputFileRecursive( xmlProblemNode );

  // Open mesh levels
  DomainPartition * domain  = problemManager->getDomainPartition();
  MeshManager * meshManager = problemManager->getGroup< MeshManager >
                                ( problemManager->groupKeys.meshManager );
  meshManager->generateMeshLevels( domain );
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * elementManager = mesh.getElemManager();
  xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child( elementManager->getName().c_str() );
  elementManager->processInputFileRecursive( topLevelNode );
  elementManager->postProcessInputRecursive();
  problemManager->problemSetup();

  ConformingVirtualElementOrder1 vemElement;
  CellElementRegion const & cellRegion =
    *elementManager->getRegion< CellElementRegion >( 0 );
  CellElementSubRegion const & cellSubRegion =
    *cellRegion.getSubRegion< CellElementSubRegion >( 0 );
  NodeManager const & nodeManager = *mesh.getNodeManager();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > nodesCoords =
    nodeManager.referencePosition();
  CellElementSubRegion::NodeMapType const & cellToNodes = cellSubRegion.nodeList();
  arrayView2d< real64 const > cellCenters = cellSubRegion.getElementCenter();
  localIndex const numCells = cellSubRegion.getElementVolume().size();
  for( localIndex cellIndex = 0; cellIndex < numCells; ++cellIndex )
  {
    vemElement.computeProjectors( mesh, 0, 0, cellIndex );
    checkIntegralMeanConsistency( vemElement.getNumSupportPoints(),
                                  vemElement.m_basisFunctionsIntegralMean );
    checkIntegralMeanDerivativesConsistency( vemElement.getNumSupportPoints(),
                                             vemElement.m_basisDerivativesIntegralMean );


    checkStabilizationMatrixConsistency( nodesCoords, cellIndex, cellToNodes, cellCenters,
                                         vemElement.m_stabilizationMatrix );
  }
  delete problemManager;
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
