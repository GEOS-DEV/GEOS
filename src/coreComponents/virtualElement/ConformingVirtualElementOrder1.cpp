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
/**
 * @file ConformingVirtualElementOrder1.cpp
 */

#include "ConformingVirtualElementOrder1.hpp"

namespace geosx
{
namespace virtualElement
{
void ConformingVirtualElementOrder1::ComputeProjectors( MeshLevel const & mesh,
                                                        localIndex const & regionIndex,
                                                        localIndex const & subRegionIndex,
                                                        localIndex const & cellIndex )
{
  // Get geometry managers
  NodeManager const & nodeManager = *mesh.getNodeManager();
  FaceManager const & faceManager = *mesh.getFaceManager();
  ElementRegionManager const & elementManager = *mesh.getElemManager();

  // Get pre-computed maps
  CellElementRegion const & cellRegion =
    *elementManager.GetRegion< CellElementRegion >( regionIndex );
  CellElementSubRegion const & cellSubRegion =
    *cellRegion.GetSubRegion< CellElementSubRegion >( subRegionIndex );
  CellElementSubRegion::NodeMapType const & cellToNodes = cellSubRegion.nodeList();
  FixedOneToManyRelation const & elementToFaceMap = cellSubRegion.faceList();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > nodesCoords =
    nodeManager.referencePosition();

  // Get pre-computed geometrical properties.
  arrayView2d< real64 const > faceCenters = faceManager.faceCenter();
  arrayView2d< real64 const > faceNormals = faceManager.faceNormal();
  arrayView2d< real64 const > cellCenters = cellSubRegion.getElementCenter();
  arrayView1d< real64 const > cellVolumes = cellSubRegion.getElementVolume();
  arraySlice1d< real64 const > cellCenter = cellCenters[cellIndex];
  localIndex const numCellFaces = elementToFaceMap[cellIndex].size();
  localIndex const numCellPoints = cellToNodes[cellIndex].size();
  m_numSupportPoints = numCellPoints;

  // Compute other geometrical properties.
  //  - compute map used to locate local point position by global id (used in the computation of
  //    basis functions boundary integrals.
  map< localIndex, localIndex > cellPointsPosition;
  for( localIndex numVertex = 0; numVertex < numCellPoints; ++numVertex )
    cellPointsPosition.insert( std::pair< localIndex, localIndex >
                                 ( cellToNodes( cellIndex, numVertex ), numVertex ));
  // - compute cell diameter.
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

  // Compute basis functions and scaled monomials integrals on the boundary.
  array1d< real64 > basisBoundaryIntegrals( numCellPoints );
  array2d< real64 > basisTimesNormalBoundaryInt( 3, numCellPoints );
  array2d< real64 > basisTimesMonomNormalDerBoundaryInt( numCellPoints, 3 );
  array1d< real64 > monomBoundaryIntegrals( 4 );
  // - initialize vectors
  LvArray::tensorOps::fill< 4 >( monomBoundaryIntegrals, 0.0 );
  for( localIndex numBasisFunction = 0; numBasisFunction < numCellPoints; ++numBasisFunction )
  {
    basisBoundaryIntegrals[numBasisFunction] = 0.0;
    for( localIndex pos = 0; pos < 3; ++pos )
      basisTimesNormalBoundaryInt( pos, numBasisFunction ) = 0.0;
    LvArray::tensorOps::fill< 3 >( basisTimesMonomNormalDerBoundaryInt[numBasisFunction], 0.0 );
  }
  for( localIndex numFace = 0; numFace < numCellFaces; ++numFace )
  {
    localIndex const faceIndex = elementToFaceMap[cellIndex][numFace];
    arraySlice1d< localIndex const > faceToNodes = faceManager.nodeList()[faceIndex];
    // - compute integrals calling auxiliary method
    array1d< real64 > faceBasisIntegrals;
    array1d< real64 > threeDMonomialIntegrals;
    ComputeFaceIntegrals( mesh, faceIndex, invCellDiameter, cellCenter,
                          faceBasisIntegrals, threeDMonomialIntegrals );
    // - get outward face normal
    array1d< real64 > faceNormal( 3 );
    LvArray::tensorOps::copy< 3 >( faceNormal, faceNormals[faceIndex] );
    array1d< real64 > signTestVector( 3 );
    LvArray::tensorOps::copy< 3 >( signTestVector, faceCenters[faceIndex] );
    LvArray::tensorOps::subtract< 3 >( signTestVector, cellCenter );
    if( LvArray::tensorOps::AiBi< 3 >( signTestVector, faceNormal ) < 0 )
      LvArray::tensorOps::scale< 3 >( faceNormal, -1.0 );
    // - add contributions to integrals of monomials
    monomBoundaryIntegrals[0] += faceManager.faceArea()[faceIndex];
    for( localIndex monomInd = 1; monomInd < 4; ++monomInd )
      monomBoundaryIntegrals[monomInd] += threeDMonomialIntegrals[monomInd-1];
    // - add contributions to integrals of basis functions
    for( localIndex numFaceBasisFunction = 0; numFaceBasisFunction < faceToNodes.size();
         ++numFaceBasisFunction )
    {
      localIndex basisFunctionIndex = cellPointsPosition
                                        .find( faceToNodes[numFaceBasisFunction] )->second;
      basisBoundaryIntegrals[basisFunctionIndex] += faceBasisIntegrals[numFaceBasisFunction];
      for( localIndex pos = 0; pos < 3; ++pos )
      {
        basisTimesNormalBoundaryInt( pos, basisFunctionIndex ) +=
          faceNormal( pos )*faceBasisIntegrals( numFaceBasisFunction );
        basisTimesMonomNormalDerBoundaryInt( basisFunctionIndex, pos ) +=
          faceNormal( pos )*faceBasisIntegrals( numFaceBasisFunction )*invCellDiameter;
      }
    }
  }

  // Compute non constant scaled monomials' integrals on the polyhedron.
  array1d< real64 > monomInternalIntegrals( 3 );
  LvArray::tensorOps::fill< 3 >( monomInternalIntegrals, 0.0 );
  m_numQuadraturePoints = 0;
  for( localIndex numFace = 0; numFace < numCellFaces; ++numFace )
  {
    localIndex const faceIndex = elementToFaceMap[cellIndex][numFace];
    arraySlice1d< localIndex const > faceToNodes = faceManager.nodeList()[faceIndex];
    arraySlice1d< real64 const > faceCenter = faceCenters[faceIndex];
    localIndex const numFaceVertices = faceToNodes.size();
    for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
    {
      ++m_numQuadraturePoints;
      localIndex numNextVertex = (numVertex+1)%numFaceVertices;
      // compute value of 3D monomials at the quadrature point on the sub-tetrahedron (the
      // barycenter).
      // The result is ((v0 + v1 + faceCenter + cellCenter)/4 - cellCenter) / cellDiameter =
      // = (v0 + v1 + faceCenter - 3*cellcenter)/(4*cellDiameter).
      array1d< real64 > monomialValues( 3 );
      for( localIndex pos = 0; pos < 3; ++pos )
      {
        monomialValues( pos ) = (nodesCoords( faceToNodes( numVertex ), pos ) +
                                 nodesCoords( faceToNodes( numNextVertex ), pos ) +
                                 faceCenter( pos ) - 3*cellCenter( pos ))*invCellDiameter/4.0;
      }
      // compute quadrature weight (the volume of the sub-tetrahedron).
      array2d< real64 > edgeTangentsMatrix( 3, 3 );
      for( localIndex pos = 0; pos < 3; ++pos )
      {
        edgeTangentsMatrix( 0, pos ) = faceCenter( pos ) - cellCenter( pos );
        edgeTangentsMatrix( 1, pos ) = nodesCoords( faceToNodes( numVertex ), pos ) - cellCenter( pos );
        edgeTangentsMatrix( 2, pos ) = nodesCoords( faceToNodes( numNextVertex ), pos ) -
                                       cellCenter( pos );
      }
      real64 subTetVolume = LvArray::math::abs
                              ( LvArray::tensorOps::determinant< 3 >( edgeTangentsMatrix )) / 6.0;
      for( localIndex pos = 0; pos < 3; ++pos )
        monomInternalIntegrals( pos ) += monomialValues( pos )*subTetVolume;
    }
  }

  // Compute integral mean of basis functions and of derivatives of basis functions.
  // Compute VEM degrees of freedom of the piNabla projection minus the identity (used for
  // m_stabilizationMatrix).
  real64 const invCellVolume = 1.0/cellVolumes[cellIndex];
  real64 const monomialDerivativeInverse = cellDiameter*cellDiameter*invCellVolume;
  m_basisFunctionsIntegralMean.resize( numCellPoints );
  m_basisDerivativesIntegralMean = basisTimesNormalBoundaryInt;
  array2d< real64 > piNablaVemDofsMinusIdentity( numCellPoints, numCellPoints );
  // - compute values of scaled monomials at the vertices (used for piNablaVemDofs)
  array2d< real64 > monomialVemDofs( 3, numCellPoints );
  for( localIndex numVertex = 0; numVertex < numCellPoints; ++numVertex )
  {
    for( localIndex pos = 0; pos < 3; ++pos )
      monomialVemDofs( pos, numVertex ) = invCellDiameter*
                                          (nodesCoords( cellToNodes( cellIndex, numVertex ), pos ) - cellCenter( pos ));
  }
  for( localIndex numBasisFunction = 0; numBasisFunction < numCellPoints; ++numBasisFunction )
  {
    // - solve linear system to obtain dofs of piNabla proj wrt monomial basis
    array1d< real64 > piNablaDofs( 4 );
    piNablaDofs[1] = monomialDerivativeInverse *
                     basisTimesMonomNormalDerBoundaryInt( numBasisFunction, 0 );
    piNablaDofs[2] = monomialDerivativeInverse *
                     basisTimesMonomNormalDerBoundaryInt( numBasisFunction, 1 );
    piNablaDofs[3] = monomialDerivativeInverse *
                     basisTimesMonomNormalDerBoundaryInt( numBasisFunction, 2 );
    piNablaDofs[0] = (basisBoundaryIntegrals[numBasisFunction] -
                      piNablaDofs[1]*monomBoundaryIntegrals[1] -
                      piNablaDofs[2]*monomBoundaryIntegrals[2] -
                      piNablaDofs[3]*monomBoundaryIntegrals[3] )/monomBoundaryIntegrals[0];
    // - integrate piNabla proj and compute integral means
    m_basisFunctionsIntegralMean( numBasisFunction ) = piNablaDofs[0] + invCellVolume *
                                                       (piNablaDofs[1]*monomInternalIntegrals[0] + piNablaDofs[2]*monomInternalIntegrals[1]
                                                        + piNablaDofs[3] * monomInternalIntegrals[2]);
    // - compute integral means of derivatives
    for( localIndex pos = 0; pos < 3; ++pos )
      m_basisDerivativesIntegralMean( pos, numBasisFunction ) =
        -invCellVolume *basisTimesNormalBoundaryInt( pos, numBasisFunction );
    // - compute VEM dofs of piNabla projection
    for( localIndex numVertex = 0; numVertex < numCellPoints; ++numVertex )
    {
      piNablaVemDofsMinusIdentity( numVertex, numBasisFunction ) = piNablaDofs[0] +
                                                                   piNablaDofs[1]*monomialVemDofs( 0, numVertex ) +
                                                                   piNablaDofs[2]*monomialVemDofs( 1, numVertex ) +
                                                                   piNablaDofs[3]*monomialVemDofs( 2, numVertex );
    }
    piNablaVemDofsMinusIdentity( numBasisFunction, numBasisFunction ) -= 1;
  }

  // Compute stabilization matrix.
  // The result is piNablaVemDofsMinusIdentity^T * piNablaVemDofsMinusIdentity.
  m_stabilizationMatrix.resize( numCellPoints, numCellPoints );
  for( localIndex i = 0; i < numCellPoints; ++i )
  {
    for( localIndex j = 0; j < numCellPoints; ++j )
    {
      real64 rowColProd = 0;
      for( localIndex k = 0; k < numCellPoints; ++k )
        rowColProd += piNablaVemDofsMinusIdentity( k, i )*piNablaVemDofsMinusIdentity( k, j );
      m_stabilizationMatrix( i, j ) = cellDiameter*rowColProd;
    }
  }
}

void
ConformingVirtualElementOrder1::
  ComputeFaceIntegrals( MeshLevel const & mesh,
                        localIndex const & faceIndex,
                        real64 const & invCellDiameter,
                        arraySlice1d< real64 const > const & cellCenter,
                        array1d< real64 > & basisIntegrals,
                        array1d< real64 > & threeDMonomialIntegrals )
{
  // Get geometry managers.
  FaceManager const & faceManager = *mesh.getFaceManager();
  NodeManager const & nodeManager = *mesh.getNodeManager();
  EdgeManager const & edgeManager = *mesh.getEdgeManager();

  // Get pre-computed maps.
  arraySlice1d< localIndex const > faceToNodes = faceManager.nodeList()[faceIndex];
  FaceManager::EdgeMapType faceToEdges = faceManager.edgeList();
  EdgeManager::NodeMapType edgeToNodes = edgeManager.nodeList();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > nodesCoords =
    nodeManager.referencePosition();

  // Get pre-computed geometrical properties.
  arrayView2d< real64 const > faceCenters = faceManager.faceCenter();
  arrayView2d< real64 const > faceNormals = faceManager.faceNormal();
  arrayView1d< real64 const > faceAreas = faceManager.faceArea();
  localIndex const numFaceVertices = faceToNodes.size();     // also equal to n. face's edges.

  // Compute other geometrical properties.
  //  - compute rotation matrix.
  array2d< real64 > faceRotationMatrix( 3, 3 );
  computationalGeometry::RotationMatrix_3D( faceNormals[faceIndex], faceRotationMatrix );
  //  - below we compute the diameter, the rotated vertices and the rotated center.
  array2d< real64 > faceRotatedVertices( numFaceVertices, 2 );
  array1d< real64 > faceRotatedCentroid( 2 );
  real64 faceDiameter = 0;
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    // apply the transpose (that is the inverse) of the rotation matrix to face vertices.
    // NOTE:
    // the second and third rows of the transpose of the rotation matrix rotate on the 2D face.
    faceRotatedVertices[numVertex][0] =
      faceRotationMatrix( 0, 1 )*nodesCoords( faceToNodes( numVertex ), 0 ) +
      faceRotationMatrix( 1, 1 )*nodesCoords( faceToNodes( numVertex ), 1 ) +
      faceRotationMatrix( 2, 1 )*nodesCoords( faceToNodes( numVertex ), 2 );
    faceRotatedVertices[numVertex][1] =
      faceRotationMatrix( 0, 2 )*nodesCoords( faceToNodes( numVertex ), 0 ) +
      faceRotationMatrix( 1, 2 )*nodesCoords( faceToNodes( numVertex ), 1 ) +
      faceRotationMatrix( 2, 2 )*nodesCoords( faceToNodes( numVertex ), 2 );
    for( localIndex numOthVertex = 0; numOthVertex < numVertex; ++numOthVertex )
    {
      array1d< real64 > vertDiff( 2 );
      vertDiff( 0 ) = faceRotatedVertices( numVertex, 0 ) - faceRotatedVertices( numOthVertex, 0 );
      vertDiff( 1 ) = faceRotatedVertices( numVertex, 1 ) - faceRotatedVertices( numOthVertex, 1 );
      real64 const candidateDiameter = LvArray::tensorOps::l2NormSquared< 2 >( vertDiff );
      if( faceDiameter < candidateDiameter )
        faceDiameter = candidateDiameter;
    }
  }
  faceDiameter = LvArray::math::sqrt< real64 >( faceDiameter );
  real64 const invFaceDiameter = 1.0/faceDiameter;
  // - rotate the face centroid as done for the vertices.
  faceRotatedCentroid( 0 ) =
    faceRotationMatrix( 0, 1 )*faceCenters( faceIndex, 0 ) +
    faceRotationMatrix( 1, 1 )*faceCenters( faceIndex, 1 ) +
    faceRotationMatrix( 2, 1 )*faceCenters( faceIndex, 2 );
  faceRotatedCentroid( 1 ) =
    faceRotationMatrix( 0, 2 )*faceCenters( faceIndex, 0 ) +
    faceRotationMatrix( 1, 2 )*faceCenters( faceIndex, 1 ) +
    faceRotationMatrix( 2, 2 )*faceCenters( faceIndex, 2 );
  // - compute edges' lengths, outward pointing normals and local edge-to-nodes map.
  array2d< real64 > edgeOutwardNormals( numFaceVertices, 2 );
  array1d< real64 > edgeLengths( numFaceVertices );
  array2d< localIndex > localEdgeToNodes( numFaceVertices, 2 );
  for( localIndex numEdge = 0; numEdge < numFaceVertices; ++numEdge )
  {
    if( edgeToNodes( faceToEdges( faceIndex, numEdge ), 0 ) == faceToNodes( numEdge ))
    {
      localEdgeToNodes( numEdge, 0 ) = numEdge;
      localEdgeToNodes( numEdge, 1 ) = (numEdge+1)%numFaceVertices;
    }
    else
    {
      localEdgeToNodes( numEdge, 0 ) = (numEdge+1)%numFaceVertices;
      localEdgeToNodes( numEdge, 1 ) = numEdge;
    }
    array1d< real64 > edgeTangent( 2 );
    LvArray::tensorOps::copy< 2 >( edgeTangent, faceRotatedVertices[(numEdge+1)%numFaceVertices] );
    LvArray::tensorOps::subtract< 2 >( edgeTangent, faceRotatedVertices[numEdge] );
    edgeOutwardNormals[numEdge][0] = edgeTangent[1];
    edgeOutwardNormals[numEdge][1] = -edgeTangent[0];
    array1d< real64 > signTestVector( 2 );
    LvArray::tensorOps::copy< 2 >( signTestVector, faceRotatedVertices[numEdge] );
    LvArray::tensorOps::subtract< 2 >( signTestVector, faceRotatedCentroid );
    if( LvArray::tensorOps::AiBi< 2 >( signTestVector, edgeOutwardNormals[numEdge] ) < 0 )
      LvArray::tensorOps::scale< 2 >( edgeOutwardNormals[numEdge], -1 );
    LvArray::tensorOps::normalize< 2 >( edgeOutwardNormals[numEdge] );
    edgeLengths[numEdge] = LvArray::tensorOps::l2Norm< 2 >( edgeTangent );
  }

  // Compute boundary quadrature weights (also equal to the integrals of basis functions on the
  // boundary).
  array1d< real64 > boundaryQuadratureWeights( numFaceVertices );
  boundaryQuadratureWeights.setValues< serialPolicy >( 0.0 );
  for( localIndex numEdge = 0; numEdge < numFaceVertices; ++numEdge )
  {
    boundaryQuadratureWeights[localEdgeToNodes( numEdge, 0 )] += 0.5*edgeLengths[numEdge];
    boundaryQuadratureWeights[localEdgeToNodes( numEdge, 1 )] += 0.5*edgeLengths[numEdge];
  }

  // Compute scaled monomials' integrals on edges.
  array1d< real64 > monomBoundaryIntegrals( 3 );
  LvArray::tensorOps::fill< 3 >( monomBoundaryIntegrals, 0.0 );
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    monomBoundaryIntegrals( 0 ) += boundaryQuadratureWeights( numVertex );
    monomBoundaryIntegrals( 1 ) += (faceRotatedVertices( numVertex, 0 ) - faceRotatedCentroid( 0 )) *
                                   invFaceDiameter*boundaryQuadratureWeights( numVertex );
    monomBoundaryIntegrals( 2 ) += (faceRotatedVertices( numVertex, 1 ) - faceRotatedCentroid( 1 )) *
                                   invFaceDiameter*boundaryQuadratureWeights( numVertex );
  }

  // Compute non constant 2D and 3D scaled monomials' integrals on the face.
  array1d< real64 > monomInternalIntegrals( 2 );
  threeDMonomialIntegrals.resize( 3 );
  LvArray::tensorOps::fill< 2 >( monomInternalIntegrals, 0.0 );
  LvArray::tensorOps::fill< 3 >( threeDMonomialIntegrals, 0.0 );
  for( localIndex numSubTriangle = 0; numSubTriangle < numFaceVertices; ++numSubTriangle )
  {
    // - compute value of 2D monomials at the quadrature point on the sub-triangle (the
    //   barycenter).
    //   The result is ((v(0)+v(1)+faceCenter)/3 - faceCenter) / faceDiameter =
    //   = (v(0) + v(1) - 2*faceCenter)/(3*faceDiameter).
    array1d< real64 > monomialValues( 2 );
    LvArray::tensorOps::copy< 2 >( monomialValues, faceRotatedCentroid ); // val = faceCenter
    LvArray::tensorOps::scale< 2 >( monomialValues, -2.0 ); // val = -2*faceCenter
    LvArray::tensorOps::add< 2 >( monomialValues,  // val = v(0) - 2*faceCenter
                                  faceRotatedVertices[numSubTriangle] );
    LvArray::tensorOps::add< 2 >( monomialValues,  // val = v(0) + v(1) - 2*faceCenter
                                  faceRotatedVertices[(numSubTriangle+1)%numFaceVertices] );
    // val = (v(0) + v(1) - 2*faceCenter)/(3*faceDiameter)
    LvArray::tensorOps::scale< 2 >( monomialValues, 1.0/(3.0*faceDiameter));
    // compute value of 3D monomials at the quadrature point on the sub-triangle (the
    // barycenter).  The result is
    // ((v(0) + v(1) + faceCenter)/3 - cellCenter)/cellDiameter.
    array1d< real64 > threeDMonomialValues( 3 );
    LvArray::tensorOps::copy< 3 >( threeDMonomialValues,
                                   faceCenters[faceIndex] ); // val = faceCenter
    LvArray::tensorOps::add< 3 >( threeDMonomialValues,  // val = v(0) + faceCenter
                                  nodesCoords[faceToNodes( numSubTriangle )] );
    LvArray::tensorOps::add< 3 >( threeDMonomialValues,  // val = v(0) + v(1) + faceCenter
                                  nodesCoords[faceToNodes((numSubTriangle+1)%numFaceVertices )] );
    // val = (v(0) + v(1) + faceCenter)/3
    LvArray::tensorOps::scale< 3 >( threeDMonomialValues, 1.0/3.0 );
    // val = (v(0) + v(1) + faceCenter)/3 - cellCenter
    LvArray::tensorOps::subtract< 3 >( threeDMonomialValues, cellCenter );
    // val = ((v(0) + v(1) + faceCenter)/3 - cellCenter)/cellDiameter
    LvArray::tensorOps::scale< 3 >( threeDMonomialValues, invCellDiameter );
    // compute quadrature weight associated to the quadrature point (the area of the
    // sub-triangle).
    array2d< real64 > edgesTangents( 2, 2 ); // used to compute the area of the sub-triangle
    LvArray::tensorOps::copy< 2 >( edgesTangents[0],
                                   faceRotatedVertices[numSubTriangle] );
    LvArray::tensorOps::subtract< 2 >( edgesTangents[0], faceRotatedCentroid );
    LvArray::tensorOps::copy< 2 >( edgesTangents[1],
                                   faceRotatedVertices[(numSubTriangle+1)%numFaceVertices] );
    LvArray::tensorOps::subtract< 2 >( edgesTangents[1], faceRotatedCentroid );
    real64 subTriangleArea = 0.5*LvArray::math::abs
                               ( LvArray::tensorOps::determinant< 2 >( edgesTangents ));
    // compute the integrals on the sub-triangle and add it to the global integrals
    LvArray::tensorOps::scale< 2 >( monomialValues, subTriangleArea );
    LvArray::tensorOps::add< 2 >( monomInternalIntegrals, monomialValues );
    LvArray::tensorOps::scale< 3 >( threeDMonomialValues, subTriangleArea );
    LvArray::tensorOps::add< 3 >( threeDMonomialIntegrals, threeDMonomialValues );
  }

  // Compute integral of basis functions times normal derivative of monomials on the boundary.
  array2d< real64 > basisTimesMonomNormalDerBoundaryInt( numFaceVertices, 2 );
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
    LvArray::tensorOps::fill< 2 >( basisTimesMonomNormalDerBoundaryInt[numVertex], 0.0 );
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    array1d< real64 > thisEdgeIntTimesNormal( 2 );
    LvArray::tensorOps::copy< 2 >( thisEdgeIntTimesNormal, edgeOutwardNormals[numVertex] );
    LvArray::tensorOps::scale< 2 >( thisEdgeIntTimesNormal, edgeLengths[numVertex] );
    LvArray::tensorOps::add< 2 >
      ( basisTimesMonomNormalDerBoundaryInt[localEdgeToNodes( numVertex, 0 )],
      thisEdgeIntTimesNormal );
    LvArray::tensorOps::add< 2 >
      ( basisTimesMonomNormalDerBoundaryInt[localEdgeToNodes( numVertex, 1 )],
      thisEdgeIntTimesNormal );
  }
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
    LvArray::tensorOps::scale< 2 >( basisTimesMonomNormalDerBoundaryInt[numVertex],
                                    0.5*invFaceDiameter );

  // Compute integral mean of basis functions on this face.
  real64 invFaceArea = 1.0/faceAreas[faceIndex];
  real64 monomialDerivativeInverse = (faceDiameter*faceDiameter)*invFaceArea;
  basisIntegrals.resize( numFaceVertices );
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    array1d< real64 > piNablaDofs( 3 );
    piNablaDofs[1] = monomialDerivativeInverse *
                     basisTimesMonomNormalDerBoundaryInt( numVertex, 0 );
    piNablaDofs[2] = monomialDerivativeInverse *
                     basisTimesMonomNormalDerBoundaryInt( numVertex, 1 );
    piNablaDofs[0] = (boundaryQuadratureWeights[numVertex] -
                      piNablaDofs[1]*monomBoundaryIntegrals[1] -
                      piNablaDofs[2]*monomBoundaryIntegrals[2])/monomBoundaryIntegrals[0];
    basisIntegrals[numVertex] = piNablaDofs[0]*faceAreas[faceIndex] +
                                (piNablaDofs[1]*monomInternalIntegrals[0] + piNablaDofs[2]*monomInternalIntegrals[1]);
  }
}
}
}
