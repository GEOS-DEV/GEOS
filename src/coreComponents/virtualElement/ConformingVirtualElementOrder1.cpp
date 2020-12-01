#include "ConformingVirtualElementOrder1.hpp"

namespace geosx
{
  namespace virtualElement
  {
    void ConformingVirtualElementOrder1::ComputeProjectors( MeshLevel const & mesh,
                                                            localIndex const & regionIndex,
                                                            localIndex const & subRegionIndex,
                                                            localIndex const & cellIndex)
    {
      // Get geometry managers
      NodeManager const & nodeManager = *mesh.getNodeManager();
      FaceManager const & faceManager = *mesh.getFaceManager();
      ElementRegionManager const & elementManager = *mesh.getElemManager();

      // Get pre-computed maps
      CellElementRegion const & cellRegion =
        *elementManager.GetRegion<CellElementRegion>(regionIndex);
      CellElementSubRegion const & cellSubRegion =
        *cellRegion.GetSubRegion<CellElementSubRegion>(subRegionIndex);
      arraySlice1d< localIndex const> cellToNodes = cellSubRegion.nodeList()[cellIndex];
      FixedOneToManyRelation const & elementToFaceMap = cellSubRegion.faceList();
      arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > nodesCoords =
        nodeManager.referencePosition();

      // Get pre-computed geometrical properties.
      arrayView2d< real64 const > faceCenters = faceManager.faceCenter();
      arrayView2d< real64 const > faceNormals = faceManager.faceNormal();
      arrayView2d< real64 const > cellCenters = cellSubRegion.getElementCenter();
      arrayView1d< real64 const > cellVolumes =	cellSubRegion.getElementVolume();
      arraySlice1d< real64 const > cellCenter = cellCenters[cellIndex];
      localIndex const numCellFaces = elementToFaceMap[cellIndex].size();
      localIndex const numCellPoints = cellToNodes.size();
      numSupportPoints = numCellPoints;

      // Compute other geometrical properties.
      //  - compute map used to locate local point position by global id (used in the computation of
      //    basis functions boundary integrals.
      map<localIndex, localIndex> cellPointsPosition;
      for(unsigned int numVertex = 0; numVertex < numCellPoints; ++numVertex)
        cellPointsPosition.insert(std::pair<localIndex, localIndex>
                                  (cellToNodes[numVertex], numVertex));
      // - compute cell diameter.
      real64 cellDiameter = 0;
      for(unsigned int numVertex = 0; numVertex < numCellPoints; ++numVertex)
      {
        for(localIndex numOthVertex = 0; numOthVertex < numVertex; ++numOthVertex)
        {
          array1d<real64> vertDiff(3);
          LvArray::tensorOps::copy<3>(vertDiff, nodesCoords[cellToNodes(numVertex)]);
          LvArray::tensorOps::subtract<3>(vertDiff, nodesCoords[cellToNodes(numOthVertex)]);
          real64 const candidateDiameter = LvArray::tensorOps::l2NormSquared<3>(vertDiff);
          if(cellDiameter < candidateDiameter)
            cellDiameter = candidateDiameter;
        }
      }
      cellDiameter = LvArray::math::sqrt<real64>(cellDiameter);
      real64 const invCellDiameter = 1.0/cellDiameter;

      // Compute basis functions and scaled monomials integrals on the boundary.
      array1d<real64> basisBoundaryIntegrals(numCellPoints);
      array2d<real64> basisNormalBoundaryInt(numCellPoints, 3);
      array1d<real64> monomBoundaryIntegrals(4);
      // - initialize vectors
      LvArray::tensorOps::fill<4>(monomBoundaryIntegrals, 0.0);
      for(localIndex numBasisFunction = 0; numBasisFunction < numCellPoints; ++numBasisFunction)
      {
        basisBoundaryIntegrals[numBasisFunction] = 0;
        LvArray::tensorOps::fill<3>(basisNormalBoundaryInt[numBasisFunction], 0.0);
      }
      for(localIndex numFace = 0; numFace < numCellFaces; ++numFace)
      {
        localIndex const faceIndex = elementToFaceMap[cellIndex][numFace];
        arraySlice1d< localIndex const > faceToNodes = faceManager.nodeList()[faceIndex];
        // - compute integrals calling auxiliary method
        array1d<real64> faceBasisIntegrals;
        real64 faceDiameter;
        array1d<real64> threeDMonomialIntegrals;
        ComputeFaceIntegrals(mesh, faceIndex, invCellDiameter, cellCenter,
                             faceDiameter, faceBasisIntegrals, threeDMonomialIntegrals);
        // - get outward face normal
        array1d<real64> faceNormal(3);
        LvArray::tensorOps::copy<3>(faceNormal, faceNormals[faceIndex]);
        array1d<real64> signTestVector(3);
        LvArray::tensorOps::copy<3>(signTestVector, faceCenters[faceIndex]);
        LvArray::tensorOps::subtract<3>(signTestVector, cellCenter);
        if(LvArray::tensorOps::AiBi<3>(signTestVector, faceNormal) < 0)
          LvArray::tensorOps::scale<3>(faceNormal, -1.0);
        // - add contributions to integrals of monomials
        monomBoundaryIntegrals[0] += faceManager.faceArea()[faceIndex];
        for(localIndex monomInd = 1; monomInd < 4; ++monomInd)
          monomBoundaryIntegrals[monomInd] += threeDMonomialIntegrals[monomInd-1];
        // - add contributions to integrals of basis functions
        real64 const invFaceDiameter = 1.0/faceDiameter; // derivative of monomials
        for(localIndex numFaceBasisFunction = 0; numFaceBasisFunction < faceToNodes.size();
            ++numFaceBasisFunction)
        {
          localIndex basisFunctionIndex = cellPointsPosition
            .find(faceToNodes[numFaceBasisFunction])->second;
          basisBoundaryIntegrals[basisFunctionIndex] += faceBasisIntegrals[numFaceBasisFunction];
          array1d<real64> normalIntegralContribution = faceNormal;
          LvArray::tensorOps::scale<3>(normalIntegralContribution,
                                       faceBasisIntegrals[numFaceBasisFunction]*invFaceDiameter);
          LvArray::tensorOps::add<3>(basisNormalBoundaryInt[basisFunctionIndex],
                                     normalIntegralContribution);
        }
      }

      // Compute non constant scaled monomials' integrals on the polyhedron.
      array1d<real64> monomInternalIntegrals(3);
      LvArray::tensorOps::fill<3>(monomInternalIntegrals, 0.0);
      for(localIndex numFace = 0; numFace < numCellFaces; ++numFace)
      {
        localIndex const faceIndex = elementToFaceMap[cellIndex][numFace];
        arraySlice1d< localIndex const > faceToNodes = faceManager.nodeList()[faceIndex];
        arraySlice1d< real64 const > faceCenter = faceCenters[faceIndex];
        localIndex const numFaceVertices = faceToNodes.size();
        for(localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex)
        {
          localIndex numNextVertex = (numVertex+1)%numFaceVertices;
          // compute value of 3D monomials at the quadrature point on the sub-tetrahedron (the
          // barycenter).
          // The result is ((v0 + v1 + faceCenter + cellCenter)/4 - cellCenter) / cellDiameter =
          // = (v0 + v1 + faceCenter - 3*cellcenter)/(4*cellDiameter).
          array1d<real64> monomialValues(3);
          for(localIndex pos = 0; pos < 3; ++pos)
          {
            monomialValues(pos) = (nodesCoords(faceToNodes(numVertex), pos) +
                                   nodesCoords(faceToNodes(numNextVertex), pos) +
                                   faceCenter(pos) - 3*cellCenter(pos))*invCellDiameter/4.0;
          }
          // compute quadrature weight (the volume of the sub-tetrahedron).
          array2d<real64> edgeTangentsMatrix(3, 3);
          for(localIndex pos = 0; pos < 3; ++pos)
          {
            edgeTangentsMatrix(0, pos) = faceCenter(pos) - cellCenter(pos);
            edgeTangentsMatrix(1, pos) = nodesCoords(faceToNodes(numVertex), pos) - cellCenter(pos);
            edgeTangentsMatrix(2, pos) = nodesCoords(faceToNodes(numNextVertex), pos) -
              cellCenter(pos);
          }
          real64 subTetVolume = LvArray::math::abs
            (LvArray::tensorOps::determinant<3>(edgeTangentsMatrix)) / 6.0;
          for(localIndex pos = 0; pos < 3; ++pos)
            monomInternalIntegrals(pos) += monomialValues(pos)*subTetVolume;
        }
      }

      // Compute integral mean of basis functions.
      real64 const monomialDerivativeInverse = cellDiameter*cellDiameter/cellVolumes[cellIndex];
      basisFunctionsIntegralMean.resize(numCellPoints);
      for(localIndex numVertex = 0; numVertex < numCellPoints; ++numVertex)
      {
        array1d<real64> piNablaDofs(4);
        piNablaDofs[1] = monomialDerivativeInverse * basisNormalBoundaryInt(numVertex, 0);
        piNablaDofs[2] = monomialDerivativeInverse * basisNormalBoundaryInt(numVertex, 1);
        piNablaDofs[3] = monomialDerivativeInverse * basisNormalBoundaryInt(numVertex, 2);
        piNablaDofs[0] = (basisBoundaryIntegrals[numVertex] -
                          piNablaDofs[1]*monomBoundaryIntegrals[1] -
                          piNablaDofs[2]*monomBoundaryIntegrals[2] -
                          piNablaDofs[3]*monomBoundaryIntegrals[3] )/monomBoundaryIntegrals[0];
        basisFunctionsIntegralMean(numVertex) = piNablaDofs[0] + (1/cellVolumes[cellIndex]) *
          (piNablaDofs[1]*monomInternalIntegrals[0] + piNablaDofs[2]*monomInternalIntegrals[1]
           + piNablaDofs[3] * monomInternalIntegrals[2]);
      }
    }

    void
    ConformingVirtualElementOrder1::ComputeFaceIntegrals( MeshLevel const & mesh,
                                                          localIndex const & faceId,
                                                          real64 const & invCellDiameter,
                                                          arraySlice1d<real64 const> const & cellCenter,
                                                          real64 & faceDiameter,
                                                          array1d<real64> & basisIntegrals,
                                                          array1d<real64> & threeDMonomialIntegrals )
    {
      // Get geometry managers.
      FaceManager const & faceManager = *mesh.getFaceManager();
      NodeManager const & nodeManager = *mesh.getNodeManager();
      EdgeManager const & edgeManager = *mesh.getEdgeManager();

      // Get pre-computed maps.
      arraySlice1d< localIndex const > faceToNodes = faceManager.nodeList()[faceId];
      FaceManager::EdgeMapType faceToEdges = faceManager.edgeList();
      EdgeManager::NodeMapType edgeToNodes = edgeManager.nodeList();
      arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > nodesCoords =
        nodeManager.referencePosition();

      // Get pre-computed geometrical properties.
      arrayView3d< real64 const > faceRotationMatrices = faceManager.faceRotationMatrix();
      arrayView2d< real64 const > faceCenters = faceManager.faceCenter();
      arrayView2d< real64 const > faceNormals = faceManager.faceNormal();
      arrayView1d< real64 const > faceAreas = faceManager.faceArea();
      localIndex const numFaceVertices = faceToNodes.size(); // also equal to n. face's edges.

      // Compute other geometrical properties.
      //  - below we compute the diameter, the rotated vertices and the rotated center.
      array2d<real64> faceRotatedVertices(numFaceVertices, 2);
      array1d<real64> faceRotatedCentroid(2);
      faceDiameter = 0;
      for(unsigned int numVertex = 0; numVertex < numFaceVertices; ++numVertex)
      {
        // apply the transpose (that is the inverse) of the rotation matrix to face vertices.
        // NOTE:
        // the second and third rows of the transpose of the rotation matrix rotate on the 2D face.
        faceRotatedVertices[numVertex][0] =
          faceRotationMatrices(faceId, 0, 1)*nodesCoords(faceToNodes(numVertex), 0) +
          faceRotationMatrices(faceId, 1, 1)*nodesCoords(faceToNodes(numVertex), 1) +
          faceRotationMatrices(faceId, 2, 1)*nodesCoords(faceToNodes(numVertex), 2);
        faceRotatedVertices[numVertex][1] =
          faceRotationMatrices(faceId, 0, 2)*nodesCoords(faceToNodes(numVertex), 0) +
          faceRotationMatrices(faceId, 1, 2)*nodesCoords(faceToNodes(numVertex), 1) +
          faceRotationMatrices(faceId, 2, 2)*nodesCoords(faceToNodes(numVertex), 2);
        for(localIndex numOthVertex = 0; numOthVertex < numVertex; ++numOthVertex)
        {
          array1d<real64> vertDiff(2);
          vertDiff(0) = faceRotatedVertices(numVertex, 0) - faceRotatedVertices(numOthVertex, 0);
          vertDiff(1) = faceRotatedVertices(numVertex, 1) - faceRotatedVertices(numOthVertex, 1);
          real64 const candidateDiameter = LvArray::tensorOps::l2NormSquared< 2 >(vertDiff);
          if(faceDiameter < candidateDiameter)
            faceDiameter = candidateDiameter;
        }
      }
      faceDiameter = LvArray::math::sqrt<real64>(faceDiameter);
      real64 const invFaceDiameter = 1.0/faceDiameter;
      // - rotate the face centroid as done for the vertices.
      faceRotatedCentroid(0) =
        faceRotationMatrices(faceId, 0, 1)*faceCenters(faceId, 0) +
        faceRotationMatrices(faceId, 1, 1)*faceCenters(faceId, 1) +
        faceRotationMatrices(faceId, 2, 1)*faceCenters(faceId, 2);
      faceRotatedCentroid(1) =
        faceRotationMatrices(faceId, 0, 2)*faceCenters(faceId, 0) +
        faceRotationMatrices(faceId, 1, 2)*faceCenters(faceId, 1) +
        faceRotationMatrices(faceId, 2, 2)*faceCenters(faceId, 2);
      // - compute edges' lengths, outward pointing normals and local edge-to-nodes map.
      array2d<real64> edgeOutwardNormals(numFaceVertices, 2);
      array1d<real64> edgeLengths(numFaceVertices);
      array2d<localIndex> localEdgeToNodes(numFaceVertices, 2);
      for(unsigned int numEdge = 0; numEdge < numFaceVertices; ++numEdge)
      {
        if(edgeToNodes(faceToEdges(faceId, numEdge), 0) == faceToNodes(numEdge))
        {
          localEdgeToNodes(numEdge, 0) = numEdge;
          localEdgeToNodes(numEdge, 1) = (numEdge+1)%numFaceVertices;
        }
        else
        {
          localEdgeToNodes(numEdge, 0) = (numEdge+1)%numFaceVertices;
          localEdgeToNodes(numEdge, 1) = numEdge;
        }
        array1d<real64> edgeTangent(2);
        LvArray::tensorOps::copy<2>(edgeTangent, faceRotatedVertices[(numEdge+1)%numFaceVertices]);
        LvArray::tensorOps::subtract<2>(edgeTangent, faceRotatedVertices[numEdge]);
        edgeOutwardNormals[numEdge][0] = edgeTangent[1];
        edgeOutwardNormals[numEdge][1] = -edgeTangent[0];
        array1d<real64> signTestVector(2);
        LvArray::tensorOps::copy<2>(signTestVector, faceRotatedVertices[numEdge]);
        LvArray::tensorOps::subtract<2>(signTestVector, faceRotatedCentroid);
        if(LvArray::tensorOps::AiBi<2>(signTestVector, edgeOutwardNormals[numEdge]) < 0)
          LvArray::tensorOps::scale<2>(edgeOutwardNormals[numEdge], -1);
        LvArray::tensorOps::normalize<2>(edgeOutwardNormals[numEdge]);
        edgeLengths[numEdge] = LvArray::tensorOps::l2Norm<2>(edgeTangent);
      }

      // Compute boundary quadrature weights (also equal to the integrals of basis functions on the
      // boundary).
      array1d<real64> boundaryQuadratureWeights(numFaceVertices);
      boundaryQuadratureWeights.setValues<serialPolicy>(0.0);
      for(localIndex numEdge = 0; numEdge < numFaceVertices; ++numEdge)
      {
        boundaryQuadratureWeights[localEdgeToNodes(numEdge, 0)] += 0.5*edgeLengths[numEdge];
        boundaryQuadratureWeights[localEdgeToNodes(numEdge, 1)] += 0.5*edgeLengths[numEdge];
      }

      // Compute scaled monomials' integrals on edges.
      array1d<real64> monomBoundaryIntegrals(3);
      LvArray::tensorOps::fill<3>(monomBoundaryIntegrals, 0.0);
      for(localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex)
      {
        monomBoundaryIntegrals(0) += boundaryQuadratureWeights(numVertex);
        monomBoundaryIntegrals(1) += (faceRotatedVertices(numVertex, 0) - faceRotatedCentroid(0)) *
          invFaceDiameter*boundaryQuadratureWeights(numVertex);
        monomBoundaryIntegrals(2) += (faceRotatedVertices(numVertex, 1) - faceRotatedCentroid(1)) *
          invFaceDiameter*boundaryQuadratureWeights(numVertex);
      }

      // Compute non constant 2D and 3D scaled monomials' integrals on the face.
      array1d<real64> monomInternalIntegrals(2);
      threeDMonomialIntegrals.resize(3);
      LvArray::tensorOps::fill<2>(monomInternalIntegrals, 0.0);
      LvArray::tensorOps::fill<3>(threeDMonomialIntegrals, 0.0);
      for(localIndex numSubTriangle = 0; numSubTriangle < numFaceVertices; ++numSubTriangle)
      {
        // compute value of 2D monomials at the quadrature point on the sub-triangle (the barycenter).
        // The result is ((v(0)+v(1)+faceCenter)/3 - faceCenter) / faceDiameter =
        //  = (v(0) + v(1) - 2*faceCenter)/(3*faceDiameter).
        array1d<real64> monomialValues(2);
        LvArray::tensorOps::copy<2>(monomialValues, faceRotatedCentroid); // val = faceCenter
        LvArray::tensorOps::scale<2>(monomialValues, -2.0); // val = -2*faceCenter
        LvArray::tensorOps::add<2>(monomialValues, // val = v(0) - 2*faceCenter
                                   faceRotatedVertices[numSubTriangle]);
        LvArray::tensorOps::add<2>(monomialValues, // val = v(0) + v(1) - 2*faceCenter
                                   faceRotatedVertices[(numSubTriangle+1)%numFaceVertices]);
        // val = (v(0) + v(1) - 2*faceCenter)/(3*faceDiameter)
        LvArray::tensorOps::scale<2>(monomialValues, 1.0/(3.0*faceDiameter));
        // compute value of 3D monomials at the quadrature point on the sub-triangle (the
        // barycenter).  The result is
        // ((v(0) + v(1) + faceCenter)/3 - cellCenter)/cellDiameter.
        array1d<real64> threeDMonomialValues(3);
        LvArray::tensorOps::copy<3>(threeDMonomialValues, faceCenters[faceId]); // val = faceCenter
        LvArray::tensorOps::add<3>(threeDMonomialValues, // val = v(0) + faceCenter
                                   nodesCoords[faceToNodes(numSubTriangle)]);
        LvArray::tensorOps::add<3>(threeDMonomialValues, // val = v(0) + v(1) + faceCenter
                                   nodesCoords[faceToNodes((numSubTriangle+1)%numFaceVertices)]);
        // val = (v(0) + v(1) + faceCenter)/3
        LvArray::tensorOps::scale<3>(threeDMonomialValues, 1.0/3.0);
        // val = (v(0) + v(1) + faceCenter)/3 - cellCenter
        LvArray::tensorOps::subtract<3>(threeDMonomialValues, cellCenter);
        // val = ((v(0) + v(1) + faceCenter)/3 - cellCenter)/cellDiameter
        LvArray::tensorOps::scale<3>(threeDMonomialValues, invCellDiameter);
        // compute quadrature weight associated to the quadrature point (the area of the
        // sub-triangle).
        array2d<real64> edgesTangents(2,2); // used to compute the area of the sub-triangle
        LvArray::tensorOps::copy<2>(edgesTangents[0],
                                    faceRotatedVertices[numSubTriangle]);
        LvArray::tensorOps::subtract<2>(edgesTangents[0], faceRotatedCentroid);
        LvArray::tensorOps::copy<2>(edgesTangents[1],
                                    faceRotatedVertices[(numSubTriangle+1)%numFaceVertices]);
        LvArray::tensorOps::subtract<2>(edgesTangents[1], faceRotatedCentroid);
        real64 subTriangleArea = 0.5*LvArray::math::abs
          (LvArray::tensorOps::determinant<2>(edgesTangents));
        // compute the integrals on the sub-triangle and add it to the global integrals
        LvArray::tensorOps::scale<2>(monomialValues, subTriangleArea);
        LvArray::tensorOps::add<2>(monomInternalIntegrals, monomialValues);
        LvArray::tensorOps::scale<3>(threeDMonomialValues, subTriangleArea);
        LvArray::tensorOps::add<3>(threeDMonomialIntegrals, threeDMonomialValues);
      }

      // Compute integral of basis functions times normal derivative of monomials on the boundary.
      array2d<real64> basisNormalBoundaryInt(numFaceVertices, 2);
      for(localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex)
        LvArray::tensorOps::fill<2>(basisNormalBoundaryInt[numVertex], 0.0);
      for(localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex)
      {
        array1d<real64> thisEdgeIntTimesNormal(2);
        LvArray::tensorOps::copy<2>(thisEdgeIntTimesNormal, edgeOutwardNormals[numVertex]);
        LvArray::tensorOps::scale<2>(thisEdgeIntTimesNormal, edgeLengths[numVertex]);
        LvArray::tensorOps::add<2>
          (basisNormalBoundaryInt[localEdgeToNodes(numVertex, 0)], thisEdgeIntTimesNormal);
        LvArray::tensorOps::add<2>
          (basisNormalBoundaryInt[localEdgeToNodes(numVertex, 1)], thisEdgeIntTimesNormal);
      }
      for(localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex)
        LvArray::tensorOps::scale<2>(basisNormalBoundaryInt[numVertex], 0.5*invFaceDiameter);

      // Compute integral mean of basis functions on this face.
      real64 invFaceArea = 1.0/faceAreas[faceId];
      real64 monomialDerivativeInverse = (faceDiameter*faceDiameter)*invFaceArea;
      basisIntegrals.resize(numFaceVertices);
      for(localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex)
      {
        array1d<real64> piNablaDofs(3);
        piNablaDofs[1] = monomialDerivativeInverse * basisNormalBoundaryInt(numVertex, 0);
        piNablaDofs[2] = monomialDerivativeInverse * basisNormalBoundaryInt(numVertex, 1);
        piNablaDofs[0] = (boundaryQuadratureWeights[numVertex] -
                          piNablaDofs[1]*monomBoundaryIntegrals[1] -
                          piNablaDofs[2]*monomBoundaryIntegrals[2])/monomBoundaryIntegrals[0];
        basisIntegrals[numVertex] = piNablaDofs[0]*faceAreas[faceId] +
          (piNablaDofs[1]*monomInternalIntegrals[0] + piNablaDofs[2]*monomInternalIntegrals[1]);
      }
    }
  }
}
