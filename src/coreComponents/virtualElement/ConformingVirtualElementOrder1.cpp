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
      FaceManager const & faceManager = *mesh.getFaceManager();
      ElementRegionManager const & elementManager = *mesh.getElemManager();
      CellElementRegion const & cellRegion =
        *elementManager.GetRegion<CellElementRegion>(regionIndex);
      CellElementSubRegion const & cellSubRegion =
        *cellRegion.GetSubRegion<CellElementSubRegion>(subRegionIndex);
      CellElementSubRegion::NodeMapType const & elementToNodeMap = cellSubRegion.nodeList();
      FixedOneToManyRelation const & elementToFaceMap = cellSubRegion.faceList();
      localIndex const numCellFaces = elementToFaceMap[cellIndex].size();
      localIndex const numCellPoints = elementToNodeMap[cellIndex].size();
      numSupportPoints = elementToNodeMap[cellIndex].size();
      numQuadraturePoints = 0;
      map<localIndex, localIndex> cellPointsPosition;
      for(unsigned int numVertex = 0; numVertex < numCellPoints; ++numVertex)
        cellPointsPosition.insert(std::pair<localIndex, localIndex>(elementToNodeMap[cellIndex][numVertex], numVertex));
      array1d<real64> basisBoundaryIntegrals(numCellPoints);
      for(localIndex numBasisFunction = 0; numBasisFunction < numCellPoints; ++numBasisFunction)
        basisBoundaryIntegrals[numBasisFunction] = 0;
      for(localIndex numFace = 0; numFace < numCellFaces; ++numFace)
      {
        localIndex const faceIndex = elementToFaceMap[cellIndex][numFace];
        numQuadraturePoints += faceManager.nodeList()[elementToFaceMap[cellIndex][numFace]].size();
        array1d<real64> faceBasisProjections;
        ComputeFaceProjectors(mesh, faceIndex, faceBasisProjections);
        arraySlice1d< localIndex const > faceToNodes = faceManager.nodeList()[faceIndex];
        for(localIndex numFaceBasisFunctions = 0; numFaceBasisFunctions < faceToNodes.size(); ++numFaceBasisFunctions)
        {
          localIndex basisFunctionIndex = cellPointsPosition.find(faceToNodes[numFaceBasisFunctions])->second;
          basisBoundaryIntegrals[basisFunctionIndex] += faceManager.faceArea()[faceIndex] * faceBasisProjections[numFaceBasisFunctions];
        }
      }
    }

    void
    ConformingVirtualElementOrder1::ComputeFaceProjectors( MeshLevel const & mesh,
                                                           localIndex const & faceId,
                                                           array1d<real64> & basisProjections)
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
      arrayView1d< real64 const > faceAreas = faceManager.faceArea();
      localIndex const numFaceVertices = faceToNodes.size(); // also equal to n. face's edges.

      // Compute other geometrical properties.
      //  - below we compute the diameter, the rotated vertices and the rotated center.
      array2d<real64> faceRotatedVertices(numFaceVertices, 2);
      array1d<real64> faceRotatedCentroid(2);
      real64 faceDiameter = 0;
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
      for(unsigned int numVertex = 0; numVertex < numFaceVertices; ++numVertex)
      {
        boundaryQuadratureWeights[numVertex] =
          0.5*(edgeLengths[numVertex] + edgeLengths[(numVertex+1)%numFaceVertices]);
      }

      // Compute scaled monomials' boundary integrals.
      array1d<real64> monomBoundaryIntegrals(3);
      LvArray::tensorOps::fill<3>(monomBoundaryIntegrals, 0.0);
      for(unsigned int numVertex = 0; numVertex < numFaceVertices; ++numVertex)
      {
        monomBoundaryIntegrals(0) += boundaryQuadratureWeights(numVertex);
        monomBoundaryIntegrals(1) += (faceRotatedVertices(numVertex, 0) - faceRotatedCentroid(0)) *
          invFaceDiameter*boundaryQuadratureWeights(numVertex);
        monomBoundaryIntegrals(2) += (faceRotatedVertices(numVertex, 1) - faceRotatedCentroid(1)) *
          invFaceDiameter*boundaryQuadratureWeights(numVertex);
      }

      // Compute non constant scaled monomials' integrals on the face.
      array1d<real64> monomInternalIntegrals(2);
      LvArray::tensorOps::fill<2>(monomInternalIntegrals, 0.0);
      for(unsigned int numSubTriangle = 0; numSubTriangle < numFaceVertices; ++numSubTriangle)
      {
        // compute value of monomials at the quadrature point on the sub-triangle (the barycenter)
        // the result is (v(0) + v(1) + center)/3 - center = (v(0) + v(1) - 2*center)/3
        array1d<real64> monomialValues(2);
        LvArray::tensorOps::copy<2>(monomialValues, faceRotatedCentroid); // val = center
        LvArray::tensorOps::scale<2>(monomialValues, -2.0); // val = -2*center
        LvArray::tensorOps::add<2>(monomialValues, // val = v(0) - 2*center
                                   faceRotatedVertices[localEdgeToNodes(numSubTriangle, 0)]);
        LvArray::tensorOps::add<2>(monomialValues, // val = v(0) + v(1) - 2*center
                                   faceRotatedVertices[localEdgeToNodes(numSubTriangle, 1)]);
        LvArray::tensorOps::scale<2>(monomialValues, 1.0/3.0); // val = (v(0) + v(1) - 2*center)/3
        // compute quadrature weight associated to the quadrature point (the area of the
        // sub-triangle).
        array2d<real64> edgesTangents(2,2); // used to compute the area of the sub-triangle
        LvArray::tensorOps::copy<2>(edgesTangents[0],
                                    faceRotatedVertices[localEdgeToNodes(numSubTriangle, 0)]);
        LvArray::tensorOps::subtract<2>(edgesTangents[0], faceRotatedCentroid);
        LvArray::tensorOps::copy<2>(edgesTangents[1],
                                    faceRotatedVertices[localEdgeToNodes(numSubTriangle, 1)]);
        LvArray::tensorOps::subtract<2>(edgesTangents[1], faceRotatedCentroid);
        real64 subTriangleArea = 0.5*LvArray::math::abs
          (LvArray::tensorOps::determinant<2>(edgesTangents));
        // compute the integrals on the sub-triangle and add it to the global integral
        LvArray::tensorOps::scale<2>(monomialValues, subTriangleArea);
        LvArray::tensorOps::add<2>(monomInternalIntegrals, monomialValues);
      }

      // Compute integral of basis functions times normal derivative of monomials on the boundary.
      array2d<real64> basisNormalBoundaryInt(numFaceVertices, 2);
      for(unsigned int numVertex = 0; numVertex < numFaceVertices; ++numVertex)
        LvArray::tensorOps::fill<2>(basisNormalBoundaryInt[numVertex], 0.0);
      for(unsigned int numVertex = 0; numVertex < numFaceVertices; ++numVertex)
      {
        array1d<real64> thisEdgeIntTimesNormal(2);
        LvArray::tensorOps::copy<2>(thisEdgeIntTimesNormal, edgeOutwardNormals[numVertex]);
        LvArray::tensorOps::scale<2>(thisEdgeIntTimesNormal, edgeLengths[numVertex]);
        LvArray::tensorOps::add<2>
          (basisNormalBoundaryInt[localEdgeToNodes(numVertex, 0)], thisEdgeIntTimesNormal);
        LvArray::tensorOps::add<2>
          (basisNormalBoundaryInt[localEdgeToNodes(numVertex, 1)], thisEdgeIntTimesNormal);
      }
      for(unsigned int numVertex = 0; numVertex < numFaceVertices; ++numVertex)
        LvArray::tensorOps::scale<2>(basisNormalBoundaryInt[numVertex], 0.5*invFaceDiameter);

      // Compute integral mean of basis functions on this face.
      real64 invFaceArea = 1.0/faceAreas[faceId];
      real64 monomialDerivativeInverse = (faceDiameter*faceDiameter)*invFaceArea;
      basisProjections.resize(numFaceVertices);
      for(unsigned int numVertex = 0; numVertex < numFaceVertices; ++numVertex)
      {
        array1d<real64> piNablaDofs(3);
        piNablaDofs[1] = monomialDerivativeInverse * basisNormalBoundaryInt(numVertex, 0);
        piNablaDofs[2] = monomialDerivativeInverse * basisNormalBoundaryInt(numVertex, 1);
        piNablaDofs[0] = (boundaryQuadratureWeights[numVertex] -
                          piNablaDofs[1]*monomBoundaryIntegrals[1] -
                          piNablaDofs[2]*monomBoundaryIntegrals[2])/monomBoundaryIntegrals[0];
        basisProjections[numVertex] = piNablaDofs[0] + invFaceArea *
          (piNablaDofs[1]*monomInternalIntegrals[0] + piNablaDofs[2]*monomInternalIntegrals[1]);
      }
    }
  }
}
