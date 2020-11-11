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
      numSupportPoints = elementToNodeMap[cellIndex].size();
      numQuadraturePoints = 0;
      for(localIndex numFace = 0; numFace < numCellFaces; ++numFace)
      {
        numQuadraturePoints += faceManager.nodeList()[elementToFaceMap[cellIndex][numFace]].size();
        ComputeFaceProjectors(mesh, elementToFaceMap[cellIndex][numFace]);
      }
    }

    void ConformingVirtualElementOrder1::ComputeFaceProjectors( MeshLevel const & mesh, localIndex const & faceId )
    {
      // Get geometry managers.
      FaceManager const & faceManager = *mesh.getFaceManager();
      NodeManager const & nodeManager = *mesh.getNodeManager();
      EdgeManager const & edgeManager = *mesh.getEdgeManager();

      // Get pre-compute maps.
      arraySlice1d< localIndex const > faceToNodes = faceManager.nodeList()[faceId];
      FaceManager::EdgeMapType faceToEdges = faceManager.edgeList();
      EdgeManager::NodeMapType edgeToNodes = edgeManager.nodeList();
      arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > nodesCoords = nodeManager.referencePosition();

      // Get pre-computed geometrical properties.
      arrayView3d< real64 const > faceRotationMatrices = faceManager.faceRotationMatrix();
      arrayView2d< real64 const > faceCenters = faceManager.faceCenter();
      arrayView1d< real64 const > faceAreas = faceManager.faceArea();
      localIndex const numFaceVertices = faceToNodes.size(); // also equal to n. face's edges.

      // Compute other geometrical properties.
      // Below we compute the diameter, the rotated vertices and the rotated center.
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
      // rotate the face centroid as done for the vertices.
      faceRotatedCentroid(0) =
        faceRotationMatrices(faceId, 0, 1)*faceCenters(faceId, 0) +
        faceRotationMatrices(faceId, 1, 1)*faceCenters(faceId, 1) +
        faceRotationMatrices(faceId, 2, 1)*faceCenters(faceId, 2);
      faceRotatedCentroid(1) =
        faceRotationMatrices(faceId, 0, 2)*faceCenters(faceId, 0) +
        faceRotationMatrices(faceId, 1, 2)*faceCenters(faceId, 1) +
        faceRotationMatrices(faceId, 2, 2)*faceCenters(faceId, 2);



      //
      GEOSX_UNUSED_VAR(faceAreas);
    }
  }
}
