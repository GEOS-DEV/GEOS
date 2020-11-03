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
        numQuadraturePoints += faceManager.nodeList()[elementToFaceMap[cellIndex][numFace]].size();
    }
  }
}
