#ifndef __PROJECTIONEDFMHELPER_H_
#define __PROJECTIONEDFMHELPER_H_

#include "finiteVolume/FluxApproximationBase.hpp"  // FIXME: include appropriate headers
#include <vector>

namespace geosx {

class ProjectionEDFMHelper {
 public:
  ProjectionEDFMHelper(MeshLevel const & m_mesh);
  void addNonNeighboringConnections(EmbeddedSurfaceSubRegion const & fractureSubRegion);

  virtual ~ProjectionEDFMHelper() = default;

 private:

  struct CellID
  {
    CellID(localIndex r, localIndex sr, localIndex i)
        : region(r), subRegion(sr), index(i)
    {}

    localIndex region;
    localIndex subRegion;
    localIndex index;
  };

  std::vector<localIndex> selectFaces(FixedOneToManyRelation const & subRegionFaces,
                                      CellID const & hostCellID,
                                      localIndex const fracElement,
                                      EmbeddedSurfaceSubRegion const & fractureSubRegion) const;

  bool intersection(R1Tensor const & fracOrigin,
                    R1Tensor const & fracNormal,
                    localIndex edgeIdx,
                    R1Tensor & tmp) const noexcept;

  bool isBoundaryFace(localIndex faceIdx) const noexcept;
  bool onLargerSide( localIndex faceIdx,
                     real64 signedDistanceCellCenterToFrac,
                     R1Tensor const & fracOrigin,
                     R1Tensor const & fracNormal ) const noexcept;
  real64 getSingedDistanceCellCenterToFracPlane( CellID const & hostCellID,
                                                 R1Tensor const & fracNormal,
                                                 R1Tensor const & fracOrigin,
                                                 R1Tensor & tmp) const noexcept;
  bool neighborOnSameSide( localIndex faceIdx,
                           read64 signedDistanceCellCenterToFrac,
                           CellID const & hostCellID ) const;

  MeshLevel const & m_mesh;
  ElementRegionManager const * const m_elementManager;
  FaceManager const * const m_faceManager;
  NodeManager const * const m_nodeManager;
  EdgeManager const * const m_edgeManager;
  // ArrayOfArraysView< localIndex const > const & m_faceToEdges;
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_nodesCoord;
  arrayView2d< localIndex const > const m_edgeToNodes;
  arrayView2d< localIndex const > const m_facesToCells;
  ArrayOfArraysView< localIndex const > const m_facesToNodes;
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > m_nodeReferencePosition;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const m_cellCenters;
};

}  // end namespace geosx


#endif // __PROJECTIONEDFMHELPER_H_
