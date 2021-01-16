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
  std::vector<localIndex> selectFaces(FixedOneToManyRelation const & subRegionFaces,
                                      localIndex hostCellIdx,
                                      localIndex const fracElement,
                                      EmbeddedSurfaceSubRegion const & fractureSubRegion) const;

  bool intersection(R1Tensor const & fracOrigin,
                    R1Tensor const & fracNormal,
                    localIndex edgeIdx,
                    R1Tensor & tmp) const noexcept;

  bool isBoundaryFace(localIndex faceIdx) const noexcept;

  MeshLevel const & m_mesh;
  ElementRegionManager const * const m_elementManager;
  FaceManager const * const m_faceManager;
  NodeManager const * const m_nodeManager;
  EdgeManager const * const m_edgeManager;
  // ArrayOfArraysView< localIndex const > const & m_faceToEdges;
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_nodesCoord;
  arrayView2d< localIndex const > const m_edgeToNodes;
  arrayView2d< localIndex const > const m_facesToCells;

};

}  // end namespace geosx


#endif // __PROJECTIONEDFMHELPER_H_
