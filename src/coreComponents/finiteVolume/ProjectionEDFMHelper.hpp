#ifndef __PROJECTIONEDFMHELPER_H_
#define __PROJECTIONEDFMHELPER_H_

#include "finiteVolume/FluxApproximationBase.hpp"  // FIXME: include appropriate headers

namespace geosx {

class ProjectionEDFMHelper {
 public:
  ProjectionEDFMHelper(MeshLevel const & m_mesh);
  void addNonNeighboringConnections(EmbeddedSurfaceSubRegion const & fractureSubRegion);

  virtual ~ProjectionEDFMHelper() = default;

 private:
  void selectFaces(FixedOneToManyRelation const & subRegionFaces,
                   localIndex hostCellIdx,
                   localIndex const fracElement,
                   EmbeddedSurfaceSubRegion const & fractureSubRegion) const;

  bool intersection(R1Tensor const & fracOrigin,
                    R1Tensor const & fracNormal,
                    R1Tensor const & vertex1,
                    R1Tensor const & vertex2) const noexcept;

  MeshLevel const & m_mesh;
  ElementRegionManager const * const m_elementManager;
  FaceManager const * const m_faceManager;
  NodeManager const * const m_nodeManager;
  EdgeManager const * const m_edgeManager;

};

}  // end namespace geosx


#endif // __PROJECTIONEDFMHELPER_H_
