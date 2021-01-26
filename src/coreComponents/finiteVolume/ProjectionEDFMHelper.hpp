#ifndef __PROJECTIONEDFMHELPER_H_
#define __PROJECTIONEDFMHELPER_H_

#include "finiteVolume/FluxApproximationBase.hpp"  // FIXME: include appropriate headers
#include <vector>

namespace geosx {

class ProjectionEDFMHelper {
 public:
  ProjectionEDFMHelper( MeshLevel const & mesh,
                        GeometricObjectManager const * geometricObjManager );
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

  bool intersection( arraySlice1d< real64 const > const & fracOrigin,
                     arraySlice1d< real64 const > const & fracNormal,
                     localIndex edgeIdx,
                     real64 (&tmp)[3] ) const noexcept;

  bool isBoundaryFace( localIndex faceIdx ) const noexcept;
  bool onLargerSide( localIndex faceIdx,
                     real64 signedDistanceCellCenterToFrac,
                     arraySlice1d< real64 const > const & fracOrigin,
                     arraySlice1d< real64 const > const & fracNormal ) const noexcept;
  real64 getSignedDistanceCellCenterToFracPlane( CellID const & hostCellID,
                                                 arraySlice1d< real64 const > const & fracNormal,
                                                 arraySlice1d< real64 const > const & fracOrigin,
                                                 real64 (&tmp)[3] ) const noexcept;
  bool neighborOnSameSide( localIndex faceIdx,
                           real64 signedDistanceCellCenterToFrac,
                           CellID const & hostCellID,
                           EmbeddedSurfaceSubRegion const & fractureSubRegion ) const;

  CellID otherCell( localIndex faceIdx, CellID const & hostCellID ) const;

  real64 fractureMatrixTransmissilibility( CellID const & neighborCell,
                                           localIndex fracElement,
                                           EmbeddedSurfaceSubRegion const & fractureSubRegion,
                                           localIndex faceIdx ) const;

  // void addNonNeighboringConnection(  )

  MeshLevel const & m_mesh;
  GeometricObjectManager const * m_geometricObjManager;
  ElementRegionManager const * const m_elementManager;
  FaceManager const * const m_faceManager;
  NodeManager const * const m_nodeManager;
  EdgeManager const * const m_edgeManager;
  // ArrayOfArraysView< localIndex const > const & m_faceToEdges;
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_nodesCoord;
  arrayView2d< localIndex const > const m_edgeToNodes;
  arrayView2d< localIndex const > const m_facesToCells;
  arrayView2d< localIndex const > const m_facesToRegions;
  arrayView2d< localIndex const > const m_facesToSubRegions;
  ArrayOfArraysView< localIndex const > const m_facesToNodes;
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > m_nodeReferencePosition;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const m_cellCenters;
};

}  // end namespace geosx


#endif // __PROJECTIONEDFMHELPER_H_
