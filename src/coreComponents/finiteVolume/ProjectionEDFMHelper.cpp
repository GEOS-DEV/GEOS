#include "ProjectionEDFMHelper.hpp"

namespace geosx {

ProjectionEDFMHelper::ProjectionEDFMHelper(MeshLevel const & mesh)
    : m_mesh( mesh ), m_elementManager( mesh.getElemManager() ), m_faceManager( mesh.getFaceManager() ),
      m_nodeManager( mesh.getNodeManager() ), m_edgeManager( mesh.getEdgeManager() )
{}

void ProjectionEDFMHelper::addNonNeighboringConnections(EmbeddedSurfaceSubRegion const & fractureSubRegion)
{
  arrayView1d< integer const > const ghostRank = fractureSubRegion.ghostRank();

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const cellCenters =
    m_elementManager->ConstructArrayViewAccessor< real64, 2 >( CellBlock::viewKeyStruct::elementCenterString );
  // arrayView2d< real64 const > const fractureElemCenter = fractureSubRegion.getElementCenter().toViewConst();

  FixedToManyElementRelation const & surfaceElementsToCells = fractureSubRegion.getToCellRelation();

  for( localIndex fracElement = 0; fracElement < fractureSubRegion.size(); fracElement++ )
  {
    if( ghostRank[fracElement] < 0 )
    {
      // get host cell information
      localIndex const hostCellRegionIdx  = surfaceElementsToCells.m_toElementRegion[fracElement][0];
      localIndex const hostCellSubRegionIdx = surfaceElementsToCells.m_toElementSubRegion[fracElement][0];
      localIndex const hostCellIdx = surfaceElementsToCells.m_toElementIndex[fracElement][0];

      // get cell center
      real64 cellCenter[ 3 ];
      LvArray::tensorOps::copy< 3 >(cellCenter, cellCenters[ hostCellRegionIdx ][ hostCellSubRegionIdx ][ hostCellIdx ]);

      // get host cell faces
      CellElementRegion const * cellRegion = m_elementManager->GetRegion< CellElementRegion >( hostCellRegionIdx );
      CellElementSubRegion const * cellSubRegion = cellRegion->GetSubRegion< CellElementSubRegion >(hostCellSubRegionIdx);
      // auto const & hostCellFaces = cellSubRegion->faceList()[ hostCellIdx ];
      // std::cout << "hostCellFaces = " << typeid(hostCellFaces).name() << std::endl;
      exit(0);
      // auto const & hostCellFaces = elemManager.GetRegion< CellElementRegion >( hostCellRegionIdx )->
      //                                          GetSubRegion< CellElementSubRegion >(hostCellSubRegionIdx)->
      //                                          faceList()[ hostCellIdx ];

      selectFaces(cellSubRegion->faceList(), hostCellIdx, fracElement, fractureSubRegion);
      // for( localIndex const iface : hostCellFaces )
      // {

      // }

      // for ()
      // const auto isolating_faces = select_faces_(*frac_face);

      // get face properties
      // real64 faceCenter[ 3 ], faceNormal[ 3 ], faceConormal[ 3 ], cellToFaceVec[ 3 ];
      // computationalGeometry::Centroid_3DPolygon( faceToNodes[kf], X, faceCenter, faceNormal, areaTolerance );


    }

  }

}

void ProjectionEDFMHelper::selectFaces(FixedOneToManyRelation const & subRegionFaces,
                                       localIndex hostCellIdx,
                                       localIndex fracElement,
                                       EmbeddedSurfaceSubRegion const & fractureSubRegion) const
{
  using std::ref;
  ArrayOfArraysView< localIndex const > const & faceToEdges = m_faceManager->edgeList().toViewConst();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoord =
      m_nodeManager->referencePosition();
  arrayView2d< localIndex const > const edgeToNodes = m_edgeManager->nodeList();
  auto const & n = fractureSubRegion.getNormalVector(fracElement);
  auto const & o = fractureSubRegion.getOrigin(fracElement);
  vector<localIndex> faces;
  R1Tensor vertex1, vertex2;
  for (localIndex const iface : subRegionFaces[hostCellIdx])
  {

    for ( localIndex const iEdge : faceToEdges[iface] )
    {
      LvArray::tensorOps::copy< 3 >
          ( vertex1, nodesCoord[edgeToNodes[iEdge][0]] );
      LvArray::tensorOps::copy< 3 >
          ( vertex2, nodesCoord[edgeToNodes[iEdge][1]] );

      if ( intersection(ref(vertex1), ref(vertex2), ref(n), ref(o)) )
      {
        faces.push_back( iface );
        break;
      }
    }
    // if (intersects(iface))
    // for
    // for (localIndex edgeIdx = 0; edgeIdx < )
  }

}

bool ProjectionEDFMHelper::intersection(R1Tensor const & fracOrigin,
                                        R1Tensor const & fracNormal,
                                        R1Tensor const & vertex1,
                                        R1Tensor const & vertex2) const noexcept
{
  R1Tensor ov1, ov2;
  ov1 = vertex1;
  ov1 -= fracOrigin;
  ov2 = vertex1;
  ov2 -= fracOrigin;
  return Dot(ov1, fracNormal) * Dot(ov2, fracNormal) <= 0;
}

}  // end namespace geosx
