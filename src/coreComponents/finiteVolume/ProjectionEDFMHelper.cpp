#include "ProjectionEDFMHelper.hpp"

namespace geosx {

ProjectionEDFMHelper::ProjectionEDFMHelper(MeshLevel const & mesh)
    : m_mesh( mesh ), m_elementManager( mesh.getElemManager() ), m_faceManager( mesh.getFaceManager() ),
      m_nodeManager( mesh.getNodeManager() ), m_edgeManager( mesh.getEdgeManager() ),
      m_nodesCoord( m_nodeManager->referencePosition() ),
      m_edgeToNodes( m_edgeManager->nodeList() ),
      m_facesToCells( m_faceManager->elementList() )
      // m_faceToEdges(m_faceManager->edgeList().toViewConst())
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

      auto const faces = selectFaces(cellSubRegion->faceList(), hostCellIdx, fracElement, fractureSubRegion);
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

std::vector<localIndex> ProjectionEDFMHelper::selectFaces(FixedOneToManyRelation const & subRegionFaces,
                                                          localIndex hostCellIdx,
                                                          localIndex fracElement,
                                                          EmbeddedSurfaceSubRegion const & fractureSubRegion) const
{
  using std::ref;
  // arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoord =
  //     m_nodeManager->referencePosition();
  // arrayView2d< localIndex const > const edgeToNodes = m_edgeManager->nodeList();
  auto const & n = fractureSubRegion.getNormalVector(fracElement);
  auto const & o = fractureSubRegion.getOrigin(fracElement);
  std::vector<localIndex> faces;
  R1Tensor tmp;
  // pick faces that intersect the fracture
  for (localIndex const iface : subRegionFaces[hostCellIdx])
  {
    if (isBoundaryFace(iface)) continue;
    if ( intersection( ref(o), ref(n), iface, tmp ) )
      faces.push_back( iface );

    // for ( localIndex const iEdge : faceToEdges[iface] )
    // {
    //   LvArray::tensorOps::copy< 3 >
    //       ( vertex1, nodesCoord[edgeToNodes[iEdge][0]] );
    //   LvArray::tensorOps::copy< 3 >
    //       ( vertex2, nodesCoord[edgeToNodes[iEdge][1]] );

    //   if ( intersection(ref(vertex1), ref(vertex2), ref(n), ref(o)) )
    //   {
    //     faces.push_back( iface );
    //     break;
    //   }
    // }
  }

  // pick faces that are on the larger side of the fracure
  for (localIndex const iface : subRegionFaces[hostCellIdx])
    if (std::find( faces.begin(), faces.end(), iface ) == faces.end())  // not intersection
    {

    }

  return faces;
}

bool ProjectionEDFMHelper::intersection(R1Tensor const & fracOrigin,
                                        R1Tensor const & fracNormal,
                                        localIndex edgeIdx,
                                        R1Tensor & tmp) const noexcept
{
  int const nEdgeVertices = 2;
  double sameSide = 1;
  for (int ivertex = 0; ivertex < nEdgeVertices; ivertex++)
  {
    LvArray::tensorOps::copy< 3 > ( tmp, m_nodesCoord[m_edgeToNodes[edgeIdx][ivertex]] );
    tmp -= fracOrigin;
    sameSide *= Dot(tmp, fracNormal);
  }

  return sameSide <= 0;
}

bool ProjectionEDFMHelper::isBoundaryFace(localIndex faceIdx) const noexcept
{
  if( m_facesToCells[faceIdx][0] < 0 || m_facesToCells[faceIdx][1] < 0 )
    return true;
  return false;
}

}  // end namespace geosx
