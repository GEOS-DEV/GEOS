#include "ProjectionEDFMHelper.hpp"

namespace geosx {

ProjectionEDFMHelper::ProjectionEDFMHelper(MeshLevel const & mesh)
    : m_mesh( mesh ), m_elementManager( mesh.getElemManager() ), m_faceManager( mesh.getFaceManager() ),
      m_nodeManager( mesh.getNodeManager() ), m_edgeManager( mesh.getEdgeManager() ),
      m_nodesCoord( m_nodeManager->referencePosition() ),
      m_edgeToNodes( m_edgeManager->nodeList() ),
      m_facesToCells( m_faceManager->elementList() ),
      m_facesToNodes( m_faceManager->nodeList().toViewConst() ),
      m_nodeReferencePosition( m_nodeManager->referencePosition() ),
      m_cellCenters( m_elementManager->ConstructArrayViewAccessor< real64, 2 >( CellBlock::viewKeyStruct::elementCenterString ) )
      // m_faceToEdges(m_faceManager->edgeList().toViewConst())
{}

void ProjectionEDFMHelper::addNonNeighboringConnections(EmbeddedSurfaceSubRegion const & fractureSubRegion)
{
  arrayView1d< integer const > const ghostRank = fractureSubRegion.ghostRank();

  // ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const cellCenters =
  //   m_elementManager->ConstructArrayViewAccessor< real64, 2 >( CellBlock::viewKeyStruct::elementCenterString );
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
      CellID cellID( hostCellRegionIdx, hostCellSubRegionIdx, hostCellIdx );

      // get cell center
      real64 cellCenter[ 3 ];
      LvArray::tensorOps::copy< 3 >(cellCenter, m_cellCenters[ hostCellRegionIdx ][ hostCellSubRegionIdx ][ hostCellIdx ]);

      // get host cell faces
      CellElementRegion const * cellRegion = m_elementManager->GetRegion< CellElementRegion >( hostCellRegionIdx );
      CellElementSubRegion const * cellSubRegion = cellRegion->GetSubRegion< CellElementSubRegion >(hostCellSubRegionIdx);

      auto const faces = selectFaces(cellSubRegion->faceList(), cellID, fracElement, fractureSubRegion);
      exit(0);
    }

  }

}

std::vector<localIndex> ProjectionEDFMHelper::selectFaces(FixedOneToManyRelation const & subRegionFaces,
                                                          CellID const & hostCellID,
                                                          localIndex fracElement,
                                                          EmbeddedSurfaceSubRegion const & fractureSubRegion) const
{
  using std::ref;
  auto const & n = fractureSubRegion.getNormalVector(fracElement);
  auto const & o = fractureSubRegion.getOrigin(fracElement);
  std::vector<localIndex> faces;
  R1Tensor tmp;
  real64 const distToFrac = getSingedDistanceCellCenterToFracPlane( hostCellID, n, o, tmp );

  // pick faces that intersect the fracture
  for (localIndex const iface : subRegionFaces[hostCellID.index])
  {
    if (isBoundaryFace(iface)) continue;
    // count those that are intersectedd
    if ( intersection( ref(o), ref(n), iface, tmp ) )
    {
      faces.push_back( iface );
      continue;
      // TODO: do smart branch elimination
    }

    if ( onLargerSide( iface, distToFrac ) )
    {
      faces.push_back( iface );
      continue;
    }

    // get face properties
    real64 faceCenter[ 3 ], faceNormal[ 3 ], faceConormal[ 3 ], cellToFaceVec[ 3 ];
    // real64 const lengthTolerance = m_lengthScale * m_areaRelTol;
    // real64 const areaTolerance = lengthTolerance * lengthTolerance;
    // computationalGeometry::Centroid_3DPolygon( m_facesToNodes[iface], m_nodeReferencePosition,
                                               // faceCenter, faceNormal, areaTolerance );
  }

  return faces;
}

bool ProjectionEDFMHelper::intersection(R1Tensor const & fracOrigin,
                                        R1Tensor const & fracNormal,
                                        localIndex edgeIdx,
                                        R1Tensor & tmp) const noexcept
{
  int const nEdgeVertices = 2;
  double signedDistanceProduct = 1;
  for (int ivertex = 0; ivertex < nEdgeVertices; ivertex++)
  {
    LvArray::tensorOps::copy< 3 > ( tmp, m_nodesCoord[m_edgeToNodes[edgeIdx][ivertex]] );
    tmp -= fracOrigin;
    signedDistanceProduct *= Dot(tmp, fracNormal);
  }

  return signedDistanceProduct <= 0;
}

bool ProjectionEDFMHelper::isBoundaryFace(localIndex faceIdx) const noexcept
{
  if( m_facesToCells[faceIdx][0] < 0 || m_facesToCells[faceIdx][1] < 0 )
    return true;
  return false;
}

bool ProjectionEDFMHelper::onLargerSide(localIndex faceIdx,
                                        real64 signedDistanceCellCenterToFrac) const noexcept
{

  return true;
}

real64 ProjectionEDFMHelper::getSingedDistanceCellCenterToFracPlane( CellID const & hostCellID,
                                                                     R1Tensor const & fracNormal,
                                                                     R1Tensor const & fracOrigin,
                                                                     R1Tensor & tmp) const noexcept
{
  auto const & cellCenter = m_cellCenters[ hostCellID.region ][ hostCellID.subRegion ][ hostCellID.index ];
  LvArray::tensorOps::copy< 3 >( tmp, cellCenter );
  LvArray::tensorOps::subtract< 3 >( tmp, fracOrigin );
  return Dot( tmp, fracOrigin );
}

}  // end namespace geosx
