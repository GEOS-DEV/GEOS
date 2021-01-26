#include "ProjectionEDFMHelper.hpp"
#include <limits>  // numeric_limits

namespace geosx {

using std::ref;

ProjectionEDFMHelper::ProjectionEDFMHelper( MeshLevel const & mesh,
                                            GeometricObjectManager const * geometricObjManager )
    : m_mesh( mesh ), m_geometricObjManager( geometricObjManager ),
      m_elementManager( mesh.getElemManager() ), m_faceManager( mesh.getFaceManager() ),
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
      // real64 cellCenter[ 3 ];
      // LvArray::tensorOps::copy< 3 >(cellCenter, m_cellCenters[ hostCellRegionIdx ][ hostCellSubRegionIdx ][ hostCellIdx ]);

      // get host cell faces
      CellElementRegion const * cellRegion = m_elementManager->GetRegion< CellElementRegion >( hostCellRegionIdx );
      CellElementSubRegion const * cellSubRegion = cellRegion->GetSubRegion< CellElementSubRegion >(hostCellSubRegionIdx);

      // pick faces for non-neighboring pEDFM connections
      auto const faces = selectFaces(cellSubRegion->faceList(), cellID, fracElement, fractureSubRegion);
      for (localIndex const faceIdx : faces)
      {
        CellID neighborCell = otherCell(faceIdx, cellID);
        real64 transFM = fractureMatrixTransmissilibility( neighborCell, fracElement,
                                                           fractureSubRegion, faceIdx );
        // addNonNeighboringConnection(fracElement, otherCell);
      }

      exit(0);
    }

  }

}

std::vector<localIndex> ProjectionEDFMHelper::selectFaces(FixedOneToManyRelation const & subRegionFaces,
                                                          CellID const & hostCellID,
                                                          localIndex fracElement,
                                                          EmbeddedSurfaceSubRegion const & fractureSubRegion) const
{
  arraySlice1d< real64 const > const n = fractureSubRegion.getNormalVector(fracElement);
  arraySlice1d< real64 const> const o = fractureSubRegion.getOrigin(fracElement);
  real64 tmp[ 3 ];
  real64 const distToFrac = getSignedDistanceCellCenterToFracPlane( hostCellID, n, o, tmp );

  // pick faces that intersect the fracture
  std::vector<localIndex> faces;
  for (localIndex const iface : subRegionFaces[hostCellID.index])
  {
    if (isBoundaryFace(iface)) continue;
    // face intersected by frac and non-branching
    if ( intersection( ref(o), ref(n), iface, tmp ) && neighborOnSameSide( iface, distToFrac,
                                                                           hostCellID,
                                                                           fractureSubRegion) )
    {
      faces.push_back( iface );
      continue;
    }

    if ( onLargerSide( iface, distToFrac, o, n ) )
    {
      faces.push_back( iface );
      continue;
    }
  }

  return faces;
}

bool ProjectionEDFMHelper::intersection( arraySlice1d< real64 const > const & fracOrigin,
                                         arraySlice1d< real64 const > const & fracNormal,
                                         localIndex edgeIdx,
                                         real64 (&tmp)[3] ) const noexcept
{
  real64 signedDistanceProduct = 1.f;
  int const nEdgeVertices = 2;
  for (int ivertex = 0; ivertex < nEdgeVertices; ivertex++)
  {
    LvArray::tensorOps::copy< 3 > ( tmp, m_nodesCoord[m_edgeToNodes[edgeIdx][ivertex]] );
    LvArray::tensorOps::subtract< 3 >( tmp, fracOrigin );
    signedDistanceProduct *= LvArray::tensorOps::AiBi< 3 >( tmp, fracNormal );
  }

  return signedDistanceProduct <= 0;
}

bool ProjectionEDFMHelper::isBoundaryFace(localIndex faceIdx) const noexcept
{
  if( m_facesToCells[faceIdx][0] < 0 || m_facesToCells[faceIdx][1] < 0 )
    return true;
  return false;
}

bool ProjectionEDFMHelper::onLargerSide( localIndex faceIdx,
                                         real64 signedDistanceCellCenterToFrac,
                                         arraySlice1d< real64 const > const & fracOrigin,
                                         arraySlice1d< real64 const > const & fracNormal ) const noexcept
{
  real64 const areaTolerance = 10.f * std::numeric_limits<real64>::epsilon();
  real64 faceCenter[ 3 ], faceNormal[ 3 ];
  // compute face center
  computationalGeometry::Centroid_3DPolygon( m_facesToNodes[faceIdx], m_nodeReferencePosition,
                                               faceCenter, faceNormal, areaTolerance );
  LvArray::tensorOps::subtract< 3 >( faceCenter, fracOrigin );
  return LvArray::tensorOps::AiBi< 3 >( faceCenter, fracNormal ) * signedDistanceCellCenterToFrac > 0 ;
}

real64 ProjectionEDFMHelper::getSignedDistanceCellCenterToFracPlane( CellID const & hostCellID,
                                                                     arraySlice1d< real64 const > const & fracNormal,
                                                                     arraySlice1d< real64 const > const & fracOrigin,
                                                                     real64 (&tmp)[3] ) const noexcept
{
  auto const & cellCenter = m_cellCenters[ hostCellID.region ][ hostCellID.subRegion ][ hostCellID.index ];
  LvArray::tensorOps::copy< 3 >( tmp, cellCenter );
  LvArray::tensorOps::subtract< 3 >( tmp, fracOrigin );
  return LvArray::tensorOps::AiBi< 3 >( tmp, fracNormal );
}

ProjectionEDFMHelper::CellID
ProjectionEDFMHelper::otherCell( localIndex faceIdx, CellID const & hostCellID ) const
{
  const auto& neighbors = m_facesToCells[faceIdx];
  localIndex const ineighbor = (neighbors[0] == hostCellID.index) ? 1 : 0;
  return CellID ( m_facesToRegions[faceIdx][ineighbor],
                  m_facesToSubRegions[faceIdx][ineighbor],
                  m_facesToSubRegions[faceIdx][ineighbor] );
}

bool ProjectionEDFMHelper::neighborOnSameSide( localIndex faceIdx,
                                               real64 signedDistanceCellCenterToFrac,
                                               CellID const & hostCellID,
                                               EmbeddedSurfaceSubRegion const & fractureSubRegion ) const
{
  // find the identification of the neighbor cell
  const auto& neighbors = m_facesToCells[faceIdx];
  localIndex const ineighbor = (neighbors[0] == hostCellID.index) ? 1 : 0;
  CellID neighborCellID = otherCell( faceIdx, hostCellID );

  // get cell center
  real64 cellCenter[ 3 ];
  LvArray::tensorOps::copy< 3 >(cellCenter,
                                m_cellCenters[ neighborCellID.region ]
                                [ neighborCellID.subRegion ][ neighborCellID.index ]);

  // get fracture normal and origin in the neighbor cells:
  // they might be different than those in the host cell
  CellElementRegion const * neighborRegion =
      m_elementManager->GetRegion< CellElementRegion >( neighborCellID.region );
  CellElementSubRegion const * neighborSubRegion =
      neighborRegion->GetSubRegion< CellElementSubRegion >( neighborCellID.subRegion );
  auto const & efracSurfaces = neighborSubRegion->embeddedSurfacesList();

  // embedded fracture elements hosted in the neighbor cell
  // (there could be more than one)
  auto const & efracElements = efracSurfaces[neighborCellID.index];
  // TODO: there could also be zero
  if (efracElements.size() == 0)
    return false;
  if (efracElements.size() > 1)
    throw std::invalid_argument("Write the code to support pEDFM with fracture intersection");

  localIndex const fracElement = efracElements[0];
  arraySlice1d< real64 const > const n = fractureSubRegion.getNormalVector(fracElement);
  arraySlice1d< real64 const> const o = fractureSubRegion.getOrigin(fracElement);
  real64 tmp[ 3 ];
  real64 const signedDistanceNeighborlCenterToFrac =
      getSignedDistanceCellCenterToFracPlane( neighborCellID, n, o, tmp );

  return signedDistanceCellCenterToFrac * signedDistanceNeighborlCenterToFrac > 0;
}

// void ProjectionEDFMHelper::addNonNeighboringConnection( localIndex fracElement,
//                                                         localIndex faceIdx )
// {
//   localIndex const numElems = 2;   // tpfa
//   stackArray1d< localIndex, maxElems > stencilCellsRegionIndex( numElems );
//   stackArray1d< localIndex, maxElems > stencilCellsSubRegionIndex( numElems );
//   stackArray1d< localIndex, maxElems > stencilCellsIndex( numElems );
//   stackArray1d< real64, maxElems > stencilWeights( numElems );

// }

real64 ProjectionEDFMHelper::
fractureMatrixTransmissilibility( CellID const & neighborCell,
                                  localIndex fracElement,
                                  EmbeddedSurfaceSubRegion const & fractureSubRegion,
                                  localIndex faceIdx ) const
{
  // TODO: use m_geometricObjManager to get the nodal coordinates of the fracuture bounded plane
  // m_geometricObjManager.
  // compute face center and normal
  real64 const areaTolerance = 10.f * std::numeric_limits<real64>::epsilon();
  real64 faceCenter[ 3 ], faceNormal[ 3 ];
  computationalGeometry::Centroid_3DPolygon( m_facesToNodes[faceIdx], m_nodeReferencePosition,
                                               faceCenter, faceNormal, areaTolerance );
  // compute the projection of frac element onto the face

  // get cell center
  real64 cellCenter[ 3 ];
  LvArray::tensorOps::copy< 3 >(cellCenter,
                                m_cellCenters[ neighborCell.region ]
                                [ neighborCell.subRegion ]
                                [ neighborCell.index ]);

  // get frac element center


  return 0.f;
}

}  // end namespace geosx
