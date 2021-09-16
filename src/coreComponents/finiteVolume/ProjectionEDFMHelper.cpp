#include "ProjectionEDFMHelper.hpp"
#include "mesh/MeshLevel.hpp"
#include "finiteVolume/CellElementStencilTPFA.hpp"
#include "finiteVolume/SurfaceElementStencil.hpp"
#include "finiteVolume/EmbeddedSurfaceToCellStencil.hpp"

#include <limits>  // numeric_limits

namespace geosx
{

using std::ref;

ProjectionEDFMHelper::ProjectionEDFMHelper( MeshLevel const & mesh,
                                            CellElementStencilTPFA & cellStencil,
                                            EmbeddedSurfaceToCellStencil & edfmStencil,
                                            string const & embeddedSurfaceRegionName )
  : m_elementManager( mesh.getElemManager() ),
  m_nodesCoord( mesh.getNodeManager().referencePosition() ),
  m_edgeToNodes( mesh.getEdgeManager().nodeList() ),
  m_facesToCells( mesh.getFaceManager().elementList() ),
  m_facesToRegions( mesh.getFaceManager().elementRegionList() ),
  m_facesToSubRegions( mesh.getFaceManager().elementSubRegionList() ),
  m_facesToNodes( mesh.getFaceManager().nodeList().toViewConst() ),
  m_cellCenters( m_elementManager.constructArrayViewAccessor< real64, 2 >( CellBlock::viewKeyStruct::elementCenterString() ) ),
  m_cellStencil( cellStencil ),
  m_edfmStencil( edfmStencil ),
  m_embeddedSurfaceRegionName( embeddedSurfaceRegionName )
{}

void ProjectionEDFMHelper::addNonNeighboringConnections() const
{
  SurfaceElementRegion const & fractureRegion = m_elementManager.getRegion< SurfaceElementRegion >( m_embeddedSurfaceRegionName );
  EmbeddedSurfaceSubRegion const & fractureSubRegion = fractureRegion.getSubRegion< EmbeddedSurfaceSubRegion >( "embeddedSurfaceSubRegion" );

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
      CellDescriptor cellID( hostCellRegionIdx, hostCellSubRegionIdx, hostCellIdx );

      // get host cell faces
      CellElementRegion const & cellRegion = m_elementManager.getRegion< CellElementRegion >( hostCellRegionIdx );
      CellElementSubRegion const & cellSubRegion = cellRegion.getSubRegion< CellElementSubRegion >( hostCellSubRegionIdx );

      // pick faces for non-neighboring pEDFM connections
      std::list< localIndex > const faces = selectFaces( cellSubRegion.faceList(), cellID, fracElement, fractureSubRegion );
      for( localIndex const faceIdx : faces )
      {
        CellDescriptor neighborCell = otherCell( faceIdx, cellID );
        real64 transFM[2];
        fractureMatrixTransmissilibility( neighborCell, fracElement, fractureSubRegion, faceIdx, transFM );

        addNonNeighboringConnection( fracElement, neighborCell, transFM, fractureSubRegion );

        // zero out matrix-matrix connections that are replaced by fracture-matrix connections
        // I assume that the m-m connection index is equal to the face id
        m_cellStencil.zero( faceIdx );
      }
    }
  }
}

std::list< localIndex > ProjectionEDFMHelper::selectFaces( FixedOneToManyRelation const & subRegionFaces,
                                                           CellDescriptor const & hostCellID,
                                                           localIndex const fracElement,
                                                           EmbeddedSurfaceSubRegion const & fractureSubRegion ) const
{
  arraySlice1d< real64 const > const n = fractureSubRegion.getNormalVector( fracElement );
  arrayView2d< real64 const > const & centers = fractureSubRegion.getElementCenter().toViewConst();
  real64 origin[3];
  LvArray::tensorOps::copy< 3 >( origin, centers[ fracElement ] );

  real64 tmp[ 3 ];
  real64 distToFrac = getSignedDistanceCellCenterToFracPlane( hostCellID, n, origin, tmp );

  // pick faces that intersect the fracture
  std::list< localIndex > faces;
  for( localIndex const iface : subRegionFaces[hostCellID.index] )
  {
    if( isBoundaryFace( iface ))
      continue;
    // face intersected by frac and non-branching
    if( intersection(origin,  n, iface, tmp) && !neighborOnSameSide(iface, distToFrac, hostCellID, fractureSubRegion) )
    {
      faces.push_back( iface );
      continue;
    }

    // if the frac is horizontal this can be 0 so we just pick one side.
    if (distToFrac < std::numeric_limits< real64 >::epsilon() )
      distToFrac = 1.0;

    if( onLargerSide(iface, distToFrac, origin, n) )
    {
      faces.push_back( iface );
      continue;
    }
  }

  return faces;
}

bool ProjectionEDFMHelper::intersection( real64 const (& fracOrigin)[3],
                                         arraySlice1d< real64 const > const & fracNormal,
                                         localIndex const edgeIdx,
                                         real64 (& tmp)[3] ) const
{
  real64 signedDistanceProduct = 1.f;
  int const nEdgeVertices = 2;
  for( int ivertex = 0; ivertex < nEdgeVertices; ivertex++ )
  {
    LvArray::tensorOps::copy< 3 >( tmp, m_nodesCoord[m_edgeToNodes[edgeIdx][ivertex]] );
    LvArray::tensorOps::subtract< 3 >( tmp, fracOrigin );
    signedDistanceProduct *= LvArray::tensorOps::AiBi< 3 >( tmp, fracNormal );
  }

  return signedDistanceProduct <= 0;
}

bool ProjectionEDFMHelper::isBoundaryFace( localIndex const faceIdx ) const
{
  if( m_facesToCells[faceIdx][0] < 0 || m_facesToCells[faceIdx][1] < 0 )
    return true;
  return false;
}

bool ProjectionEDFMHelper::onLargerSide( localIndex const faceIdx,
                                         real64 const signedDistanceCellCenterToFrac,
                                         real64 const (& fracOrigin)[3],
                                         arraySlice1d< real64 const > const & fracNormal ) const
{
  /* Check If face is on the same side of the fracture as the mass center of the cell. */
  real64 const areaTolerance = 10.f * std::numeric_limits< real64 >::epsilon();
  real64 faceCenter[ 3 ], faceNormal[ 3 ];
  // compute face center
  computationalGeometry::Centroid_3DPolygon( m_facesToNodes[faceIdx], m_nodesCoord,
                                             faceCenter, faceNormal, areaTolerance );
  LvArray::tensorOps::subtract< 3 >( faceCenter, fracOrigin );
  return LvArray::tensorOps::AiBi< 3 >( faceCenter, fracNormal ) * signedDistanceCellCenterToFrac > 0;
}

real64 ProjectionEDFMHelper::getSignedDistanceCellCenterToFracPlane( CellDescriptor const & hostCellID,
                                                                     arraySlice1d< real64 const > const & fracNormal,
                                                                     real64 const (&fracOrigin)[3],
                                                                     real64 (& tmp)[3] ) const
{
  auto const & cellCenter = m_cellCenters[ hostCellID.region ][ hostCellID.subRegion ][ hostCellID.index ];
  LvArray::tensorOps::copy< 3 >( tmp, cellCenter );
  LvArray::tensorOps::subtract< 3 >( tmp, fracOrigin );
  return LvArray::tensorOps::AiBi< 3 >( tmp, fracNormal );
}

CellDescriptor ProjectionEDFMHelper::otherCell( localIndex const faceIdx, CellDescriptor const & hostCellID ) const
{
  const auto & neighbors = m_facesToCells[faceIdx];
  localIndex const ineighbor = (neighbors[0] == hostCellID.index) ? 1 : 0;
  return CellDescriptor ( m_facesToRegions[faceIdx][ineighbor],
                          m_facesToSubRegions[faceIdx][ineighbor],
                          m_facesToCells[faceIdx][ineighbor] );
}

bool ProjectionEDFMHelper::neighborOnSameSide( localIndex const faceIdx,
                                               real64 const signedDistanceCellCenterToFrac,
                                               CellDescriptor const & hostCellID,
                                               EmbeddedSurfaceSubRegion const & fractureSubRegion ) const
{
  // find the identification of the neighbor cell
  CellDescriptor neighborCellID = otherCell( faceIdx, hostCellID );

  // get cell center
  real64 cellCenter[ 3 ];
  LvArray::tensorOps::copy< 3 >( cellCenter,
                                 m_cellCenters[ neighborCellID.region ]
                                 [ neighborCellID.subRegion ][ neighborCellID.index ] );

  // get fracture normal and origin in the neighbor cells:
  // they might be different than those in the host cell
  CellElementRegion const & neighborRegion =
    m_elementManager.getRegion< CellElementRegion >( neighborCellID.region );
  CellElementSubRegion const & neighborSubRegion =
    neighborRegion.getSubRegion< CellElementSubRegion >( neighborCellID.subRegion );
  auto const & efracSurfaces = neighborSubRegion.embeddedSurfacesList();

  // embedded fracture elements hosted in the neighbor cell
  // (there could be more than one)
  auto const & efracElements = efracSurfaces[neighborCellID.index];
  if( efracElements.size() == 0 )
    return false;
  if( efracElements.size() > 1 )
    GEOSX_ERROR( "pEDFM with fracture intersections is not supported yet." );

  localIndex const fracElement = efracElements[0];
  arraySlice1d< real64 const > const n = fractureSubRegion.getNormalVector( fracElement );
  arrayView2d< real64 const > const & centers = fractureSubRegion.getElementCenter().toViewConst();
  real64 origin[3];
  LvArray::tensorOps::copy< 3 >( origin, centers[ fracElement ] );

  // const auto o = [fracElement];
  real64 tmp[ 3 ];
  real64 const signedDistanceNeighborlCenterToFrac =
    getSignedDistanceCellCenterToFracPlane( neighborCellID, n, origin, tmp );

  return signedDistanceCellCenterToFrac * signedDistanceNeighborlCenterToFrac > 0;
}

void ProjectionEDFMHelper::addNonNeighboringConnection( localIndex const fracElement,
                                                        CellDescriptor const & cell,
                                                        real64 const (& transmissibility)[2],
                                                        EmbeddedSurfaceSubRegion const & fractureSubRegion ) const
{
  localIndex constexpr maxElems = EmbeddedSurfaceToCellStencil::MAX_STENCIL_SIZE;
  localIndex const numElems = 2;

  stackArray1d< localIndex, maxElems > stencilCellsRegionIndex( 2 );
  stackArray1d< localIndex, maxElems > stencilCellsSubRegionIndex( 2 );
  stackArray1d< localIndex, maxElems > stencilCellsIndex( 2 );
  stackArray1d< real64, maxElems > stencilWeights( 2 );

  // cell data
  stencilCellsRegionIndex[0] = cell.region;
  stencilCellsSubRegionIndex[0] = cell.subRegion;
  stencilCellsIndex[0] = cell.index;
  stencilWeights[0] =  transmissibility[0];

  // fracture data
  stencilCellsRegionIndex[1] = fractureSubRegion.getParent().getParent().getIndexInParent();
  stencilCellsSubRegionIndex[1] = fractureSubRegion.getIndexInParent();
  stencilCellsIndex[1] = fracElement;
  stencilWeights[1] = transmissibility[1];

  m_edfmStencil.add( numElems,
                     stencilCellsRegionIndex.data(),
                     stencilCellsSubRegionIndex.data(),
                     stencilCellsIndex.data(),
                     stencilWeights.data(),
                     m_edfmStencil.size() );
}

 void ProjectionEDFMHelper::
  fractureMatrixTransmissilibility( CellDescriptor const & neighborCell,
                                    localIndex const fracElement,
                                    EmbeddedSurfaceSubRegion const & fractureSubRegion,
                                    localIndex const faceIdx,
                                    real64 (& trans)[2] ) const
{

  // TODO: should I really compute the real projection, or assuming the whole face is occupied by the projection?

  // compute face center and normal
  real64 const areaTolerance = 10.f * std::numeric_limits< real64 >::epsilon();
  real64 faceCenter[ 3 ], faceNormal[ 3 ];
  real64 faceArea = computationalGeometry::Centroid_3DPolygon( m_facesToNodes[faceIdx], m_nodesCoord,
                                                               faceCenter, faceNormal, areaTolerance );

  // get cell center
  real64 cellCenter[ 3 ];
  LvArray::tensorOps::copy< 3 >( cellCenter,
                                 m_cellCenters[ neighborCell.region ]
                                 [ neighborCell.subRegion ]
                                 [ neighborCell.index ] );

  // get frac  center
  real64 fracCenter[ 3 ];
  arrayView2d< real64 const > const & centers = fractureSubRegion.getElementCenter().toViewConst();
  LvArray::tensorOps::copy< 3 >( fracCenter, centers[ fracElement ] );

  // Compute projection point: a point  on the face that resides  on a  line
  // connecting edfm element and a neighbor cell
  real64 tmp[3];
  LvArray::tensorOps::copy< 3 >( tmp, faceCenter );
  LvArray::tensorOps::subtract< 3 >( tmp, fracCenter );
  real64 const numerator = LvArray::tensorOps::AiBi< 3 >( tmp, faceNormal );
  LvArray::tensorOps::copy< 3 >( tmp, cellCenter );
  LvArray::tensorOps::subtract< 3 >( tmp, fracCenter );
  real64 const denom = LvArray::tensorOps::AiBi< 3 >( tmp, faceNormal );
  real64 t = numerator / denom;

  LvArray::tensorOps::copy< 3 >( tmp, cellCenter );
  LvArray::tensorOps::subtract< 3 >( tmp, fracCenter );
  real64 projectionPoint[3];
  LvArray::tensorOps::copy< 3 >( tmp, fracCenter );
  LvArray::tensorOps::scaledAdd< 3 >( projectionPoint, fracCenter, t );

  // part geometric trasmissibilities
  LvArray::tensorOps::copy< 3 >( tmp, cellCenter );
  LvArray::tensorOps::subtract< 3 >( tmp, projectionPoint );
  trans[0] = faceArea * LvArray::tensorOps::l2Norm< 3 >( tmp );

  LvArray::tensorOps::copy< 3 >( tmp, fracCenter );
  LvArray::tensorOps::subtract< 3 >( tmp, projectionPoint );
  trans[1] = faceArea * LvArray::tensorOps::l2Norm< 3 >( tmp );

}

}  // end namespace geosx
