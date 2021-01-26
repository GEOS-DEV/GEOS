#include "ProjectionEDFMHelper.hpp"
#include <limits>  // numeric_limits

namespace geosx {

using std::ref;

ProjectionEDFMHelper::ProjectionEDFMHelper( MeshLevel const & mesh,
                                            GeometricObjectManager const * geometricObjManager,
                                            std::string const & coeffName,
                                            CellElementStencilTPFA & stencil )
    : m_mesh( mesh ), m_geometricObjManager( geometricObjManager ),
      m_elementManager( mesh.getElemManager() ), m_faceManager( mesh.getFaceManager() ),
      m_nodeManager( mesh.getNodeManager() ), m_edgeManager( mesh.getEdgeManager() ),
      m_nodesCoord( m_nodeManager->referencePosition() ),
      m_edgeToNodes( m_edgeManager->nodeList() ),
      m_facesToCells( m_faceManager->elementList() ),
      m_facesToNodes( m_faceManager->nodeList().toViewConst() ),
      m_nodeReferencePosition( m_nodeManager->referencePosition() ),
      m_cellCenters( m_elementManager->ConstructArrayViewAccessor< real64, 2 >( CellBlock::viewKeyStruct::elementCenterString ) ),
      m_permTensor( m_elementManager->ConstructArrayViewAccessor< real64, 2 >( coeffName ) ),
      m_stencil( stencil )
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


        addNonNeighboringConnection(fracElement, neighborCell, transFM, fractureSubRegion);
        // I assume that the m-m connection id is equal to the face id
        m_stencil.zero( faceIdx );
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
  // arraySlice1d< real64 const> const o = fractureSubRegion.getOrigin(fracElement);
  arrayView2d< real64 const > const & centers = fractureSubRegion.getElementCenter().toViewConst();
  real64 origin[3];
  LvArray::tensorOps::copy< 3 >( origin, centers[ fracElement ] );

  real64 tmp[ 3 ];
  real64 const distToFrac = getSignedDistanceCellCenterToFracPlane( hostCellID, n, origin, tmp );

  // pick faces that intersect the fracture
  std::vector<localIndex> faces;
  for (localIndex const iface : subRegionFaces[hostCellID.index])
  {
    if (isBoundaryFace(iface)) continue;
    // face intersected by frac and non-branching
    if ( intersection( ref(origin), ref(n), iface, tmp ) && neighborOnSameSide( iface, distToFrac,
                                                                           hostCellID,
                                                                           fractureSubRegion) )
    {
      faces.push_back( iface );
      continue;
    }

    if ( onLargerSide( iface, distToFrac, ref(origin), ref(n) ) )
    {
      faces.push_back( iface );
      continue;
    }
  }

  return faces;
}

bool ProjectionEDFMHelper::intersection( real64 (&fracOrigin)[3],
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
                                         real64 (&fracOrigin)[3],
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
                                                                     real64 const (&fracOrigin)[3],
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
  // arraySlice1d< real64 const> const o = fractureSubRegion.getOrigin(fracElement);
  arrayView2d< real64 const > const & centers = fractureSubRegion.getElementCenter().toViewConst();
  real64 origin[3];
  LvArray::tensorOps::copy< 3 >( origin, centers[ fracElement ] );

  // const auto o = [fracElement];
  real64 tmp[ 3 ];
  real64 const signedDistanceNeighborlCenterToFrac =
      getSignedDistanceCellCenterToFracPlane( neighborCellID, n, origin, tmp );

  return signedDistanceCellCenterToFrac * signedDistanceNeighborlCenterToFrac > 0;
}

void ProjectionEDFMHelper::addNonNeighboringConnection( localIndex fracElement,
                                                        CellID const & cell,
                                                        real64 transmissibility,
                                                        EmbeddedSurfaceSubRegion const & fractureSubRegion )
{
  localIndex constexpr maxElems = FaceElementStencil::MAX_STENCIL_SIZE;
  localIndex constexpr numElems = 2;   // tpfa
  stackArray1d< localIndex, maxElems > stencilCellsRegionIndex( numElems );
  stackArray1d< localIndex, maxElems > stencilCellsSubRegionIndex( numElems );
  stackArray1d< localIndex, maxElems > stencilCellsIndex( numElems );
  stackArray1d< real64, maxElems > stencilWeights( numElems );

  stencilCellsRegionIndex[0] = cell.region;
  stencilCellsSubRegionIndex[0] = cell.subRegion;
  stencilCellsIndex[0] = cell.index;
  stencilWeights[0] =  transmissibility;

  /* auto * fractureRegion = fractureSubRegion.getParent(); */
  /* localIndex const fractureRegionIndex = fractureRegion->getIndexInParent(); */

  stencilCellsRegionIndex[1] = fractureSubRegion.getParent()->getIndexInParent();
  stencilCellsSubRegionIndex[1] = fractureSubRegion.getIndexInParent();
  stencilCellsIndex[1] = fracElement;
  stencilWeights[1] = -transmissibility;

  m_stencil.add(  2,
                  stencilCellsRegionIndex.data(),
                  stencilCellsSubRegionIndex.data(),
                  stencilCellsIndex.data(),
                  stencilWeights.data(),
                  m_stencil.size() );
}

real64 ProjectionEDFMHelper::
fractureMatrixTransmissilibility( CellID const & neighborCell,
                                  localIndex fracElement,
                                  EmbeddedSurfaceSubRegion const & fractureSubRegion,
                                  localIndex faceIdx ) const
{

  // TODO: should I really compute the real projection, or assuming the whole face is occupied by the projection?
  // std::string const & fractureName = fractureSubRegion.getFractureName(fracElement);
  // const BoundedPlane & plane = m_geometricObjManager->getReference< BoundedPlane >( fractureName );
  // array2d< real64 > & boundingPoints = plane

  // compute face center and normal
  real64 const areaTolerance = 10.f * std::numeric_limits<real64>::epsilon();
  real64 faceCenter[ 3 ], faceNormal[ 3 ];
  real64 faceArea = computationalGeometry::Centroid_3DPolygon( m_facesToNodes[faceIdx], m_nodeReferencePosition,
                                                               faceCenter, faceNormal, areaTolerance );

  // get cell center
  real64 cellCenter[ 3 ];
  LvArray::tensorOps::copy< 3 >(cellCenter,
                                m_cellCenters[ neighborCell.region ]
                                [ neighborCell.subRegion ]
                                [ neighborCell.index ]);

  real64 fracCenter[ 3 ];
  arrayView2d< real64 const > const & centers = fractureSubRegion.getElementCenter().toViewConst();
  LvArray::tensorOps::copy< 3 >( fracCenter, centers[ fracElement ] );


  // projection point
  // const double t =  (cf - c1).dot(n) / (c2 - c1).dot(n);
  real64 tmp[3];
  LvArray::tensorOps::copy< 3 >( tmp, faceCenter );
  LvArray::tensorOps::subtract< 3 >( tmp, fracCenter );
  real64 const numerator = LvArray::tensorOps::AiBi< 3 >( tmp, faceNormal );
  LvArray::tensorOps::copy< 3 >( tmp, cellCenter );
  LvArray::tensorOps::subtract< 3 >( tmp, fracCenter );
  real64 const denom = LvArray::tensorOps::AiBi< 3 >( tmp, faceNormal );
  real64 t = numerator / denom;

  /* const auto &c1 = frac.center; */
  /* const auto &c2 = cell.center; */
  /* const auto &cf = con.center; */
  /* const auto &n = con.normal; */
  /* const auto cp = c1 + t*(c2 - c1); */
  LvArray::tensorOps::copy< 3 >( tmp, cellCenter );
  LvArray::tensorOps::subtract< 3 >( tmp, fracCenter );
  real64 projectionPoint[3];
  LvArray::tensorOps::copy< 3 >( tmp, fracCenter );
  LvArray::tensorOps::scaledAdd< 3 >( projectionPoint, fracCenter, t );

  // fracture projected permeability: equal to one for now
  real64 const directionalPermFracutre = 1.f;
  // matrix permeability
  auto const neighborPerm =  m_permTensor[neighborCell.region][neighborCell.subRegion][neighborCell.index];
  /* const double Kp2 = (K2 * (c2 - cp).normalize()).norm(); */
  real64 direction[3];
  LvArray::tensorOps::copy< 3 >( direction, cellCenter );
  LvArray::tensorOps::subtract< 3 >( direction, projectionPoint );
  LvArray::tensorOps::normalize< 3 >(direction);
  LvArray::tensorOps::Ri_eq_symAijBj< 3 >( tmp, neighborPerm, direction );
  real64 const directionalPermCell = LvArray::tensorOps::l2Norm< 3 >( tmp );

  // part trasmissibilities
  /* const double T1 = face_area * Kp1 / (c1 - cp).norm(); */
  LvArray::tensorOps::copy< 3 >( tmp, fracCenter );
  LvArray::tensorOps::subtract< 3 >( tmp, projectionPoint );
  real64 const t1 = faceArea * directionalPermFracutre * LvArray::tensorOps::l2Norm< 3 >( tmp );

  /* const double T2 = face_area * Kp2 / (c2 - cp).norm(); */
  LvArray::tensorOps::copy< 3 >( tmp, cellCenter );
  LvArray::tensorOps::subtract< 3 >( tmp, projectionPoint );
  real64 const t2 = faceArea * directionalPermCell * LvArray::tensorOps::l2Norm< 3 >( tmp );

  real64 trans = 0.f;
  if (std::isnormal( t1 + t2 ))
    trans = t1 * t2 / ( t1 + t2 );

  return trans;
}

}  // end namespace geosx
