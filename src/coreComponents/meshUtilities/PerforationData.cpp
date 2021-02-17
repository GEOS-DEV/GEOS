/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * @file PerforationData.cpp
 */

#include "PerforationData.hpp"

#include "mpiCommunications/MpiWrapper.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/Perforation.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"


namespace geosx
{

using namespace dataRepository;

PerforationData::PerforationData( string const & name, Group * const parent )
  : ObjectManagerBase( name, parent ),
  m_numPerforationsGlobal( 0 ),
  m_location( 0, 3 )
{
  registerWrapper( viewKeyStruct::numPerforationsGlobalString(), &m_numPerforationsGlobal );
  registerWrapper( viewKeyStruct::reservoirElementRegionString(), &m_toMeshElements.m_toElementRegion );
  registerWrapper( viewKeyStruct::reservoirElementSubregionString(), &m_toMeshElements.m_toElementSubRegion );
  registerWrapper( viewKeyStruct::reservoirElementIndexString(), &m_toMeshElements.m_toElementIndex );
  registerWrapper( viewKeyStruct::wellElementIndexString(), &m_wellElementIndex );
  registerWrapper( viewKeyStruct::locationString(), &m_location );
  registerWrapper( viewKeyStruct::wellTransmissibilityString(), &m_wellTransmissibility );
}

PerforationData::~PerforationData()
{}

namespace
{

/**
 * @brief Check if the well is along the x-, y-, or z- directions.
 * @tparam VEC_TYPE type of @p vecWellElemCenterToPerf
 * @tparam PERM_TYPE type of @p perm
 * @param[in] vecWellElemCenterToPerf vector connecting the well element center to the perforation
 * @param[in] dx dimension of the element in the x-direction
 * @param[in] dy dimension of the element in the y-direction
 * @param[in] dz dimension of the element in the z-direction
 * @param[in] perm absolute permeability in the reservoir element
 * @param[out] d1 dimension of the element in the first direction
 * @param[out] d2 dimension of the element in the second direction
 * @param[out] h dimension of the element in the third direction
 * @param[out] k1 absolute permeability in the reservoir element (first direction)
 * @param[out] k2 absolute permeability in the reservoir element (second direction)
 */
template< typename VEC_TYPE, typename PERM_TYPE >
void DecideWellDirection( VEC_TYPE const & vecWellElemCenterToPerf,
                          real64 const & dx,
                          real64 const & dy,
                          real64 const & dz,
                          PERM_TYPE const & perm,
                          real64 & d1,
                          real64 & d2,
                          real64 & h,
                          real64 & k1,
                          real64 & k2 )
{
  // vertical well (approximately along the z-direction)
  if( fabs( vecWellElemCenterToPerf[2] ) > fabs( vecWellElemCenterToPerf[0] )
      && fabs( vecWellElemCenterToPerf[2] ) > fabs( vecWellElemCenterToPerf[1] ) )
  {
    d1 = dx;
    d2 = dy;
    h  = dz;
    k1 = perm[0];
    k2 = perm[1];
  }
  // well approximately along the y-direction
  else if( fabs( vecWellElemCenterToPerf[1] ) > fabs( vecWellElemCenterToPerf[0] )
           && fabs( vecWellElemCenterToPerf[1] ) > fabs( vecWellElemCenterToPerf[2] ) )
  {
    d1 = dx;
    d2 = dz;
    h  = dy;
    k1 = perm[0];
    k2 = perm[2];
  }
  // well approximately along the x-direction
  else
  {
    d1 = dy;
    d2 = dz;
    h  = dx;
    k1 = perm[1];
    k2 = perm[2];
  }
}

}

void PerforationData::computeWellTransmissibility( MeshLevel const & mesh,
                                                   WellElementSubRegion const & wellElemSubRegion,
                                                   string const & permeabilityKey )
{

  // get the permeability in the domain
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const perm =
    mesh.getElemManager().constructArrayViewAccessor< real64, 2 >( permeabilityKey );

  arrayView2d< real64 const > const wellElemCenter = wellElemSubRegion.getElementCenter();

  // for all the local perforations on this well
  for( localIndex iperf = 0; iperf < size(); ++iperf )
  {

    // if the well transmissibility has been read from the XML
    // then skip the computation of the well transmissibility carried out below
    if( m_wellTransmissibility[iperf] > 0 )
    {
      continue;
    }

    // get the indices of the reservoir element
    localIndex const er  = m_toMeshElements.m_toElementRegion   [iperf];
    localIndex const esr = m_toMeshElements.m_toElementSubRegion[iperf];
    localIndex const ei  = m_toMeshElements.m_toElementIndex    [iperf];

    real64 dx = 0;
    real64 dy = 0;
    real64 dz = 0;

    // get an approximate dx, dy, dz for the reservoir element
    // this is done by computing a bounding box
    getReservoirElementDimensions( mesh, er, esr, ei, dx, dy, dz );

    real64 d1 = 0;
    real64 d2 = 0;
    real64 h  = 0;
    real64 k1 = 0;
    real64 k2 = 0;

    // compute the vector perforation - well elem center
    // this vector will be used to decide whether this is a vectical well or not
    localIndex const wellElemIndex   = m_wellElementIndex[iperf];
    real64 vecWellElemCenterToPerf[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( wellElemCenter[wellElemIndex] );
    LvArray::tensorOps::subtract< 3 >( vecWellElemCenterToPerf, m_location[iperf] );

    // check if this is a vertical well or a horizontal well
    // assign d1, d2, h, k1, and k2 accordingly
    DecideWellDirection( vecWellElemCenterToPerf,
                         dx, dy, dz, perm[er][esr][ei],
                         d1, d2, h, k1, k2 );

    real64 const k21 = k1 > 0
                     ? k2 / k1
                     : 0;
    real64 const k12 = k2 > 0
                     ? k1 / k2
                     : 0;

    // compute the equivalent radius
    real64 const num = 0.28 * sqrt( d1*d1 * sqrt( k21 ) + d2*d2 * sqrt( k12 ) );
    real64 const den = std::pow( k12, 0.25 ) + std::pow( k21, 0.25 );
    real64 const rEq = (den > 0)
                     ? num / den
                     : 0.0;

    real64 const kh = h * sqrt( k1 * k2 );

    arrayView1d< real64 const > const & wellElemRadius =
      wellElemSubRegion.getReference< array1d< real64 > >( WellElementSubRegion::viewKeyStruct::radiusString() );

    GEOSX_ERROR_IF( rEq < wellElemRadius[wellElemIndex],
                    "The equivalent radius r_eq = " << rEq <<
                    " is smaller than the well radius (r = " << wellElemRadius[wellElemIndex] <<
                    ") in " << getName() );

    // compute the well Peaceman index
    m_wellTransmissibility[iperf] = 2 * M_PI * kh / std::log( rEq / wellElemRadius[wellElemIndex] );

    GEOSX_ERROR_IF( m_wellTransmissibility[iperf] <= 0,
                    "The well index is negative or equal to zero in " << getName() );
  }
}


void PerforationData::getReservoirElementDimensions( MeshLevel const & mesh,
                                                     localIndex const er, localIndex const esr, localIndex const ei,
                                                     real64 & dx, real64 & dy, real64 & dz ) const
{
  ElementRegionManager const & elemManager = mesh.getElemManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  CellElementRegion const & region = elemManager.getRegion< CellElementRegion >( er );
  CellBlock const & subRegion = region.getSubRegion< CellElementSubRegion >( esr );

  // compute the bounding box of the element
  real64 boxDims[ 3 ];
  computationalGeometry::GetBoundingBox( ei,
                                         subRegion.nodeList(),
                                         nodeManager.referencePosition(),
                                         boxDims );

  // dx and dz from bounding box
  dx = boxDims[ 0 ];
  dy = boxDims[ 1 ];

  // dz is computed as vol / (dx * dy)
  dz  = subRegion.getElementVolume()[ei];
  dz /= dx * dy;

  GEOSX_ERROR_IF( dx <= 0 || dy <= 0 || dz <= 0,
                  "The reservoir element dimensions (dx, dy, and dz) should be positive in " << getName() );

}

void PerforationData::connectToWellElements( InternalWellGenerator const & wellGeometry,
                                             unordered_map< globalIndex, localIndex > const & globalToLocalWellElemMap,
                                             globalIndex elemOffsetGlobal )
{
  arrayView1d< globalIndex const > const & perfElemIndexGlobal = wellGeometry.getPerfElemIndex();

  for( localIndex iperfLocal = 0; iperfLocal < size(); ++iperfLocal )
  {
    globalIndex const iwelemGlobal = perfElemIndexGlobal[m_localToGlobalMap[iperfLocal]];
    globalIndex const ielemGlobal  = elemOffsetGlobal + iwelemGlobal;
    GEOSX_ASSERT( globalToLocalWellElemMap.count( ielemGlobal ) > 0 );
    m_wellElementIndex[iperfLocal] = globalToLocalWellElemMap.at( ielemGlobal );
  }

}

} //namespace geosx
