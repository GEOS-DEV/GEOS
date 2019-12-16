/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * @file PerforationData.cpp
 *
 */

#include "PerforationData.hpp"

#include "mpiCommunications/MpiWrapper.hpp"
#include "Perforation.hpp"

#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"

namespace geosx
{

using namespace dataRepository;

PerforationData::PerforationData(string const & name, Group * const parent)
  : ObjectManagerBase(name, parent),
    m_numPerforationsGlobal(0)
{
  registerWrapper( viewKeyStruct::numPerforationsGlobalString, &m_numPerforationsGlobal, false );
 
  registerWrapper( viewKeyStruct::reservoirElementRegionString,    &m_toMeshElements.m_toElementRegion,    false );
  registerWrapper( viewKeyStruct::reservoirElementSubregionString, &m_toMeshElements.m_toElementSubRegion, false );
  registerWrapper( viewKeyStruct::reservoirElementIndexString,     &m_toMeshElements.m_toElementIndex,     false );

  registerWrapper( viewKeyStruct::wellElementIndexString, &m_wellElementIndex, false );
  registerWrapper( viewKeyStruct::locationString, &m_location, false );
  registerWrapper( viewKeyStruct::transmissibilityString, &m_transmissibility, false );
}

PerforationData::~PerforationData()
{
}

void PerforationData::ConnectToMeshElements( MeshLevel const & mesh,
                                             InternalWellGenerator const & wellGeometry )
{
  ElementRegionManager const * const elemManager = mesh.getElemManager();
  NodeManager const * const nodeManager = mesh.getNodeManager();

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor const>> 
  elemCenter = elemManager->ConstructViewAccessor<array1d<R1Tensor>, arrayView1d<R1Tensor const>>( ElementSubRegionBase::
                                                                                                   viewKeyStruct::
                                                                                                   elementCenterString );
  
  arrayView1d<R1Tensor const> const & perfCoordsGlobal = wellGeometry.GetPerfCoords();
  arrayView1d<real64 const>   const & perfTransGlobal  = wellGeometry.GetPerfTransmissibility();

  resize( perfCoordsGlobal.size() );
  localIndex iperfLocal = 0;

  // loop over all the perforations
  for ( globalIndex iperfGlobal = 0; iperfGlobal < perfCoordsGlobal.size(); ++iperfGlobal )
  {
    R1Tensor const & coords = perfCoordsGlobal[iperfGlobal];

    // TODO actually trace coords
    // TODO what if a fracture element is located

    // find the closest reservoir element
    auto ret = minLocOverElemsInMesh( &mesh, [&] ( localIndex const er,
                                                   localIndex const esr,
                                                   localIndex const ei ) -> real64
    {
      R1Tensor v = coords;
      v -= elemCenter[er][esr][ei];
      return v.L2_Norm();
    } );

    // save the region, subregion and index
    localIndex const er  = std::get<0>(ret.second);
    localIndex const esr = std::get<1>(ret.second);
    localIndex const ei  = std::get<2>(ret.second);

    // a ghostRank check may be inserted if needed
    // (not needed in the current implementation)

    // check if perforation location is indeed inside the element
    CellBlock const * const cellBlock = elemManager->GetRegion( er )->GetSubRegion<CellElementSubRegion>( esr );
    array1d<array1d<localIndex>> faceNodes( cellBlock->numFacesPerElement() );

    for (localIndex kf = 0; kf < cellBlock->numFacesPerElement(); ++kf)
    {
      cellBlock->GetFaceNodes( ei, kf, faceNodes[kf] );
    }

    if (! computationalGeometry::IsPointInsidePolyhedron( nodeManager->referencePosition(), faceNodes, coords ))
    {
      continue;
    }

    // TODO: what happens when the boundary is at the boundary of the MPI domain??

    // now construct the local data

    // store the indices of the mesh element
    m_toMeshElements.m_toElementRegion   [iperfLocal] = er;
    m_toMeshElements.m_toElementSubRegion[iperfLocal] = esr;
    m_toMeshElements.m_toElementIndex    [iperfLocal] = ei;
    
    // construct the local transmissibility and location maps
    m_transmissibility[iperfLocal] = perfTransGlobal[iperfGlobal];
    m_location[iperfLocal] = coords;

    m_localToGlobalMap[iperfLocal++] = iperfGlobal;
  }

  // set the size based on the number of perforations matched with local reservoir elements
  resize( iperfLocal );
  ConstructGlobalToLocalMap();

}  
 
void PerforationData::ConnectToWellElements( InternalWellGenerator const & wellGeometry,
                                             unordered_map<globalIndex,localIndex> const & globalToLocalWellElemMap, 
                                             globalIndex elemOffsetGlobal )
{
  arrayView1d<globalIndex const> const & perfElemIndexGlobal = wellGeometry.GetPerfElemIndex();

  for (localIndex iperfLocal = 0; iperfLocal < size(); ++iperfLocal)
  {
    globalIndex const iwelemGlobal = perfElemIndexGlobal[m_localToGlobalMap[iperfLocal]];
    globalIndex const ielemGlobal  = elemOffsetGlobal + iwelemGlobal;
    GEOSX_ASSERT( globalToLocalWellElemMap.count( ielemGlobal ) > 0 );
    m_wellElementIndex[iperfLocal] = globalToLocalWellElemMap.at( ielemGlobal ); 
  }

  //DebugLocalPerforations();
}
 
void PerforationData::InitializePostInitialConditions_PreSubGroups( Group * const GEOSX_UNUSED_ARG( problemManager ) )
{
  for (localIndex iperf = 0; iperf < size(); ++iperf)
  {
    if (m_transmissibility[iperf] < 0.0)
    {
      // TODO: compute transmissibility internally
      GEOSX_ERROR( "Invalid transmissibility value: " << m_transmissibility[iperf] );
    }
  }
}


void PerforationData::DebugLocalPerforations() const
{
  if (size() == 0)
  {
    return;
  } 

  if ( MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) != 1)
  {
    return;
  } 

  std::cout << std::endl;
  std::cout << "++++++++++++++++++++++++++" << std::endl;
  std::cout << "PerforationData = " << getName() << " of " << getParent()->getName() << std::endl;
  std::cout << "MPI rank = " << MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) << std::endl;
  std::cout << "Number of local perforations = " << size() << std::endl;

  for (localIndex iperf = 0; iperf < size(); ++iperf) 
  {
    std::cout << "m_toMeshElements.m_toElementRegion[" << iperf << "] = "
              << m_toMeshElements.m_toElementRegion[iperf]
              << std::endl;
    std::cout << "m_toMeshElements.m_toElementSubRegion[" << iperf << "] = "
              << m_toMeshElements.m_toElementSubRegion[iperf]
              << std::endl;
    std::cout << "m_toMeshElements.m_toElementIndex[" << iperf << "] = "
              << m_toMeshElements.m_toElementIndex[iperf]
              << std::endl;

    std::cout << "m_wellElementIndexLocal[" << iperf << "] = " << m_wellElementIndex[iperf] 
              << std::endl;
    std::cout << "m_location[" << iperf << "] = " << m_location[iperf]
              << std::endl;
    std::cout << "m_transmissibility[" << iperf << "] = " << m_transmissibility[iperf]
              << std::endl;
  }
}

} //namespace geosx
