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

#include "WellElementSubRegion.hpp"
#include "WellElementRegion.hpp"


#include "mesh/MeshLevel.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mpiCommunications/MpiWrapper.hpp"

namespace geosx
{

WellElementSubRegion::WellElementSubRegion( string const & name, Group * const parent ):
  ElementSubRegionBase( name, parent ),
  m_wellControlsName(""),
  m_toNodesRelation(),
  m_topWellElementIndex( -1 ),
  m_perforationData( groupKeyStruct::perforationDataString, this ),
  m_topRank( -1 )
{

  registerWrapper( viewKeyStruct::wellControlsString, &m_wellControlsName, false );
  registerWrapper( viewKeyStruct::wellNodeListString, &m_toNodesRelation, false );
  registerWrapper( viewKeyStruct::nextWellElementIndexString, &m_nextWellElementIndex, false );
  registerWrapper( viewKeyStruct::nextWellElementIndexGlobalString, &m_nextWellElementIndexGlobal, false );
  registerWrapper( viewKeyStruct::topWellElementIndexString, &m_topWellElementIndex, false );
  registerWrapper( viewKeyStruct::topRankString, &m_topRank, false );
  registerWrapper( viewKeyStruct::radiusString, &m_radius, false );

  registerWrapper( ElementSubRegionBase::viewKeyStruct::elementCenterString, &m_elementCenter, false );
  registerWrapper( ElementSubRegionBase::viewKeyStruct::elementVolumeString, &m_elementVolume, false );

  RegisterGroup( groupKeyStruct::perforationDataString, &m_perforationData, false );

  this->numNodesPerElement() = 2;
  this->numFacesPerElement() = 0;
  m_toNodesRelation.resizeDimension<1>( this->numNodesPerElement() );
  m_elementTypeString = "BEAM";
}


WellElementSubRegion::~WellElementSubRegion()
{}


void WellElementSubRegion::setupRelatedObjectsInRelations( MeshLevel const * const mesh )
{
  m_toNodesRelation.SetRelatedObject( mesh->getNodeManager() );
}

void WellElementSubRegion::Generate( MeshLevel                        & mesh, 
                                     InternalWellGenerator      const & wellGeometry,
                                     arrayView1d<integer>             & elemStatusGlobal, 
                                     globalIndex                        nodeOffsetGlobal,
                                     globalIndex                        elemOffsetGlobal )
{

  map< integer, set<globalIndex> > elemSetsByStatus;
 
  // convert elemStatus list into sets of indices
  for (localIndex iwelemGlobal = 0; iwelemGlobal < elemStatusGlobal.size(); ++iwelemGlobal)
  {
    elemSetsByStatus[elemStatusGlobal[iwelemGlobal]].insert( iwelemGlobal );
  }

  // initialize the sets using the classification of well elems
  // localElems will be enlarged once boundary elements ownership is determined
  set<globalIndex> & localElems   = elemSetsByStatus[WellElemStatus::LOCAL];
  set<globalIndex> & sharedElems  = elemSetsByStatus[WellElemStatus::SHARED];
  set<globalIndex> & unownedElems = elemSetsByStatus[WellElemStatus::UNOWNED];

  // here we make sure that there are no shared elements
  // this is enforced in the InternalWellGenerator that currently merges two perforations 
  // if they belong to the same well element. This is a temporary solution.
  // TODO: split the well elements that contain multiple perforations, so that no element is shared
  GEOSX_ERROR_IF( sharedElems.size() > 0,
                 "Well " << getName() << " contains shared well elements");


  // In Steps 1 and 2 we determine the local objects on this rank (elems and nodes)
  // Once this is done, in Steps 3, 4, and 5, we update the nodeManager and wellElementSubRegion (size, maps)

  
  // 1) First assign the unowned elements to a rank
  // this is done in two steps

  // 1.a) First assign unowned elements in the reservoir based on location
  //      ie., if the center of the well element falls in the domain owned by rank k
  //      then the well element is assigned to rank k
  AssignUnownedElementsInReservoir( mesh, 
                                    wellGeometry,
                                    unownedElems,
                                    localElems,
                                    elemStatusGlobal );
  // 1.b) Then we check that all the well elements have been assigned (and assigned once)
  //      This is needed because if the center of the well element falls on the boundary of 
  //      a reservoir element, the assignment algorithm of 1.a) can assign the same well element
  //      to two ranks, or to no rank at all (which will break the solver).
  //      In this function we also check that the resulting well partitioning is valid, that is,
  //      we make sure that if two ranks are neighbors in the well, that are also neighbors in the
  //      reservoir mesh
  CheckPartitioningValidity( wellGeometry, 
                             localElems,
                             elemStatusGlobal );

  set<globalIndex> localNodes;
  set<globalIndex> boundaryNodes;

  // 2) collect the local nodes and tag the boundary nodes using element info
  // now that all the elements have been assigned, we collected the local nodes
  // and tag the boundary nodes (i.e., the nodes in contact with both local and remote elems)
  CollectLocalAndBoundaryNodes( wellGeometry, 
                                localElems,
                                localNodes,
                                boundaryNodes );  

  //DebugNodeManager( mesh );

  // 3) size update in the nodeManager
  // this is necessary to later use the node matching procedure
  // to place ghosts in DomainPartition::SetupCommunications
  UpdateNodeManagerSize( mesh,
                         wellGeometry,
                         localNodes,
                         boundaryNodes,
                         nodeOffsetGlobal );

  //DebugNodeManager( mesh );

  // 4) resize the well element subregion 
  // and construct local to global, global to local, maps, etc
  ConstructSubRegionLocalElementMaps( mesh,
                                      wellGeometry,
                                      localElems,
                                      nodeOffsetGlobal,
                                      elemOffsetGlobal );

  //DebugWellElementSubRegions( elemStatusGlobal, elemOffsetGlobal );

  // 5) node-to-elem map update in the nodeManager
  // This map will be used by MeshLevel::GenerateAdjacencyLists
  // this assumes that the elemToNodes maps has been filled at Step 5)
  UpdateNodeManagerNodeToElementMap( mesh );

}


void WellElementSubRegion::AssignUnownedElementsInReservoir( MeshLevel                   & mesh,
                                                             InternalWellGenerator const & wellGeometry,
                                                             set<globalIndex>      const & unownedElems,
                                                             set<globalIndex>            & localElems,
                                                             arrayView1d<integer>        & elemStatusGlobal ) const
{
  NodeManager const * const nodeManager          = mesh.getNodeManager();
  ElementRegionManager const * const elemManager = mesh.getElemManager();

  // get the well and reservoir element coordinates
  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor const>> 
  resElemCoords = elemManager->ConstructViewAccessor<array1d<R1Tensor>, arrayView1d<R1Tensor const>>( ElementSubRegionBase::
                                                                                                      viewKeyStruct::
                                                                                                      elementCenterString );
  arrayView1d<R1Tensor const> const & wellElemCoordsGlobal = wellGeometry.GetElemCoords();

  // assign the well elements based on location wrt the reservoir elements
  // if the center of the well element falls in the domain owned by rank k
  // then the well element is assigned to rank k
  for (globalIndex currGlobal : unownedElems)
  {
    R1Tensor const & wellElemCoords = wellElemCoordsGlobal[currGlobal];



    // find the closest reservoir element
    auto ret = minLocOverElemsInMesh( &mesh, [&] ( localIndex const er,
                                                   localIndex const esr,
                                                   localIndex const ei ) -> real64
    {
      R1Tensor v = wellElemCoords;
      v -= resElemCoords[er][esr][ei];
      return v.L2_Norm();
    } );

    // save the region, subregion and index
    localIndex const er  = std::get<0>(ret.second);
    localIndex const esr = std::get<1>(ret.second);
    localIndex const ei  = std::get<2>(ret.second);

    // check if well element center location is indeed inside the element
    CellBlock const * const cellBlock = elemManager->GetRegion( er )->GetSubRegion<CellElementSubRegion>( esr );
    array1d<array1d<localIndex>> faceNodes( cellBlock->numFacesPerElement() );

    for (localIndex kf = 0; kf < cellBlock->numFacesPerElement(); ++kf)
    {
      cellBlock->GetFaceNodes( ei, kf, faceNodes[kf] );
    }

    if (! computationalGeometry::IsPointInsidePolyhedron( nodeManager->referencePosition(), faceNodes, wellElemCoords ))
    {
      continue;
    }

    // the well element is in the reservoir element (er,esr,ei), so tag it as local
    localElems.insert( currGlobal );
    elemStatusGlobal[currGlobal] = WellElemStatus::LOCAL;
  }
}


void WellElementSubRegion::CheckPartitioningValidity( InternalWellGenerator const & wellGeometry,
                                                      set<globalIndex>            & localElems,
                                                      arrayView1d<integer>        & elemStatusGlobal ) const
{
  arrayView1d< arrayView1d< globalIndex const > const > const & prevElemIdsGlobal = wellGeometry.GetPrevElemIndices();

  // we are going to make sure that the partitioning is good, 
  // well element per well element, starting from the bottom of the well
  for (globalIndex iwelemGlobal = wellGeometry.GetNumElements()-1; iwelemGlobal >= 0; --iwelemGlobal)
  {

    // communicate the status of this element
    array1d<integer> thisElemStatusGlobal;
    MpiWrapper::allGather( elemStatusGlobal[iwelemGlobal],
                                   thisElemStatusGlobal );
    // group the ranks by well element status
    map< integer, set<globalIndex> > rankSetsByStatus;
    for (globalIndex irank = 0; irank < thisElemStatusGlobal.size(); ++irank)
    {
      rankSetsByStatus[thisElemStatusGlobal[irank]].insert( irank );
    }
    globalIndex const numLocalRanks = rankSetsByStatus[WellElemStatus::LOCAL].size();

    // in this case, this element has not been assigned
    //    => we assign it to the rank that owns 
    //       the well element below iwelemGlobal (prevGlobal, already assigned and checked)
    if (numLocalRanks == 0)
    { 
      globalIndex const numBranches = prevElemIdsGlobal[iwelemGlobal].size();
      globalIndex const prevGlobal  = prevElemIdsGlobal[iwelemGlobal][numBranches-1];

      GEOSX_ERROR_IF( prevGlobal <= iwelemGlobal || prevGlobal < 0,
                    "Invalid partitioning in well " << getName() );

      if (elemStatusGlobal[prevGlobal] == WellElemStatus::LOCAL)
      { 
        localElems.insert(iwelemGlobal);
        elemStatusGlobal[iwelemGlobal] = WellElemStatus::LOCAL;
      }
      else
      {
        elemStatusGlobal[iwelemGlobal] = WellElemStatus::REMOTE;
      }
    }  

    // in this case, everything is fine, 
    // we just update the elemStatusGlobal array for all ranks
    else if (numLocalRanks == 1)
    {

      for (globalIndex iownerRank : rankSetsByStatus[WellElemStatus::LOCAL])
      {
        if (MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) != iownerRank)
        {
          elemStatusGlobal[iwelemGlobal] = WellElemStatus::REMOTE;
        }
      }

    }     

    // in this case, this element has been assigned to more than rank
    //    => the smallest rank keeps it
    else // (numLocalRanks > 1)
    {

      localIndex rankCount = 0;
      for (globalIndex iownerRank : rankSetsByStatus[WellElemStatus::LOCAL])
      {
        if (rankCount == 0)
        {
          // update the elemStatusGlobal array for all ranks
          if (MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) != iownerRank)
          {
            elemStatusGlobal[iwelemGlobal] = WellElemStatus::REMOTE;
          }
        }
        else // (rankCount > 0)
        {
          // remove the duplicate elements
          if (MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) == iownerRank)
          {
            localElems.erase(iwelemGlobal);
          }
        } 
        rankCount++;
      }

    }

    // TODO: check neighbor rank
  }
}


void WellElementSubRegion::CollectLocalAndBoundaryNodes( InternalWellGenerator const & wellGeometry, 
                                                         set<globalIndex>      const & localElems,
                                                         set<globalIndex>            & localNodes,
                                                         set<globalIndex>            & boundaryNodes ) const
{
  // get the well connectivity
  arrayView1d< globalIndex const >                      const & nextElemIdGlobal  = wellGeometry.GetNextElemIndex();
  arrayView1d< arrayView1d< globalIndex const > const > const & prevElemIdsGlobal = wellGeometry.GetPrevElemIndices();
  arrayView2d< globalIndex const >                      const & elemToNodesGlobal = wellGeometry.GetElemToNodesMap();

  // loop over the local elements and collect the local and boundary nodes
  for (globalIndex currGlobal : localElems)
  {

    // if the element is local, its two nodes are also local
    globalIndex const inodeTopGlobal    = elemToNodesGlobal[currGlobal][InternalWellGenerator::NodeLocation::TOP];
    globalIndex const inodeBottomGlobal = elemToNodesGlobal[currGlobal][InternalWellGenerator::NodeLocation::BOTTOM];
    localNodes.insert( inodeTopGlobal );
    localNodes.insert( inodeBottomGlobal );

    localIndex const nextGlobal = 
      integer_conversion<localIndex>( nextElemIdGlobal[ integer_conversion<localIndex>(currGlobal) ] );

    // if the next well elem is not local, add the node in between curr and next to boundaryNodes
    if (nextGlobal >= 0 && !localElems.contains( nextGlobal ))
    {
      boundaryNodes.insert( inodeTopGlobal );
    } 

    // if the prev well elem is not local, add the node in between curr and prev to boundaryNodes (relevant for branches)
    for (localIndex iwelem = 0; iwelem < prevElemIdsGlobal[currGlobal].size(); ++iwelem)
    {
      globalIndex const prevGlobal = prevElemIdsGlobal[currGlobal][iwelem];
      if (prevGlobal >= 0 && !localElems.contains( prevGlobal ))
      {
        boundaryNodes.insert( inodeBottomGlobal );
      } 
    }
  }
}

void WellElementSubRegion::UpdateNodeManagerSize( MeshLevel                   & mesh,
                                                  InternalWellGenerator const & wellGeometry,
                                                  set<globalIndex>      const & localNodes,
                                                  set<globalIndex>      const & boundaryNodes, 
                                                  globalIndex                   nodeOffsetGlobal )
{

 // get the node manager to compute the total number of mesh nodes
  NodeManager * const nodeManager       = mesh.getNodeManager();
  localIndex    const numWellNodesLocal = localNodes.size();
  localIndex    const oldNumNodesLocal  = nodeManager->size();

  // resize nodeManager to account for the new well nodes and update the properties
  nodeManager->resize( oldNumNodesLocal + numWellNodesLocal );

  array1d<integer> & 
  isDomainBoundary = nodeManager->getReference<integer_array>(m_ObjectManagerBaseViewKeys.domainBoundaryIndicator);
  arrayView1d<R1Tensor const> const & nodeCoordsGlobal = wellGeometry.GetNodeCoords();

  // local *well* index
  localIndex iwellNodeLocal = 0; 
  // loop over global *well* indices
  for ( globalIndex iwellNodeGlobal : localNodes ) 
  {
    // local *nodeManager* index
    localIndex const inodeLocal = oldNumNodesLocal + iwellNodeLocal; 
 
    // update node manager maps and position 
    nodeManager->m_localToGlobalMap[inodeLocal]  = nodeOffsetGlobal + iwellNodeGlobal; // global *nodeManager* index
    nodeManager->referencePosition()[inodeLocal] = nodeCoordsGlobal[iwellNodeGlobal];

    // mark the boundary nodes for ghosting in DomainPartition::SetupCommunications
    if ( boundaryNodes.contains( iwellNodeGlobal ) )
    {
      isDomainBoundary[inodeLocal] = 1;
    } 

    iwellNodeLocal++;
  }

  // now with update the relevant node indices in nodeManager->globalToLocalMap
  // this is to avoid a call to nodeManager->ConstructGlobalToLocalMap everytime we add a well
  for (iwellNodeLocal = 0; iwellNodeLocal < numWellNodesLocal; ++iwellNodeLocal)
  {
    // local *nodeManager* index
    localIndex const inodeLocal = oldNumNodesLocal + iwellNodeLocal; 

    nodeManager->m_globalToLocalMap[nodeManager->m_localToGlobalMap[inodeLocal]] = inodeLocal;
  }
}

void WellElementSubRegion::ConstructSubRegionLocalElementMaps( MeshLevel                   & mesh,
                                                               InternalWellGenerator const & wellGeometry,
                                                               set<globalIndex>      const & localElems, 
                                                               globalIndex                   nodeOffsetGlobal,
                                                               globalIndex                   elemOffsetGlobal )
{
  // get the well geometry
  arrayView1d<globalIndex const> const & nextElemIdGlobal  = wellGeometry.GetNextElemIndex();
  arrayView1d<R1Tensor const>    const & elemCoordsGlobal  = wellGeometry.GetElemCoords();
  arrayView2d<globalIndex const> const & elemToNodesGlobal = wellGeometry.GetElemToNodesMap();
  arrayView1d<real64 const>      const & elemVolumeGlobal  = wellGeometry.GetElemVolume();

  NodeManager const * const nodeManager = mesh.getNodeManager();

  resize( localElems.size() );

  // create local elem numbering

  // local well elem ordering
  localIndex iwelemLocal = 0;
  // loop over global well elem indices
  for ( globalIndex iwelemGlobal : localElems )
  {
    // create a global *elemManager* index
    m_localToGlobalMap[iwelemLocal++] = elemOffsetGlobal + iwelemGlobal; 
  }  
  ConstructGlobalToLocalMap();
  
  // recreate local wellbore tree by connecting locally relevant elems
  for (iwelemLocal = 0; iwelemLocal < size(); ++iwelemLocal)
  {
    globalIndex const ielemGlobal      = m_localToGlobalMap[iwelemLocal];     // global index in elemManager ordering
    globalIndex const iwelemGlobal     = ielemGlobal - elemOffsetGlobal;      // global index in well ordering
    globalIndex const iwelemNextGlobal = nextElemIdGlobal[iwelemGlobal];      // global index in well ordering
    globalIndex const ielemNextGlobal  = elemOffsetGlobal + iwelemNextGlobal; // global index in elemManager ordering

    if (iwelemNextGlobal < 0)
    {
      m_nextWellElementIndexGlobal[iwelemLocal] = -1; // wellhead
      m_nextWellElementIndex[iwelemLocal]       = -1;
      m_topWellElementIndex = iwelemLocal;
    }
    else 
    {
      m_nextWellElementIndexGlobal[iwelemLocal] = ielemNextGlobal; // wellhead
      
      if (m_globalToLocalMap.count( ielemNextGlobal ) > 0)
      {
        m_nextWellElementIndex[iwelemLocal] = m_globalToLocalMap.at( ielemNextGlobal );
      }
      else
      {
        m_nextWellElementIndex[iwelemLocal] = -2; // remote elem
      }
    }
    m_elementCenter[iwelemLocal] = elemCoordsGlobal[iwelemGlobal];
    m_elementVolume[iwelemLocal] = elemVolumeGlobal[iwelemGlobal];
    m_radius[iwelemLocal] = wellGeometry.GetElementRadius();

    // update local well elem to node map (note: nodes are in nodeManager ordering)

    // first get the global node indices in nodeManager ordering
    globalIndex const inodeTopGlobal    = nodeOffsetGlobal + elemToNodesGlobal[iwelemGlobal][InternalWellGenerator::NodeLocation::TOP];
    globalIndex const inodeBottomGlobal = nodeOffsetGlobal + elemToNodesGlobal[iwelemGlobal][InternalWellGenerator::NodeLocation::BOTTOM];

    // then get the local node indices in nodeManager ordering
    localIndex const inodeTopLocal    = nodeManager->m_globalToLocalMap.at( inodeTopGlobal );
    localIndex const inodeBottomLocal = nodeManager->m_globalToLocalMap.at( inodeBottomGlobal );

    m_toNodesRelation[iwelemLocal][InternalWellGenerator::NodeLocation::TOP]    = inodeTopLocal;
    m_toNodesRelation[iwelemLocal][InternalWellGenerator::NodeLocation::BOTTOM] = inodeBottomLocal;
  }

}

void WellElementSubRegion::UpdateNodeManagerNodeToElementMap( MeshLevel & mesh ) 
{
  ElementRegionManager const * const elemManager = mesh.getElemManager();
  NodeManager * const nodeManager = mesh.getNodeManager();

  // at this point, NodeManager::SetElementMaps has already been called for the mesh nodes
  // we have to update the following maps for the well nodes
  ArrayOfArrays<localIndex> & toElementRegionList    = nodeManager->elementRegionList();
  ArrayOfArrays<localIndex> & toElementSubRegionList = nodeManager->elementSubRegionList();
  ArrayOfArrays<localIndex> & toElementList          = nodeManager->elementList();

  // we get the region and subregion indices in the elemManager
  WellElementRegion const * const elemRegion = this->getParent()->getParent()->group_cast<WellElementRegion*>();
  string const elemRegionName = elemRegion->getName();
      
  localIndex const iregion    = elemManager->GetRegions().getIndex( elemRegionName );
  localIndex const isubRegion = elemRegion->GetSubRegions().getIndex( getName() );

  // for each (new) well element
  for(localIndex iwelemLocal = 0;  iwelemLocal < size() ; ++iwelemLocal)
  {
     for(localIndex a=0; a < numNodesPerElement(); ++a)
     { 
       // get the local node index (in nodeManager ordering) using the elem-to-nodes maps constructed above
       localIndex const inodeLocal = m_toNodesRelation[iwelemLocal][a];

       // update the reverse map from well node to well element
       // this is needed to generate the adjacency list in communication setup phase
       toElementRegionList.appendToArray( inodeLocal, iregion );
       toElementSubRegionList.appendToArray( inodeLocal, isubRegion );
       toElementList.appendToArray( inodeLocal, iwelemLocal );
     }
  }
  setupRelatedObjectsInRelations( &mesh );
}


void WellElementSubRegion::ReconstructLocalConnectivity()
{
  // here we reconstruct the array m_nextWellElementIndexGlobal
  // this is needed after the addition of ghost well elements

  for (localIndex iwelemLocal = 0; iwelemLocal < size(); ++iwelemLocal)
  {
    globalIndex const nextGlobal = m_nextWellElementIndexGlobal[iwelemLocal];

    if ( nextGlobal < 0 ) // well head
    {
      m_nextWellElementIndex[iwelemLocal] = -1;
      m_topWellElementIndex = iwelemLocal; // reset this is case top element was added as ghost
    }
    else if ( m_globalToLocalMap.count( nextGlobal ) == 0 ) // next is remote
    {
      m_nextWellElementIndex[iwelemLocal] = -2;
    }
    else // local
    {
      m_nextWellElementIndex[iwelemLocal] = this->m_globalToLocalMap[nextGlobal];
    }
  }
}


bool WellElementSubRegion::IsLocallyOwned() const
{
  return m_topRank == MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
}

void WellElementSubRegion::ViewPackingExclusionList( set<localIndex> & exclusionList ) const
{ 
  ObjectManagerBase::ViewPackingExclusionList(exclusionList);
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::nodeListString));
}

localIndex WellElementSubRegion::PackUpDownMapsSize( arrayView1d<localIndex const> const & packList ) const
{ 
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate<false>( junk, packList );
}

localIndex WellElementSubRegion::PackUpDownMaps( buffer_unit_type * & buffer,
                                                 arrayView1d<localIndex const> const & packList ) const
{ 
  return PackUpDownMapsPrivate<true>( buffer, packList );
}

template< bool DOPACK >
localIndex WellElementSubRegion::PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                                        arrayView1d<localIndex const> const & packList ) const
{
  localIndex packedSize = 0;

  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         nodeList().Base().toViewConst(),
                                         m_unmappedGlobalIndicesInNodelist,
                                         packList,
                                         this->m_localToGlobalMap,
                                         nodeList().RelatedObjectLocalToGlobal() );

  return packedSize;
}

localIndex WellElementSubRegion::UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                                   localIndex_array & packList,
                                                   bool const GEOSX_UNUSED_ARG( overwriteUpMaps ),
                                                   bool const GEOSX_UNUSED_ARG( overwriteDownMaps ) )
{ 
  localIndex unPackedSize = 0;

  unPackedSize += bufferOps::Unpack( buffer,
                                     nodeList().Base().toView(),
                                     packList,
                                     m_unmappedGlobalIndicesInNodelist,
                                     this->m_globalToLocalMap,
                                     nodeList().RelatedObjectGlobalToLocal() );

  return unPackedSize;
}

void WellElementSubRegion::FixUpDownMaps( bool const clearIfUnmapped )
{ 
  ObjectManagerBase::FixUpDownMaps( nodeList(),
                                    m_unmappedGlobalIndicesInNodelist,
                                    clearIfUnmapped );
}

void WellElementSubRegion::DebugNodeManager( MeshLevel const & mesh ) const
{
  NodeManager const * const nodeManager = mesh.getNodeManager();

  if ( MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) != 1)
  {
    return;
  } 
  
  std::cout << std::endl;
  std::cout << "++++++++++++++++++++++++++" << std::endl;
  std::cout << "Node manager from = " << getName() << std::endl;
  std::cout << "MPI rank = " << MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) << std::endl;
  std::cout << "Number of local node elements = " << nodeManager->size() << std::endl;

  if (nodeManager->size() > 0)
  {
    return;
  }

  for (localIndex inodeLocal = 0; inodeLocal < nodeManager->size(); ++inodeLocal)
  {
    std::cout << "nodeManager->localToGlobalMap["    << inodeLocal << "] = " << nodeManager->m_localToGlobalMap[inodeLocal]  << std::endl;
    std::cout << "nodeManager->referencePosition()[" << inodeLocal << "] = " << nodeManager->referencePosition()[inodeLocal] << std::endl;
  }
}
 
void WellElementSubRegion::DebugWellElementSubRegions( arrayView1d<integer const> const & elemStatusGlobal, globalIndex elemOffsetGlobal ) const 
{
  if (size() == 0)
  {
    return;
  } 

  if ( MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) < 1)
  {
    return;
  } 

  std::cout << std::endl;
  std::cout << "++++++++++++++++++++++++++" << std::endl;
  std::cout << "WellElementSubRegion = " << getName() << std::endl;
  std::cout << "MPI rank = " << MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) << std::endl;
  std::cout << "Number of local well elements = " << size() << std::endl;
  
  for (localIndex iwelem = 0; iwelem < size(); ++iwelem) 
  {

    std::cout << "m_nextWellElementIndex[" << iwelem << "] = " 
              << m_nextWellElementIndex[iwelem] 
              << std::endl;

    std::cout << "m_nextWellElementIndexGlobal[" << iwelem << "] = " 
              << m_nextWellElementIndexGlobal[iwelem] 
              << std::endl;

    std::cout << "m_elementCenter[" << iwelem << "] = " 
              << m_elementCenter[iwelem]
              << std::endl;

    std::cout << "m_localToGlobalMap[" << iwelem << "] = " 
              << m_localToGlobalMap[iwelem]
              << std::endl;

    std::cout << "elemStatusGlobal[" << m_localToGlobalMap[iwelem] << "] = " 
              << elemStatusGlobal[m_localToGlobalMap[iwelem] - elemOffsetGlobal] 
              << std::endl;
     
    std::cout << "m_topElementIndex = " << m_topWellElementIndex << std::endl;          
 
  }
}

void WellElementSubRegion::DebugWellElementSubRegionsAfterSetupCommunications() const 
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
  std::cout << "WellElementSubRegion = " << getName() << std::endl;
  std::cout << "MPI rank = " << MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) << std::endl;
  std::cout << "Number of local well elements = " << size() << std::endl;
  std::cout << "Number of ghost well elements = " << this->GetNumberOfGhosts() << std::endl;
  
  for (localIndex iwelem = 0; iwelem < size(); ++iwelem) 
  {

    std::cout << "m_nextWellElementIndex[" << iwelem << "] = " 
              << m_nextWellElementIndex[iwelem] 
              << std::endl;

    std::cout << "m_nextWellElementIndexGlobal[" << iwelem << "] = " 
              << m_nextWellElementIndexGlobal[iwelem] 
              << std::endl;

    std::cout << "m_elementCenter[" << iwelem << "] = " 
              << m_elementCenter[iwelem]
              << std::endl;

    std::cout << "m_localToGlobalMap[" << iwelem << "] = " 
              << m_localToGlobalMap[iwelem]
              << std::endl;

    std::cout << "m_ghostRank[" << iwelem << "] = " 
              << m_ghostRank[iwelem] 
              << std::endl;
     
    std::cout << "m_topElementIndex = " << m_topWellElementIndex << std::endl;          
 
  }
}

}
