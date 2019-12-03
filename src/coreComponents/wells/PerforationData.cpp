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
  registerWrapper( viewKeyStruct::wellPeacemanIndexString, &m_wellPeacemanIndex, false );
}

PerforationData::~PerforationData()
{
}

void PerforationData::ConnectToMeshElements( MeshLevel const & mesh,
                                             InternalWellGenerator const & wellGeometry )
{
  arrayView1d<R1Tensor const> const & perfCoordsGlobal    = wellGeometry.GetPerfCoords();
  arrayView1d<real64 const>   const & perfWellIndexGlobal = wellGeometry.GetPerfPeacemanIndex();

  resize( perfCoordsGlobal.size() );
  localIndex iperfLocal = 0;

  // loop over all the perforations
  for ( globalIndex iperfGlobal = 0; iperfGlobal < perfCoordsGlobal.size(); ++iperfGlobal )
  {
    R1Tensor const & location = perfCoordsGlobal[iperfGlobal];

    localIndex erMatched  = -1;
    localIndex esrMatched = -1;
    localIndex eiMatched  = -1;

    localIndex erInit  = -1;
    localIndex esrInit = -1;
    localIndex eiInit  = -1;
    
    if (iperfLocal > 0)
    {
      // get the info of the element matched with the previous perforation  
      erInit  = m_toMeshElements.m_toElementRegion[iperfLocal-1];
      esrInit = m_toMeshElements.m_toElementSubRegion[iperfLocal-1];
      eiInit  = m_toMeshElements.m_toElementIndex[iperfLocal-1];
    }	
    else
    {
      // start from the element that is the closest from the perforation  
      InitializeLocalSearch( mesh, location,
			     erInit, esrInit, eiInit );
    }	
	
    // Strategy #1: search locally, starting from the location of the previous perforation
    // the assumption here is that perforations have been entered in order of depth
    bool resElemFound = false;

    CellElementRegion const * region = Group::group_cast<CellElementRegion const *>(mesh.getElemManager()->GetRegion(erInit));
    CellBlock const * subRegion      = Group::group_cast<CellBlock const *>(region->GetSubRegion(esrInit));
           
    set<localIndex>  nodes;
    set<globalIndex> elements;

    // collect the nodes of the current element
    // they will be used to access the neighbors and check if they contain the perforation
    CollectNodes( subRegion, eiInit, nodes); 

    // enlarge the neighborhoods four times (if no match is found)
    // before switching to the other strategy
    localIndex const depth = 4;
    for (localIndex d = 0; d < depth; ++d)
    {
      // search the elements that can be access from the set "nodes"
      // stop if an element containing the perforation is found
      // if not, enlarge the set "nodes"
      resElemFound = SearchLocalElements( mesh, location, nodes, elements,
                                          erMatched, esrMatched, eiMatched );
      if (resElemFound)
      {
        break;
      }
    }  

    // Strategy #2: brute force 
    // if the local search was not sucessful, let's search in the entire domain
    if (resElemFound == false)  
    {
      resElemFound = SearchEntireDomain( mesh, location,
                                         erMatched, esrMatched, eiMatched );
    }
    
    // if the element was found
    if (resElemFound)
    {
      ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor const>> 
      elemCenter = mesh.getElemManager()->ConstructViewAccessor<array1d<R1Tensor>, arrayView1d<R1Tensor const>>( ElementSubRegionBase::
                                                                                                                 viewKeyStruct::
                                                                                                                 elementCenterString );
      std::cout << "CellElementCenter = "   << elemCenter[erMatched][esrMatched][eiMatched] << std::endl;
      std::cout << "PerforationLocation = " << location                << std::endl;
      
      // set the indices for the matched element
      m_toMeshElements.m_toElementRegion   [iperfLocal] = erMatched;
      m_toMeshElements.m_toElementSubRegion[iperfLocal] = esrMatched;
      m_toMeshElements.m_toElementIndex    [iperfLocal] = eiMatched;

      // construct the local wellPeacemanIndex and location maps
      m_wellPeacemanIndex[iperfLocal] = perfWellIndexGlobal[iperfGlobal];
      m_location[iperfLocal] = location;

      // increment the local to global map
      m_localToGlobalMap[iperfLocal++] = iperfGlobal;
    }
  }

  // set the size based on the number of perforations matched with local reservoir elements
  resize( iperfLocal );
  ConstructGlobalToLocalMap();

}

void PerforationData::InitializeLocalSearch( MeshLevel const  & mesh,
                                             R1Tensor  const  & location,
                                             localIndex       & erInit,
                                             localIndex       & esrInit,
                                             localIndex       & eiInit) const
{
  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor const>> 
  resElemCenter = mesh.getElemManager()->ConstructViewAccessor<array1d<R1Tensor>, arrayView1d<R1Tensor const>>( ElementSubRegionBase::
                                                                                                                viewKeyStruct::
                                                                                                                elementCenterString );

  // find the closest reservoir element
  auto ret = minLocOverElemsInMesh( &mesh, [&] ( localIndex const er,
                                                 localIndex const esr,
                                                 localIndex const ei ) -> real64
  {
    R1Tensor v = location;
    v -= resElemCenter[er][esr][ei];
    return v.L2_Norm();
  });

  // save the region, subregion and index
  erInit  = std::get<0>(ret.second);
  esrInit = std::get<1>(ret.second);
  eiInit  = std::get<2>(ret.second);
}

bool PerforationData::SearchLocalElements( MeshLevel const  & mesh,
                                           R1Tensor  const  & location,
                                           set<localIndex>  & nodes,
                                           set<globalIndex> & elements,
                                           localIndex       & erMatched,
                                           localIndex       & esrMatched,
                                           localIndex       & eiMatched ) const
{
  ElementRegionManager const * const elemManager = mesh.getElemManager();
  NodeManager const * const nodeManager          = mesh.getNodeManager();

  ArrayOfArraysView<localIndex const> const & toElementRegionList    = nodeManager->elementRegionList();
  ArrayOfArraysView<localIndex const> const & toElementSubRegionList = nodeManager->elementSubRegionList();
  ArrayOfArraysView<localIndex const> const & toElementList          = nodeManager->elementList();

  bool matched = false;

  // we will enlarge the set of nodes in the loop below
  // to do this we have to create a new set that will
  // contain the already visited nodes only
  set<localIndex> currNodes = nodes;
  
  // for all the nodes already visited
  for (localIndex currNode : currNodes)
  {
    // collect the elements that have not been visited yet
    for( localIndex b=0 ; b<toElementRegionList.sizeOfArray(currNode) ; ++b )
    {
      localIndex  const er      = toElementRegionList[currNode][b];
      localIndex  const esr     = toElementSubRegionList[currNode][b];
      localIndex  const eiLocal = toElementList[currNode][b];

      CellElementRegion const * const region    = Group::group_cast<CellElementRegion const *>(elemManager->GetRegion(er));
      CellBlock         const * const subRegion = Group::group_cast<CellElementSubRegion const *>(region->GetSubRegion(esr));
      globalIndex               const eiGlobal  = subRegion->m_localToGlobalMap[eiLocal];

      // if this element has not been visited yet, save it
      if (!elements.contains(eiGlobal))
      {
        elements.insert(eiGlobal);
        
        // perform the test to see if the point is in this element
        // if the point is in the element, save the indices and stop the search
        if (IsPointInsideElement( nodeManager, location, subRegion, eiLocal ))
        {  
          erMatched  = er;
          esrMatched = esr;
          eiMatched  = eiLocal;
          matched    = true;
          break;
        }
        // otherwise add the nodes of this element to the set of new nodes to visit
        else
        {
          CollectNodes( subRegion, eiLocal, nodes); 
        }
      }
    }
    
    if (matched)
    {
      break;
    }
  }

  // if not matched, insert the new nodes  
  return matched;
}

bool PerforationData::SearchEntireDomain( MeshLevel  const  & mesh,
                                          R1Tensor   const  & location,
                                          localIndex        & erMatched,
                                          localIndex        & esrMatched,
                                          localIndex        & eiMatched ) const
{
  ElementRegionManager const * elemManager = mesh.getElemManager();
  NodeManager const * const nodeManager    = mesh.getNodeManager();
  
  bool matched = false;

  // loop over the cell element regions
  elemManager->forElementSubRegionsComplete<CellElementSubRegion>( [&]( localIndex const er,
                                                                        localIndex const esr,
                                                                        ElementRegionBase const * GEOSX_UNUSED_ARG( elemRegion ), 
                                                                        CellElementSubRegion const * subRegion )
  {
    // loop over the elements of the subregion 
    for (localIndex ei = 0; ei < subRegion->size(); ++ei) 
    {                              
      // check if the element of the subregion contains the perforations
      if (IsPointInsideElement( nodeManager, location, subRegion, ei ))
      {
        // store the indices of the mesh element
        erMatched  = er;
        esrMatched = esr;
        eiMatched  = ei;
        matched = true;
        break;
      }
    }
  });
  return matched;
}

bool PerforationData::IsPointInsideElement( NodeManager const * const nodeManager,
                                            R1Tensor    const & location,
                                            CellBlock   const * subRegion,
                                            localIndex          ei ) const
{
  bool isInsideElement = false;

  array1d<array1d<localIndex>> faceNodes( subRegion->numFacesPerElement() );

  // collect the faces for this element
  for (localIndex kf = 0; kf < subRegion->numFacesPerElement(); ++kf)
  {
    subRegion->GetFaceNodes( ei, kf, faceNodes[kf] );
  }

  // if the point is in the element, save the indices and stop the search              
  if (computationalGeometry::IsPointInsidePolyhedron( nodeManager->referencePosition(),
                                                      faceNodes,
                                                      location ))
  {
    isInsideElement = true;
  }
  return isInsideElement;
} 

void PerforationData::CollectNodes( CellBlock const * subRegion,
                                    localIndex        ei,
                                    set<localIndex> & nodes ) const
{
  // get all the nodes belonging to this element 
  for (localIndex a = 0; a < subRegion->numNodesPerElement(); ++a)
  {
    localIndex const inode = subRegion->nodeList(ei,a);

    // if not already visited, store the newly found node
    if (!nodes.contains(inode))
    {
      nodes.insert(inode);
    }
  }  
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
    GEOS_ASSERT( globalToLocalWellElemMap.count( ielemGlobal ) > 0 );
    m_wellElementIndex[iperfLocal] = globalToLocalWellElemMap.at( ielemGlobal ); 
  }

  //DebugLocalPerforations();
}
 
void PerforationData::InitializePostInitialConditions_PreSubGroups( Group * const GEOSX_UNUSED_ARG( problemManager ) )
{
  for (localIndex iperf = 0; iperf < size(); ++iperf)
  {
    if (m_wellPeacemanIndex[iperf] < 0.0)
    {
      // TODO: compute wellPeacemanIndex internally
      GEOS_ERROR( "Invalid well Peaceman index value: " << m_wellPeacemanIndex[iperf] );
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
    std::cout << "m_wellPeacemanIndex[" << iperf << "] = " << m_wellPeacemanIndex[iperf]
              << std::endl;
  }
}

} //namespace geosx
