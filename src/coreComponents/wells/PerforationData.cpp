/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * @file PerforationData.cpp
 *
 */

#include "PerforationData.hpp"
#include "Perforation.hpp"
#include "Well.hpp"

#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"

namespace geosx
{

using namespace dataRepository;

PerforationData::PerforationData(string const & name, ManagedGroup * const parent)
  : ObjectManagerBase(name, parent)
{
  RegisterViewWrapper( viewKeyStruct::reservoirElementRegionString, &m_reservoirElementRegion, false );
  RegisterViewWrapper( viewKeyStruct::reservoirElementSubregionString, &m_reservoirElementSubregion, false );
  RegisterViewWrapper( viewKeyStruct::reservoirElementIndexString, &m_reservoirElementIndex, false );

  RegisterViewWrapper( viewKeyStruct::wellElementIndexString, &m_wellElementIndex, false );

  RegisterViewWrapper( viewKeyStruct::locationString, &m_location, false );
  RegisterViewWrapper( viewKeyStruct::transmissibilityString, &m_transmissibility, false );
  RegisterViewWrapper( viewKeyStruct::gravityDepthString, &m_gravityDepth, false );
}

PerforationData::~PerforationData()
{
}

const string PerforationData::getCatalogName() const
{
  return keys::perforationData;
}

ManagedGroup * PerforationData::CreateChild(string const & childKey, string const & childName)
{
  return nullptr;
}
  
void PerforationData::InitializePreSubGroups( ManagedGroup * const problemManager )
{
  DomainPartition const * domain
    = problemManager->GetGroup<DomainPartition>( keys::domain );
  MeshLevel const * mesh = domain->getMeshBody(0)->getMeshLevel(0);
  
  // TODO: trace the reservoir cells that are traversed by the well
  ConnectToCells( mesh );
}

localIndex PerforationData::numPerforationsGlobal() const
{
  Well const * parent = getParent()->group_cast<Well const *>();
  PerforationManager const * perforationManager
    = parent->GetGroup<PerforationManager>( Well::groupKeyStruct::perforationsString );
  return perforationManager->numPerforationsGlobal();
}
  
void PerforationData::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager )
{
}

void PerforationData::ConnectToCells( MeshLevel const * mesh )
{
  ElementRegionManager const * elemManager = mesh->getElemManager();

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor const>> elemCenter =
    elemManager->ConstructViewAccessor<array1d<R1Tensor>, arrayView1d<R1Tensor const>>( ElementSubRegionBase::
                                                                                        viewKeyStruct::
                                                                                        elementCenterString );
  // TODO: trace the reservoir cells that are traversed by the well
  // TODO: rewrite this function entirely

  // get the managers for perforations and well elements 
  PerforationManager const * perforationManager
    = getParent()->GetGroup<PerforationManager>( Well::
						 groupKeyStruct::
						 perforationsString );
  WellElementManager const * wellElementManager
    = getParent()->GetGroup<WellElementManager>( Well::
						 groupKeyStruct::
						 wellElementsString );

  // currently, assume that all perforation are on this rank
  resize( perforationManager->numPerforationsGlobal() );
  
  // TODO Until we can properly trace perforations to cells,
  // just connect to the nearest cell center (this is NOT correct in general)
  localIndex num_conn_local  = 0;
  localIndex num_conn_global = 0;

  // loop over all the perforations
  for ( globalIndex iperf = 0; iperf < perforationManager->numPerforationsGlobal(); ++iperf )
  {
    // get the physical location of the current perforation
    Perforation const * const perforation = perforationManager->getPerforation( iperf );
    R1Tensor const & loc = perforation->getLocation();

    // find the closest reservoir element
    auto ret = minLocOverElemsInMesh( mesh, [&] ( localIndex er,
                                                  localIndex esr,
                                                  localIndex ei ) -> real64
    {
      R1Tensor v = loc;
      v -= elemCenter[er][esr][ei];
      return v.L2_Norm();
    });
    
    // save the region, subregion and index 
    m_reservoirElementRegion[num_conn_local]    = std::get<0>(ret.second);
    m_reservoirElementSubregion[num_conn_local] = std::get<1>(ret.second);
    m_reservoirElementIndex[num_conn_local]     = std::get<2>(ret.second);

    // find the well element of this perforation by name
    string const wellElementName = perforation->getWellElementName();
    m_wellElementIndex[num_conn_local] = -1;
 
    for (localIndex iwelem = 0; iwelem < wellElementManager->numWellElementsGlobal(); ++iwelem)
    {
      WellElement const * const wellElem = wellElementManager->getWellElement( iwelem );
 
      // save the index if the names match
      if (wellElem->getName() == wellElementName)
      {
	m_wellElementIndex[num_conn_local] = iwelem;
	break;
      }
    }

    // error message if the well element is not found
    if (m_wellElementIndex[num_conn_local] == -1)
    {
      GEOS_ERROR("Invalid well element name: " << wellElementName
                 << " for perforation " << perforation->getName() );
    }    
    
    // increment the index, save location and transmissibility
    m_location[num_conn_local]         = perforation->getLocation();
    m_transmissibility[num_conn_local] = perforation->getTransmissibility();
    
    num_conn_local++;
  }
  resize( num_conn_local );
}

} //namespace geosx
