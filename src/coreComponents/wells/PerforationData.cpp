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
  RegisterViewWrapper( viewKeyStruct::perforationIndexString, &m_perforationIndex, false );
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

Perforation const * PerforationData::getPerforation( localIndex iperf ) const
{
  Well const * parent = getParent()->group_cast<Well const *>();
  PerforationManager const * perforationManager
    = parent->GetGroup<PerforationManager>( Well::groupKeyStruct::perforationsString );
  Perforation const * perforation
    = perforationManager->getPerforation( m_perforationIndex[iperf] );
  return perforation;
}

Perforation * PerforationData::getPerforation( localIndex iperf ) 
{
  Well * parent = getParent()->group_cast<Well *>();
  PerforationManager * perforationManager
    = parent->GetGroup<PerforationManager>( Well::groupKeyStruct::perforationsString );
  Perforation * perforation
    = perforationManager->getPerforation( m_perforationIndex[iperf] );
  return perforation;
}
  
void PerforationData::InitializePreSubGroups( ManagedGroup * const problemManager )
{
  std::cout << "PerforationData::InitializePreSubGroups started" << std::endl;
  DomainPartition const * domain
    = problemManager->GetGroup<DomainPartition>( keys::domain );
  MeshLevel const * mesh = domain->getMeshBody(0)->getMeshLevel(0);
  ConnectToCells( mesh );
  std::cout << "PerforationData::InitializePreSubGroups complete" << std::endl;
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
  R1Tensor const & gravity = getParent()->group_cast<Well *>()->getGravityVector();
  arrayView1d<real64> & gravDepth = getReference<array1d<real64>>( viewKeyStruct::gravityDepthString );

  for (localIndex iperf = 0; iperf < size(); ++iperf)
  {
    Perforation const * perforation = getPerforation( iperf );
    gravDepth[iperf] = Dot( perforation->getLocation(), gravity );
  }
}

void PerforationData::ConnectToCells( MeshLevel const * mesh )
{
  ElementRegionManager const * elemManager = mesh->getElemManager();

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor const>> elemCenter =
    elemManager->ConstructViewAccessor<array1d<R1Tensor>, arrayView1d<R1Tensor const>>( ElementSubRegionBase::
                                                                                        viewKeyStruct::
                                                                                        elementCenterString );

  // initially allocate enough memory for all (global) perforations
  PerforationManager const * perforationManager
    = getParent()->GetGroup<PerforationManager>( Well::
						 groupKeyStruct::
						 perforationsString );
  WellElementManager const * wellElementManager
    = getParent()->GetGroup<WellElementManager>( Well::
						 groupKeyStruct::
						 wellElementsString );

  resize( perforationManager->numPerforationsGlobal() );
  
  // TODO Until we can properly trace perforations to cells,
  // just connect to the nearest cell center (this is NOT correct in general)
  localIndex num_conn_local  = 0;
  localIndex num_conn_global = 0;

  for ( globalIndex iperf = 0; iperf < perforationManager->numPerforationsGlobal(); ++iperf )
  {
    Perforation const * perforation = perforationManager->getPerforation( iperf );
    R1Tensor const & loc = perforation->getLocation();

    auto ret = minLocOverElemsInMesh( mesh, [&] ( localIndex er,
                                                  localIndex esr,
                                                  localIndex ei ) -> real64
    {
      R1Tensor v = loc;
      v -= elemCenter[er][esr][ei];
      return v.L2_Norm();
    });
    
    m_reservoirElementRegion[num_conn_local]    = std::get<0>(ret.second);
    m_reservoirElementSubregion[num_conn_local] = std::get<1>(ret.second);
    m_reservoirElementIndex[num_conn_local]     = std::get<2>(ret.second);

    string const wellElementName = perforation->getWellElementName();

    // TODO: rewrite this entirely
    m_wellElementIndex[num_conn_local] = -1;
    for (localIndex iwelem = 0; iwelem < wellElementManager->numWellElementsGlobal(); ++iwelem)
    {
      WellElement const * wellElem = wellElementManager->getWellElement( iwelem );
      if (wellElem->getName() == wellElementName)
      {
	m_wellElementIndex[num_conn_local] = iwelem;
	break;
      }
    }
    if (m_wellElementIndex[num_conn_local] == -1)
      GEOS_ERROR("Invalid well element name: " << wellElementName);
    
    std::cout << "I match perforation " << perforation->getName()
	      << " with segment " << wellElementManager->getWellElement( m_wellElementIndex[num_conn_local] )->getName()
	      << std::endl;
    
    m_perforationIndex[num_conn_local] = num_conn_global++;
    m_transmissibility[num_conn_local] = perforation->getTransmissibility();
    
    num_conn_local++;
  }
  resize( num_conn_local );
}

} //namespace geosx
