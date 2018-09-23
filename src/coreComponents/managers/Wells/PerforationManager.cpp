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
 * @file PerforationManager.cpp
 *
 */

#include "PerforationManager.hpp"
#include "Perforation.hpp"
#include "WellBase.hpp"

#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"

namespace geosx
{

using namespace dataRepository;

class WellConnectionManager
{

};

PerforationManager::PerforationManager(string const & name, ManagedGroup * const parent)
  : ObjectManagerBase(name, parent)
{
  RegisterViewWrapper( viewKeysPerfManager.connectionElementRegion.Key(), &m_connectionElementRegion, false );
  RegisterViewWrapper( viewKeysPerfManager.connectionElementSubregion.Key(), &m_connectionElementSubregion, false );
  RegisterViewWrapper( viewKeysPerfManager.connectionElementIndex.Key(), &m_connectionElementIndex, false );
  RegisterViewWrapper( viewKeysPerfManager.connectionPerforationIndex.Key(), &m_connectionPerforationIndex, false );
}

PerforationManager::~PerforationManager()
{

}

void PerforationManager::FillDocumentationNode()
{
  ObjectManagerBase::FillDocumentationNode();

  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName( keys::perforations );
  docNode->setSchemaType( "UniqueNode" );

  docNode->AllocateChildNode( viewKeysPerfManager.gravityDepth.Key(),
                              viewKeysPerfManager.gravityDepth.Key(),
                              -1,
                              "real64_array",
                              "real64_array",
                              "Connection gravity-depth product",
                              "Connection gravity-depth product",
                              "",
                              "",
                              1,
                              0,
                              1 );
}

void PerforationManager::CreateChild(string const & childKey, string const & childName)
{
  if ( childKey == groupKeyStruct::perforationString )
  {
    this->RegisterGroup<Perforation>( childName );
    m_allPerfList.push_back( childName );
  }
  else
  {
    GEOS_ERROR("Unrecognized node: " << childKey);
  }
}

const string PerforationManager::getCatalogName() const
{
  return keys::perforations;
}

void PerforationManager::FinalInitialization(ManagedGroup * const problemManager)
{
  DomainPartition const * domain = problemManager->GetGroup<DomainPartition>( keys::domain );

  // initially allocate enough memory for all (global) perforations
  resize( numConnectionsGlobal() );
  ConnectToCells( domain );
  resize( numConnectionsLocal() );

  // generate the "all" set to enable application of BC
  ManagedGroup * sets = GetGroup( keys::sets );
  set<localIndex> & setAll = sets->RegisterViewWrapper<set<localIndex>>("all")->reference();
  for (localIndex iconn = 0; iconn < numConnectionsLocal(); ++iconn)
  {
    setAll.insert( iconn );
  }

  PrecomputeData( domain );
}

void PerforationManager::PrecomputeData(DomainPartition const * domain)
{
  R1Tensor const & gravity = getParent()->group_cast<WellBase *>()->getGravityVector();
  array1d<real64> & gravDepth = getReference<array1d<real64>>(viewKeysPerfManager.gravityDepth);

  for (localIndex iconn = 0; iconn < size(); ++iconn)
  {
    string const & perfName = m_allPerfList[m_connectionPerforationIndex[iconn]];
    Perforation const * perf = this->GetGroup<Perforation>( perfName );
    gravDepth[iconn] = Dot( perf->getLocation(), gravity );
  }
}

void PerforationManager::ConnectToCells(DomainPartition const * domain)
{
  MeshLevel const * mesh = domain->getMeshBody(0)->getMeshLevel(0);
  ElementRegionManager const * elemManager = mesh->getElemManager();

  auto elemCenter = elemManager->ConstructViewAccessor<array1d<R1Tensor>>( CellBlockSubRegion::
                                                                           viewKeyStruct::
                                                                           elementCenterString );

  // TODO Until we can properly trace perforations to cells,
  // just connect to the nearest cell center (this is NOT correct in general)

  localIndex iconn_local = 0;
  localIndex iconn_global = 0;
  for ( string const & perfName : m_allPerfList )
  {
    Perforation const * perf = GetGroup<Perforation>( perfName );
    R1Tensor const & loc = perf->getLocation();

    auto ret = minLocOverElemsInMesh( mesh, [&] ( localIndex er,
                                                  localIndex esr,
                                                  localIndex ei ) -> real64
    {
      R1Tensor v = loc;
      v -= elemCenter[er][esr][ei];
      return v.L2_Norm();
    });

    m_connectionElementRegion[iconn_local]    = std::get<0>(ret.second);
    m_connectionElementSubregion[iconn_local] = std::get<1>(ret.second);
    m_connectionElementIndex[iconn_local]     = std::get<2>(ret.second);
    m_connectionPerforationIndex[iconn_local] = iconn_global++;

    // This will not be correct in parallel until we can actually check that
    // the perforation belongs to local mesh partition
    ++iconn_local;
  }
}

Perforation const * PerforationManager::getPerforation(localIndex iperf) const
{
  return this->GetGroup<Perforation>( m_allPerfList[iperf] );
}

Perforation * PerforationManager::getPerforation(localIndex iperf)
{
  return this->GetGroup<Perforation>( m_allPerfList[iperf] );
}


} //namespace geosx
