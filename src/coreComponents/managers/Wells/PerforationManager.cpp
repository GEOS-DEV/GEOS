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
  RegisterViewWrapper( viewKeyStruct::connectionElementRegionString, &m_connectionElementRegion, false );
  RegisterViewWrapper( viewKeyStruct::connectionElementSubregionString, &m_connectionElementSubregion, false );
  RegisterViewWrapper( viewKeyStruct::connectionElementIndexString, &m_connectionElementIndex, false );
  RegisterViewWrapper( viewKeyStruct::connectionPerforationIndexString, &m_connectionPerforationIndex, false );

  RegisterViewWrapper( viewKeyStruct::gravityDepthString, &m_gravityDepth, false );
}

PerforationManager::~PerforationManager()
{

}

ManagedGroup * PerforationManager::CreateChild(string const & childKey, string const & childName)
{
  if ( childKey == groupKeyStruct::perforationString )
  {
    m_allPerfList.push_back( childName );
    return RegisterGroup<Perforation>( childName );
  }
  else
  {
    GEOS_ERROR( "Unrecognized node: " << childKey );
  }
  return nullptr;
}

const string PerforationManager::getCatalogName() const
{
  return keys::perforations;
}

void PerforationManager::InitializePreSubGroups( ManagedGroup * const problemManager )
{
  DomainPartition const * domain = problemManager->GetGroup<DomainPartition>( keys::domain );
  MeshLevel const * mesh = domain->getMeshBody(0)->getMeshLevel(0);

  ConnectToCells( mesh );

  // generate the "all" set to enable application of BC
  set<localIndex> & setAll = this->sets()->RegisterViewWrapper<set<localIndex>>("all")->reference();
  for (localIndex iconn = 0; iconn < numConnectionsLocal(); ++iconn)
  {
    setAll.insert( iconn );
  }
}

void PerforationManager::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager )
{
  DomainPartition const * domain = problemManager->GetGroup<DomainPartition>( keys::domain );
  MeshLevel const * mesh = domain->getMeshBody(0)->getMeshLevel(0);
  PrecomputeData( mesh );
}

void PerforationManager::PrecomputeData( MeshLevel const * mesh )
{
  R1Tensor const & gravity = getParent()->group_cast<WellBase *>()->getGravityVector();
  arrayView1d<real64> & gravDepth = getReference<array1d<real64>>( viewKeyStruct::gravityDepthString );

  for (localIndex iconn = 0; iconn < size(); ++iconn)
  {
    string const & perfName = m_allPerfList[m_connectionPerforationIndex[iconn]];
    Perforation const * perf = this->GetGroup<Perforation>( perfName );
    gravDepth[iconn] = Dot( perf->getLocation(), gravity );
  }
}

void PerforationManager::ConnectToCells( MeshLevel const * mesh )
{
  ElementRegionManager const * elemManager = mesh->getElemManager();

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor const>> elemCenter =
    elemManager->ConstructViewAccessor<array1d<R1Tensor>, arrayView1d<R1Tensor const>>( CellBlockSubRegion::
                                                                                        viewKeyStruct::
                                                                                        elementCenterString );

  // initially allocate enough memory for all (global) perforations
  resize( numConnectionsGlobal() );

  // TODO Until we can properly trace perforations to cells,
  // just connect to the nearest cell center (this is NOT correct in general)
  localIndex num_conn_local = 0;
  localIndex num_conn_global = 0;
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

    m_connectionElementRegion[num_conn_local]    = std::get<0>(ret.second);
    m_connectionElementSubregion[num_conn_local] = std::get<1>(ret.second);
    m_connectionElementIndex[num_conn_local]     = std::get<2>(ret.second);
    m_connectionPerforationIndex[num_conn_local] = num_conn_global++;

    // This will not be correct in parallel until we can actually check that
    // the perforation belongs to local mesh partition
    ++num_conn_local;
  }

  resize( num_conn_local );
}

Perforation const * PerforationManager::getPerforation( localIndex iperf ) const
{
  return this->GetGroup<Perforation>( m_allPerfList[iperf] );
}

Perforation * PerforationManager::getPerforation( localIndex iperf )
{
  return this->GetGroup<Perforation>( m_allPerfList[iperf] );
}


} //namespace geosx
