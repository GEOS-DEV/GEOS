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
 * @file ConnectionData.cpp
 *
 */

#include "ConnectionData.hpp"
#include "Connection.hpp"
#include "Well.hpp"

#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"

namespace geosx
{

using namespace dataRepository;

ConnectionData::ConnectionData(string const & name, ManagedGroup * const parent)
  : ObjectManagerBase(name, parent)
{
  RegisterViewWrapper( viewKeyStruct::connectionIndexString, &m_connectionIndex, false );
}

ConnectionData::~ConnectionData()
{

}

const string ConnectionData::getCatalogName() const
{
  return keys::connectionData;
}
  
ManagedGroup * ConnectionData::CreateChild(string const & childKey, string const & childName)
{
  if ( childKey == groupKeyStruct::connectionString )
  {
    m_connectionList.push_back( childName );
    return RegisterGroup<Connection>( childName );
  }
  else
  {
    GEOS_ERROR( "Unrecognized node: " << childKey );
  }
  return nullptr;
}

void ConnectionData::InitializePreSubGroups( ManagedGroup * const problemManager )
{
  // for now, assume that numConnectionsGlobal == numConnectionsLocal
  resize( m_connectionList.size() );

  // assume single-segmented well
  m_nextWellElementIndex[0] = 0;
  m_prevWellElementIndex[0] = -1;
}

void ConnectionData::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager )
{
  DomainPartition const * domain = problemManager->GetGroup<DomainPartition>( keys::domain );
  MeshLevel const * mesh = domain->getMeshBody(0)->getMeshLevel(0);
  PrecomputeData( mesh );
}

void ConnectionData::PrecomputeData( MeshLevel const * mesh )
{

}

Connection const * ConnectionData::getConnection( localIndex iconn ) const
{
  return this->GetGroup<Connection>( m_connectionList[iconn] );
}

Connection * ConnectionData::getConnection( localIndex iconn )
{
  return this->GetGroup<Connection>( m_connectionList[iconn] );
}


} //namespace geosx
