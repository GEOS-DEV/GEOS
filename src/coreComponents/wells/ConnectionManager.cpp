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
 * @file ConnectionManager.cpp
 *
 */

#include "ConnectionManager.hpp"
#include "Connection.hpp"
#include "WellBase.hpp"

#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"

namespace geosx
{

using namespace dataRepository;

ConnectionManager::ConnectionManager(string const & name, ManagedGroup * const parent)
  : ObjectManagerBase(name, parent)
{
  RegisterViewWrapper( viewKeyStruct::connIndexString, &m_connIndex, false );
}

ConnectionManager::~ConnectionManager()
{

}

ManagedGroup * ConnectionManager::CreateChild(string const & childKey, string const & childName)
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

const string ConnectionManager::getCatalogName() const
{
  return keys::connections;
}

void ConnectionManager::InitializePreSubGroups( ManagedGroup * const problemManager )
{

}

void ConnectionManager::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager )
{
  DomainPartition const * domain = problemManager->GetGroup<DomainPartition>( keys::domain );
  MeshLevel const * mesh = domain->getMeshBody(0)->getMeshLevel(0);
  PrecomputeData( mesh );
}

void ConnectionManager::PrecomputeData( MeshLevel const * mesh )
{

}

Connection const * ConnectionManager::getConnection( localIndex iconn ) const
{
  return this->GetGroup<Connection>( m_connectionList[iconn] );
}

Connection * ConnectionManager::getConnection( localIndex iconn )
{
  return this->GetGroup<Connection>( m_connectionList[iconn] );
}


} //namespace geosx
