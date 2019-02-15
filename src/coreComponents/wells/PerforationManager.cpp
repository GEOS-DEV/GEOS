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

namespace geosx
{

using namespace dataRepository;

PerforationManager::PerforationManager(string const & name, ManagedGroup * const parent)
  : ObjectManagerBase(name, parent)
{
}

PerforationManager::~PerforationManager()
{
}

const string PerforationManager::getCatalogName() const
{
  return keys::perforations;
}

ManagedGroup * PerforationManager::CreateChild(string const & childKey, string const & childName)
{
  if ( childKey == groupKeyStruct::perforationString )
  {
    m_globalPerforationList.push_back( childName );
    return RegisterGroup<Perforation>( childName );
  }
  else
  {
    GEOS_ERROR( "Unrecognized node: " << childKey );
  }
  return nullptr;
}

Perforation const * PerforationManager::getPerforation( globalIndex iperf ) const
{
  return GetGroup<Perforation>( m_globalPerforationList[iperf] );
}

Perforation * PerforationManager::getPerforation( globalIndex iperf )
{
  return GetGroup<Perforation>( m_globalPerforationList[iperf] );
}
  
void PerforationManager::InitializePreSubGroups( ManagedGroup * const problemManager )
{
}
  
void PerforationManager::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager )
{
}

} //namespace geosx
