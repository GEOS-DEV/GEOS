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
 * @file WellElementManager.cpp
 *
 */

#include "WellElementManager.hpp"
#include "WellElement.hpp"

namespace geosx
{

using namespace dataRepository;

WellElementManager::WellElementManager(string const & name, ManagedGroup * const parent)
  : ObjectManagerBase(name, parent)
{
}

WellElementManager::~WellElementManager()
{
}

const string WellElementManager::getCatalogName() const
{
  return keys::wellElements;
}

ManagedGroup * WellElementManager::CreateChild(string const & childKey, string const & childName)
{
  if ( childKey == groupKeyStruct::wellElementString )
  {
    m_globalWellElementList.push_back( childName );
    return RegisterGroup<WellElement>( childName );
  }
  else
  {
    GEOS_ERROR( "Unrecognized node: " << childKey );
  }
  return nullptr;
}

WellElement const * WellElementManager::getWellElement( globalIndex iwelem ) const
{
  return GetGroup<WellElement>( m_globalWellElementList[iwelem] );
}

WellElement * WellElementManager::getWellElement( globalIndex iwelem )
{
  return GetGroup<WellElement>( m_globalWellElementList[iwelem] );
}

  
void WellElementManager::InitializePreSubGroups( ManagedGroup * const problemManager )
{
}

void WellElementManager::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager )
{
}

} //namespace geosx
