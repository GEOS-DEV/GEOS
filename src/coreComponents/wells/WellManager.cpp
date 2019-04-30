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
 * @file WellManager.cpp
 *
 */

#include "WellManager.hpp"
#include "Well.hpp"

namespace geosx
{

using namespace dataRepository;


WellManager::WellManager(string const & name,
                         dataRepository::ManagedGroup * const parent)
  : dataRepository::ManagedGroup(name, parent)
{

  RegisterViewWrapper( viewKeyStruct::materialListString, &m_materialList, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("List of materials present in the well region"); 

}

WellManager::~WellManager()
{
}
  
ManagedGroup * WellManager::CreateChild(string const & childKey, string const & childName)
{
  if ( childKey == "Well")
  {
    return RegisterGroup<Well>( childName );
  }
  else
  {
    GEOS_ERROR( "Unrecognized node: " << childKey );
  }
  return nullptr;
}

} //namespace geosx
