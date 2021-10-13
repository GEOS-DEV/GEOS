/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TasksManager.cpp
 */

#include "TasksManager.hpp"
#include "TaskBase.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace LvArray;

TasksManager::TasksManager( string const & name,
                            Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
}

TasksManager::~TasksManager()
{ }

Group * TasksManager::createChild( string const & childKey, string const & childName )
{
  std::unique_ptr< TaskBase > task = TaskBase::CatalogInterface::factory( childKey, childName, this );
  return &this->registerGroup< TaskBase >( childName, std::move( task ) );
}

void TasksManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from TaskBase here
  for( auto & catalogIter: TaskBase::getCatalog() )
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}

} /* namespace geosx */
