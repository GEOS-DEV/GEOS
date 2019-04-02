/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

#include "TasksManager.hpp"
#include "TaskBase.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

TasksManager::TasksManager( std::string const & name,
                            ManagedGroup * const parent ):
  ManagedGroup( name, parent)
{
}

TasksManager::~TasksManager()
{}


ManagedGroup * TasksManager::CreateChild( string const & childKey, string const & childName )
{
  std::unique_ptr<TaskBase> tool = TaskBase::CatalogInterface::Factory( childKey, childName, this );
  return this->RegisterGroup( childName, std::move( tool ) );
}


void TasksManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from SolverBase here
  for (auto& catalogIter: TaskBase::GetCatalog())
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}


} /* namespace geosx */
