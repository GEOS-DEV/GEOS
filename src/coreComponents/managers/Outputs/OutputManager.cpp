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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file OutputManager.cpp
 */

#include "OutputManager.hpp"
#include "DocumentationNode.hpp"
#include "SiloOutput.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

OutputManager::OutputManager( std::string const & name,
                              ManagedGroup * const parent ):
  ManagedGroup( name, parent)
{
}

OutputManager::~OutputManager()
{}



ManagedGroup * OutputManager::CreateChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0("Adding Output: " << childKey << ", " << childName);
  std::unique_ptr<OutputBase> output = OutputBase::CatalogInterface::Factory( childKey, childName, this );
  return this->RegisterGroup<OutputBase>( childName, std::move(output) );
}


} /* namespace geosx */
