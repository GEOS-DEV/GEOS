/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file GraphPartitionEngine.cpp
 */

#include "GraphPartitionEngine.hpp"
#ifdef GEOS_USE_PARMETIS
  #include "ParMetisEngine.hpp"
#endif
namespace geos
{
using namespace dataRepository;

GraphPartitionEngine::GraphPartitionEngine( string const & name,
                                            dataRepository::Group * const parent )
  : Group( name, parent )
{
  // No wrappers to register
  // Only set input flags if we have a parent (i.e., when used in production code)
  // For unit tests with nullptr parent, skip this
  if( parent != nullptr )
  {
    setInputFlags( dataRepository::InputFlags::FALSE );
    setRestartFlags( RestartFlags::NO_WRITE );
  }
}

GraphPartitionEngine::~GraphPartitionEngine()
{
  // Default destructor
}

GraphPartitionEngine::CatalogInterface::CatalogType &
GraphPartitionEngine::getCatalog()
{
  static GraphPartitionEngine::CatalogInterface::CatalogType catalog;
  return catalog;
}

} // namespace geos
