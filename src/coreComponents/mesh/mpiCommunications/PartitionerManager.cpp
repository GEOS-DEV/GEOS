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


#include "PartitionerManager.hpp"


namespace geos
{

using namespace dataRepository;

PartitionerManager::PartitionerManager( string const & name,
                                        Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
}

PartitionerManager::~PartitionerManager()
{}

Group * PartitionerManager::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( GEOS_FMT( "{}: adding {} {}", getName(), childKey, childName ) );
  std::unique_ptr< PartitionerBase > partitioner =
    PartitionerBase::CatalogInterface::factory( childKey, getDataContext(), childName, this );
  return &this->registerGroup< PartitionerBase >( childName, std::move( partitioner ) );
}


void PartitionerManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from PartitionerBase here
  for( auto & catalogIter: PartitionerBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}

bool PartitionerManager::hasPartitioner() const
{
  return this->numSubGroups() > 0;
}

PartitionerBase & PartitionerManager::getPartitioner()
{
  // Get the first (and should be only) partitioner defined in the XML
  GEOS_ERROR_IF( this->numSubGroups() == 0,
                 "No partitioner defined in XML. Please add a <Partitioner> node with a partitioner type (e.g., <Cartesian name=\"partitioner\"/>)" );

  GEOS_ERROR_IF( this->numSubGroups() > 1,
                 "Multiple partitioners defined. Only one partitioner can be active at a time." );

  // Return the first child which should be a PartitionerBase
  return this->getGroup< PartitionerBase >( 0 );
}

PartitionerBase const & PartitionerManager::getPartitioner() const
{
  // Get the first (and should be only) partitioner defined in the XML
  GEOS_ERROR_IF( this->numSubGroups() == 0,
                 "No partitioner defined in XML. Please add a <Partitioner> node with a partitioner type (e.g., <Cartesian name=\"partitioner\"/>)" );

  GEOS_ERROR_IF( this->numSubGroups() > 1,
                 "Multiple partitioners defined. Only one partitioner can be active at a time." );

  // Return the first child which should be a PartitionerBase
  return this->getGroup< PartitionerBase >( 0 );
}

} // namespace geos
