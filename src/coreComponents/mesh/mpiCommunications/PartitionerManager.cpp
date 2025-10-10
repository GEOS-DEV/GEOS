/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 TheBoard of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PartitionerManager.cpp
 */

#include "PartitionerManager.hpp"

namespace geos
{

using namespace dataRepository;

PartitionerManager::PartitionerManager( string const & name,
                                        Group * const parent )
  : Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
}

PartitionerManager::~PartitionerManager()
{}

Group * PartitionerManager::createChild( string const & childKey,
                                         string const & childName )
{
  GEOS_LOG_RANK_0( GEOS_FMT( "PartitionerManager: adding {} '{}'", childKey, childName ) );

  // Create partitioner from DomainPartitioner catalog
  std::unique_ptr< DomainPartitioner > partitioner =
    DomainPartitioner::CatalogInterface::factory( childKey,
                                                  getDataContext(),
                                                  childName,
                                                  this );

  return &this->registerGroup< DomainPartitioner >( childName, std::move( partitioner ) );
}

void PartitionerManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from DomainPartitioner
  for( auto & catalogIter : DomainPartitioner::getCatalog() )
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}

bool PartitionerManager::hasPartitioner() const
{
  return this->numSubGroups() > 0;
}

DomainPartitioner & PartitionerManager::getPartitioner()
{
  // Get the first (and should be only) partitioner defined in the XML
  GEOS_ERROR_IF( !hasPartitioner(),
                 "No partitioner defined in XML. Please add a partitioner inside <Partitioner>.\n"
                 "Examples:\n"
                 "  <CellGraphPartitioner name=\"partitioner\" engine=\"parmetis\"/>\n"
                 "  <CartesianPartitioner name=\"partitioner\" partitionCounts=\"{4, 4, 2}\"/>" );

  GEOS_ERROR_IF( this->numSubGroups() > 1,
                 "Multiple partitioners defined. Only one partitioner can be active at a time." );

  // Return the first child which should be a DomainPartitioner
  return this->getGroup< DomainPartitioner >( 0 );
}

DomainPartitioner const & PartitionerManager::getPartitioner() const
{
  // Get the first (and should be only) partitioner defined in the XML
  GEOS_ERROR_IF( this->numSubGroups() == 0,
                 "No partitioner defined in XML. Please add a partitioner inside <Partitioner>.\n"
                 "Examples:\n"
                 "  <CellGraphPartitioner name=\"partitioner\" engine=\"parmetis\"/>\n"
                 "  <CartesianPartitioner name=\"partitioner\" />" );

  GEOS_ERROR_IF( this->numSubGroups() > 1,
                 "Multiple partitioners defined. Only one partitioner can be active at a time." );

  // Return the first child which should be a DomainPartitioner
  return this->getGroup< DomainPartitioner >( 0 );
}

} // namespace geos
