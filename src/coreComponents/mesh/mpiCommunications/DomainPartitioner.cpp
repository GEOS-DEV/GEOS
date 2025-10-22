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
 * @file DomainPartitioner.cpp
 */

#include "DomainPartitioner.hpp"

namespace geos
{
using namespace dataRepository;

DomainPartitioner::DomainPartitioner( string const & name,
                                      Group * const parent )
  : Group( name, parent ),
  m_neighborsRank(),
  m_numPartitions( 1 ),
  m_numColors( 1 ),
  m_color( 0 )
{
  {
    setInputFlags( InputFlags::OPTIONAL_NONUNIQUE ); // needed for schema generation
    setRestartFlags( RestartFlags::NO_WRITE ); // partitioners are configuration-only, not simulation state
  }

}

DomainPartitioner::~DomainPartitioner()
{}

void DomainPartitioner::postInputInitialization()
{
  Group::postInputInitialization();
}

void DomainPartitioner::setPartitionCounts( unsigned int const xPartitions,
                                            unsigned int const yPartitions,
                                            unsigned int const zPartitions )
{
  m_numPartitions = xPartitions * yPartitions * zPartitions;

  GEOS_ERROR_IF( m_numPartitions < 1,
                 "Total partition count must be at least 1" );
}

void DomainPartitioner::setNeighborsRank( std::vector< int > const & neighborsRank )
{
  m_neighborsRank = neighborsRank;
}

DomainPartitioner::CatalogInterface::CatalogType &
DomainPartitioner::getCatalog()
{
  static DomainPartitioner::CatalogInterface::CatalogType catalog;
  return catalog;
}

void DomainPartitioner::printInfo() const
{
  GEOS_LOG_RANK_0( getInfoString() );
}

string DomainPartitioner::getInfoString() const
{
  return GEOS_FMT( "{} '{}': partitions={}, neighbors={}, colors={}",
                   getCatalogName(),
                   getName(),
                   m_numPartitions,
                   m_neighborsRank.size(),
                   m_numColors );
}

} // namespace geos
