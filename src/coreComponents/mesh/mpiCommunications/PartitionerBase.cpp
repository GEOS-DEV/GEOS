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
 * @file PartitionerBase.cpp
 */


#include "PartitionerBase.hpp"


namespace geos
{
using namespace dataRepository;

PartitionerBase::PartitionerBase( string const & name,
                                  Group * const parent ):
  Group( name, parent ),
  m_numPartitions( 1 ), // Default to a single partition
  m_numColors( 1 ),     // Default to a single color
  m_color( 0 )          // Default to color 0
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE ); // needed for schema generation
}

PartitionerBase::CatalogInterface::CatalogType & PartitionerBase::getCatalog()
{
  static PartitionerBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

} // namespace geos
