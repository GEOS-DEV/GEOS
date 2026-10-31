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
 * @file MixedMimeticDiscretizationManager.cpp
 */

#include "MixedMimeticDiscretizationManager.hpp"

#include "mixedMimetic/MixedMimeticDiscretization.hpp"

namespace geos
{

using namespace dataRepository;

MixedMimeticDiscretizationManager::MixedMimeticDiscretizationManager( string const & name, Group * const parent )
  : Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
}

MixedMimeticDiscretizationManager::~MixedMimeticDiscretizationManager()
{}

Group * MixedMimeticDiscretizationManager::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( GEOS_FMT( "{}: adding {} {}", getName(), childKey, childName ) );
  std::unique_ptr< MixedMimeticDiscretization > disc =
    MixedMimeticDiscretization::CatalogInterface::factory( childKey, getDataContext(), childName, this );
  return &this->registerGroup< MixedMimeticDiscretization >( childName, std::move( disc ) );
}

void MixedMimeticDiscretizationManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from MixedMimeticDiscretization here
  for( auto & catalogIter: MixedMimeticDiscretization::getCatalog() )
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}

MixedMimeticDiscretization const &
MixedMimeticDiscretizationManager::getMixedMimeticDiscretization( string const & name ) const
{
  return getGroup< MixedMimeticDiscretization >( name );
}

MixedMimeticDiscretization &
MixedMimeticDiscretizationManager::getMixedMimeticDiscretization( string const & name )
{
  return getGroup< MixedMimeticDiscretization >( name );
}

} // namespace geos
