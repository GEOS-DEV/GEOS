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
 * @file MixedVEMManager.cpp
 */

#include "mixedVEM/MixedVEMManager.hpp"

#include "mixedVEM/MixedVEMDiscretization.hpp"

namespace geos
{

using namespace dataRepository;

MixedVEMManager::MixedVEMManager( string const & name, Group * const parent )
  : Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
}

MixedVEMManager::~MixedVEMManager()
{}

Group * MixedVEMManager::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( GEOS_FMT( "{}: adding {} {}", getName(), childKey, childName ) );
  GEOS_ERROR_IF_NE_MSG( childKey, MixedVEMDiscretization::catalogName(),
                        getDataContext() << ": unknown child " << childKey );

  std::unique_ptr< MixedVEMDiscretization > discretization =
    std::make_unique< MixedVEMDiscretization >( childName, this );

  return &this->registerGroup< MixedVEMDiscretization >( childName, std::move( discretization ) );
}

void MixedVEMManager::expandObjectCatalogs()
{
  createChild( MixedVEMDiscretization::catalogName(), MixedVEMDiscretization::catalogName() );
}

MixedVEMDiscretization const & MixedVEMManager::getDiscretization( string const & name ) const
{
  return getGroup< MixedVEMDiscretization >( name );
}

MixedVEMDiscretization & MixedVEMManager::getDiscretization( string const & name )
{
  return getGroup< MixedVEMDiscretization >( name );
}

} // namespace geos
