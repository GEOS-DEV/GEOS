/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * @file FiniteVolumeManager.cpp
 *
 */

#include "FiniteVolumeManager.hpp"

#include "finiteVolume/FluxApproximationBase.hpp"
#include "finiteVolume/HybridMimeticDiscretization.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{

using namespace dataRepository;


FiniteVolumeManager::FiniteVolumeManager( std::string const & name, Group * const parent )
  : Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
}

FiniteVolumeManager::~FiniteVolumeManager()
{}

Group * FiniteVolumeManager::createChild( std::string const & childKey, std::string const & childName )
{
  if( childKey == HybridMimeticDiscretization::catalogName() )
  {
    std::unique_ptr< HybridMimeticDiscretization > hm = std::make_unique< HybridMimeticDiscretization >( childName, this );
    return this->registerGroup< HybridMimeticDiscretization >( childName, std::move( hm ) );
  }
  else
  {
    std::unique_ptr< FluxApproximationBase > approx = FluxApproximationBase::CatalogInterface::factory( childKey, childName, this );
    return this->registerGroup< FluxApproximationBase >( childName, std::move( approx ));
  }
}


void FiniteVolumeManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from FluxApproximationBase here
  for( auto & catalogIter: FluxApproximationBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
  // Then do the same thing for the HybridMimeticDiscretization
  for( auto & catalogIter: HybridMimeticDiscretization::getCatalog())
  {
    std::string const childName = catalogIter.first;
    std::unique_ptr< HybridMimeticDiscretization > hm = std::make_unique< HybridMimeticDiscretization >( childName, this );
    this->registerGroup< HybridMimeticDiscretization >( childName, std::move( hm ) );
  }
}


FluxApproximationBase const & FiniteVolumeManager::getFluxApproximation( std::string const & name ) const
{
  return getGroupReference< FluxApproximationBase >( name );
}

FluxApproximationBase & FiniteVolumeManager::getFluxApproximation( std::string const & name )
{
  return getGroupReference< FluxApproximationBase >( name );
}

HybridMimeticDiscretization const & FiniteVolumeManager::getHybridMimeticDiscretization( std::string const & name ) const
{
  return getGroupReference< HybridMimeticDiscretization >( name );
}

HybridMimeticDiscretization & FiniteVolumeManager::getHybridMimeticDiscretization( std::string const & name )
{
  return getGroupReference< HybridMimeticDiscretization >( name );
}


} // namespace geosx
