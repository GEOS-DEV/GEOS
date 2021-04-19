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
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

using namespace dataRepository;


FiniteVolumeManager::FiniteVolumeManager( string const & name, Group * const parent )
  : Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
}

FiniteVolumeManager::~FiniteVolumeManager()
{}

Group * FiniteVolumeManager::createChild( string const & childKey, string const & childName )
{
  if( childKey == HybridMimeticDiscretization::catalogName() )
  {
    std::unique_ptr< HybridMimeticDiscretization > hm = std::make_unique< HybridMimeticDiscretization >( childName, this );
    return &this->registerGroup< HybridMimeticDiscretization >( childName, std::move( hm ) );
  }
  else
  {
    std::unique_ptr< FluxApproximationBase > approx = FluxApproximationBase::CatalogInterface::factory( childKey, childName, this );
    return &this->registerGroup< FluxApproximationBase >( childName, std::move( approx ));
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
    string const childName = catalogIter.first;
    std::unique_ptr< HybridMimeticDiscretization > hm = std::make_unique< HybridMimeticDiscretization >( childName, this );
    this->registerGroup< HybridMimeticDiscretization >( childName, std::move( hm ) );
  }
}


FluxApproximationBase const & FiniteVolumeManager::getFluxApproximation( string const & name ) const
{
  return getGroup< FluxApproximationBase >( name );
}

FluxApproximationBase & FiniteVolumeManager::getFluxApproximation( string const & name )
{
  return getGroup< FluxApproximationBase >( name );
}

HybridMimeticDiscretization const & FiniteVolumeManager::getHybridMimeticDiscretization( string const & name ) const
{
  return getGroup< HybridMimeticDiscretization >( name );
}

HybridMimeticDiscretization & FiniteVolumeManager::getHybridMimeticDiscretization( string const & name )
{
  return getGroup< HybridMimeticDiscretization >( name );
}


} // namespace geosx
