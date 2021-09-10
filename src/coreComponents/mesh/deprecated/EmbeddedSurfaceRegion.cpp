/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file EmbeddedSurfaceRegion.cpp
 */

#include "EmbeddedSurfaceRegion.hpp"

#include "EdgeManager.hpp"

namespace geosx
{
using namespace dataRepository;

EmbeddedSurfaceRegion::EmbeddedSurfaceRegion( string const & name, Group * const parent ):
  ElementRegionBase( name, parent )
{
  this->getGroup( viewKeyStruct::elementSubRegions() ).registerGroup< EmbeddedSurfaceSubRegion >( "default" );

  registerWrapper( viewKeyStruct::defaultApertureString(), &m_defaultAperture ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The default aperture of for new embedded surface Elements." );
}

EmbeddedSurfaceRegion::~EmbeddedSurfaceRegion()
{}


void EmbeddedSurfaceRegion::initializePreSubGroups()
{
  this->forElementSubRegions< EmbeddedSurfaceSubRegion >( [&] ( EmbeddedSurfaceSubRegion & subRegion )
  {
    subRegion.getWrapper< array1d< real64 > >( EmbeddedSurfaceSubRegion::viewKeyStruct::elementApertureString() ).
      setApplyDefaultValue( m_defaultAperture );
  } );
}



REGISTER_CATALOG_ENTRY( ObjectManagerBase, EmbeddedSurfaceRegion, string const &, Group * const )

} /* namespace geosx */
