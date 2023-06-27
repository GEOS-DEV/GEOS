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

#include "WellGeneratorBase.hpp"
#include "mesh/Perforation.hpp"

namespace geos
{
using namespace dataRepository;

WellGeneratorBase::WellGeneratorBase( string const & name, Group * const parent ):
  Group( name, parent ),
  m_numPerforations( 0 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
}

Group * WellGeneratorBase::createChild( string const & childKey, string const & childName )
{
  if( childKey == viewKeyStruct::perforationString() )
  {
    ++m_numPerforations;

    // keep track of the perforations that have been added
    m_perforationList.emplace_back( childName );

    return &registerGroup< Perforation >( childName );
  }
  else
  {
    GEOS_THROW( "Unrecognized node: " << childKey, InputError );
  }
  return nullptr;
}


void WellGeneratorBase::expandObjectCatalogs()
{
  createChild( viewKeyStruct::perforationString(), viewKeyStruct::perforationString() );
}

WellGeneratorBase::CatalogInterface::CatalogType & WellGeneratorBase::getCatalog()
{
  static WellGeneratorBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

}
