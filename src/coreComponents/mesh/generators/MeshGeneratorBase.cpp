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

#include "MeshGeneratorBase.hpp"

namespace geosx
{
using namespace dataRepository;

MeshGeneratorBase::MeshGeneratorBase( string const & name, Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
}

Group * MeshGeneratorBase::createChild( string const & childKey, string const & childName )
{
  // Mesh generators generally don't have child XML nodes, must override this method to enable
  GEOSX_THROW( GEOSX_FMT( "Mesh '{}': invalid child XML node '{}' of type {}", getName(), childName, childKey ),
               InputError );
}

MeshGeneratorBase::CatalogInterface::CatalogType & MeshGeneratorBase::getCatalog()
{
  static MeshGeneratorBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

}
