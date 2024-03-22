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

#include "MeshBase.hpp"

namespace geos
{
using namespace dataRepository;

MeshBase::MeshBase( string const & name, Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
}

// Group * MeshBase::createChild( string const & GEOS_UNUSED_PARAM(childKey), string const & GEOS_UNUSED_PARAM(childName) )
// {
//   return nullptr;
// }

// void MeshBase::expandObjectCatalogs()
// {

// }

MeshBase::CatalogInterface::CatalogType & MeshBase::getCatalog()
{
  static MeshBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

}
