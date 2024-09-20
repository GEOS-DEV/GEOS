/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file Region.cpp
 */

#include "Region.hpp"


namespace geos
{

using namespace dataRepository;

Region::Region( string const & name,
                Group * const parent )
  : MeshComponentBase( name, parent )
{
  registerWrapper( viewKeyStruct::idString(), &m_id ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Interval region identifier" );

  registerWrapper( viewKeyStruct::pathInRepositoryString(), &m_pathInRepository ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Path of the dataset in the repository" );
}

Region::CatalogInterface::CatalogType & Region::getCatalog()
{
  static Region::CatalogInterface::CatalogType catalog;
  return catalog;
}

REGISTER_CATALOG_ENTRY( MeshComponentBase, Region, string const &, Group * const )

}
