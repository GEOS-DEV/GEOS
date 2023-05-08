/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MultiscaleMeshObjectManager.cpp
 */

#include "MeshObjectManager.hpp"

namespace geos
{

using namespace dataRepository;

namespace multiscale
{

MeshObjectManager::MeshObjectManager( string const & name,
                                      dataRepository::Group * const parent )
  : ObjectManagerBase( name, parent ),
  m_numOwnedObjects( 0 )
{
  registerWrapper( viewKeyStruct::dualObjectString(), &m_toDualRelation );

  excludeWrappersFromPacking( { viewKeyStruct::dualObjectString() } );
}

void MeshObjectManager::setNumOwnedObjects( localIndex const n )
{
  GEOS_ERROR_IF_GT_MSG( n, size(), "Number of owned objects cannot exceed the total size. This is an internal logic error." );
  m_numOwnedObjects = n;
}

} // namespace multiscale
} // namespace geos
