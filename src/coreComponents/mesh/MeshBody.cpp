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

/**
 * @file MeshBody.cpp
 */

#include "MeshBody.hpp"
#include "MeshLevel.hpp"

namespace geosx
{
using namespace dataRepository;

MeshBody::MeshBody( string const & name,
                    Group * const parent ):
  Group( name, parent ),
  m_globalLengthScale( 0 )
{
  registerWrapper< integer >( viewKeys.meshLevels );
}

MeshBody::~MeshBody()
{
  // TODO Auto-generated destructor stub
}



MeshLevel & MeshBody::createMeshLevel( localIndex const GEOSX_UNUSED_PARAM( newLevel ) )
{
  return this->registerGroup< MeshLevel >( "Level0" );
}

void MeshBody::setGlobalLengthScale( real64 scale )
{
  m_globalLengthScale = scale;
}

} /* namespace geosx */
