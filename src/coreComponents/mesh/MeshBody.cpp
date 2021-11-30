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



MeshLevel & MeshBody::createMeshLevel( localIndex const newLevel )
{
  string const name = "Level" + std::to_string(newLevel);
  return this->registerGroup< MeshLevel >( name );
}

MeshLevel & MeshBody::createMeshLevel( localIndex const newLevel, localIndex const sourceLevel )
{
  string const sourceName = "Level" + std::to_string(sourceLevel);
  string const newName = "Level" + std::to_string(newLevel);
  MeshLevel const & sourceMeshLevel = this->getGroupByPath<MeshLevel>( sourceName );

//  MeshLevel newMeshLevel( newName, this, sourceMeshLevel, 1 );


  return this->registerGroup( newName, std::make_unique<MeshLevel>( newName, this, sourceMeshLevel, 1  ) );
}


void MeshBody::setGlobalLengthScale( real64 scale )
{
  m_globalLengthScale = scale;
}

} /* namespace geosx */
