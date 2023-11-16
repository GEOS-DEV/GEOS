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
 * @file MaterialPropertiesUpdate.cpp
 */

#include "MaterialPropertiesUpdate.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mainInterface/GeosxState.hpp"

namespace geos
{

using namespace dataRepository;

MaterialPropertiesUpdate::MaterialPropertiesUpdate( const string & name,
                                                    Group * const parent ):
  TaskBase( name, parent )
{
  enableLogLevelInput();
}

MaterialPropertiesUpdate::~MaterialPropertiesUpdate()
{}

void MaterialPropertiesUpdate::postProcessInput()
{}

bool MaterialPropertiesUpdate::execute( real64 const time_n,
                                        real64 const dt,
                                        integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                        integer const GEOS_UNUSED_PARAM( eventCounter ),
                                        real64 const GEOS_UNUSED_PARAM( eventProgress ),
                                        DomainPartition & domain )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  domain.forMeshBodies( [&]( MeshBody & meshBody )
  {
    meshBody.forMeshLevels( [&]( MeshLevel & mesh )
    {
      fsManager.forSubGroups< FieldSpecificationBase >( [&] ( FieldSpecificationBase const & fs )
      {
        if( fs.isPropertyUpdate() )
        {
          fs.apply< dataRepository::Group >( mesh,
                                             [&]( FieldSpecificationBase const & fieldSpec,
                                                  string const &,
                                                  SortedArrayView< localIndex const > const & targetSet,
                                                  Group & targetGroup,
                                                  string const fieldName )
          {
            fieldSpec.applyFieldValue< FieldSpecificationEqual >( targetSet, time_n+dt, targetGroup, fieldName );
          } );
        }
      } );
    } );
  } );

  return false;
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        MaterialPropertiesUpdate,
                        string const &, dataRepository::Group * const )

} /* namespace geos */
