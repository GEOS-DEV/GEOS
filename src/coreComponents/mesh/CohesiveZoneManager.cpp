/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include <vector>

#include "CohesiveZoneManager.hpp"

#include "common/TimingMacros.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsMPM.hpp"
#include "constitutive/ConstitutiveManager.hpp"

#include "mesh/MeshManager.hpp"
#include "schema/schemaUtilities.hpp"

namespace geos
{

using namespace dataRepository;

class SpatialPartition;

// *********************************************************************************************************************
/**
 * @return
 */

CohesiveZoneManager::CohesiveZoneManager( string const & name, Group * const parent ):
  ObjectManagerBase( name, parent )
{
  registerGroup< Group >( CohesiveZoneManager::groupKeyStruct::cohesiveZoneRegionsGroup() );
}

CohesiveZoneManager::~CohesiveZoneManager()
{
  // TODO Auto-generated destructor stub
}

void CohesiveZoneManager::resize( integer_array const & numNodes,
                                  string_array const & regionNames )
{
  localIndex const n_regions = LvArray::integerConversion< localIndex >( regionNames.size() );
  for( localIndex reg = 0; reg < n_regions; ++reg )
  {
    this->getRegion( reg ).resize( numNodes[reg] );
  }
}

Group * CohesiveZoneManager::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( GEOS_FMT( "{}: adding {} {}", getName(), childKey, childName ) );
  Group & cohesiveZoneRegions = this->getGroup( CohesiveZoneManager::groupKeyStruct::cohesiveZoneRegionsGroup() );
  return &cohesiveZoneRegions.registerGroup( childName,
                                             CatalogInterface::factory( childKey, getDataContext(),
                                                                        childName, &cohesiveZoneRegions ) );
}

void CohesiveZoneManager::expandObjectCatalogs()
{
  ObjectManagerBase::CatalogInterface::CatalogType const & catalog = ObjectManagerBase::getCatalog();
  for( ObjectManagerBase::CatalogInterface::CatalogType::const_iterator iter = catalog.begin();
       iter!=catalog.end();
       ++iter )
  {
    string const key = iter->first;
    if( key.find( "CohesiveZoneRegion" ) != string::npos )
    {
      this->createChild( key, key );
    }
  }
}

void CohesiveZoneManager::setSchemaDeviations( xmlWrapper::xmlNode schemaRoot,
                                               xmlWrapper::xmlNode schemaParent,
                                               integer documentationType )
{
  xmlWrapper::xmlNode targetChoiceNode = schemaParent.child( "xsd:choice" );
  if( targetChoiceNode.empty() )
  {
    targetChoiceNode = schemaParent.prepend_child( "xsd:choice" );
    targetChoiceNode.append_attribute( "minOccurs" ) = "0";
    targetChoiceNode.append_attribute( "maxOccurs" ) = "unbounded";
  }

  std::set< string > names;
  this->forCohesiveZoneRegions( [&]( CohesiveZoneRegionBase & cohesiveZoneRegion )
  {
    names.insert( cohesiveZoneRegion.getName() );
  } );

  for( string const & name: names )
  {
    schemaUtilities::SchemaConstruction( getRegion( name ), schemaRoot, targetChoiceNode, documentationType );
  }
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, CohesiveZoneManager, string const &, Group * const )
}
