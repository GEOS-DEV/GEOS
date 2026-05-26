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
 * @file ReferenceCohesiveZonesMPMEvent.cpp
 */

#include "ReferenceCohesiveZonesMPMEvent.hpp"

#include "constitutive/cohesiveZone/CohesiveZoneBase.hpp"

namespace geos
{

using namespace dataRepository;

ReferenceCohesiveZonesMPMEvent::ReferenceCohesiveZonesMPMEvent( const string & name,
                                                              Group * const parent ):
  MPMEventBase( name, parent ),
  m_czVolumeNormalization( 1 ),
  m_computeNormalsAndPositions( 0 ),
  m_normalsAndPositionsMethod( mpm::NormalsAndPositionsMethodOption::LogisticRegression )
{
  registerWrapper( viewKeyStruct::regionNamesString(), &m_regionNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Region names for cohesive zones");

  registerWrapper( viewKeyStruct::constitutiveModelsString(), &m_constitutiveModelNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Constitutive model names for cohesive zones");

  registerWrapper( viewKeyStruct::czTagsString(), &m_czTags ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Tag IDs for cohesive zones");

  registerWrapper( viewKeyStruct::czVolumeNormalizationString(), &m_czVolumeNormalization ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_czVolumeNormalization ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Flag to perform volume normalization of cohesive nodal area for partially filled grid cells" );

  registerWrapper( viewKeyStruct::computeNormalsAndPositionsString(), &m_computeNormalsAndPositions ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_computeNormalsAndPositions ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Flag to automicatlly compute particle surface normals and positions" );

  registerWrapper( viewKeyStruct::normalsAndPositionsMethodString(), &m_normalsAndPositionsMethod ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_normalsAndPositionsMethod ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Method for computing particle surface normals and positions" );
}

ReferenceCohesiveZonesMPMEvent::~ReferenceCohesiveZonesMPMEvent()
{}

// Group * CohesiveZoneMPMEvent::createChild( string const & childKey, string const & childName )
// {
//   GEOS_LOG_RANK_0( GEOS_FMT( "{}: adding {} {}", getName(), childKey, childName ) );
//   return &registerGroup( childName,
//                          CatalogInterface::factory( childKey, getDataContext(),
//                                                     childName, this ) );
// }

void ReferenceCohesiveZonesMPMEvent::postInputInitialization()
{
  MPMEventBase::postInputInitialization();

  GEOS_ERROR_IF( m_regionNames.size() == 0,
                 "ReferenceCohesiveZone event must specify at least one region name." );

  bool const hasModels = m_constitutiveModelNames.size() > 0;
  bool const hasTags = m_czTags.size() > 0;

  GEOS_ERROR_IF( hasModels && hasTags,
                 "ReferenceCohesiveZones event must specify both constitutiveModels and czTags" );

  GEOS_ERROR_IF( m_regionNames.size() != m_constitutiveModelNames.size() ||
                 m_regionNames.size() != static_cast< size_t >( m_czTags.size() ),
                 "ReferenceCohesiveZones event regionNames, constitutiveModels, and czTags "
                 "must have the same length." );

  GEOS_ERROR_IF( !( m_czVolumeNormalization == 0 || m_czVolumeNormalization == 1 ),
                 "czVolumeNormalization can only be 0 or 1." );
}

REGISTER_CATALOG_ENTRY( MPMEventBase, ReferenceCohesiveZonesMPMEvent, string const &, Group * const )

} /* namespace geos */
