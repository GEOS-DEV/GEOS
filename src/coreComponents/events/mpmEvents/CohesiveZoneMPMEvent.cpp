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
 * @file CohesiveZoneMPMEvent.cpp
 */

#include "CohesiveZoneMPMEvent.hpp"

#include "constitutive/cohesiveZone/CohesiveZoneBase.hpp"

namespace geos
{

using namespace dataRepository;

CohesiveZoneMPMEvent::CohesiveZoneMPMEvent( const string & name,
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

CohesiveZoneMPMEvent::~CohesiveZoneMPMEvent()
{}

// Group * CohesiveZoneMPMEvent::createChild( string const & childKey, string const & childName )
// {
//   GEOS_LOG_RANK_0( GEOS_FMT( "{}: adding {} {}", getName(), childKey, childName ) );
//   return &registerGroup( childName,
//                          CatalogInterface::factory( childKey, getDataContext(),
//                                                     childName, this ) );
// }

void CohesiveZoneMPMEvent::postInputInitialization()
{
  MPMEventBase::postInputInitialization();

  GEOS_ERROR_IF( !( ( m_regionNames.size() == m_constitutiveModelNames.size() ) & ( m_constitutiveModelNames.size() == static_cast< long unsigned int >( m_czTags.size() ) ) ), "Region names, constitutive model names, and cz tags must be the same length" );
  GEOS_ERROR_IF( m_regionNames.size() == 0, "Region names, constitutive model names, and cz tags must not be empty" );

  GEOS_ERROR_IF( !( m_czVolumeNormalization == 0 || m_czVolumeNormalization == 1 ), "czVolumeNormalization can only be 0 or 1" );
  GEOS_ERROR_IF( !( m_computeNormalsAndPositions == 0 || m_computeNormalsAndPositions == 1 ), "computeNormalsAndPositions can only be 0 or 1" );
}

REGISTER_CATALOG_ENTRY( MPMEventBase, CohesiveZoneMPMEvent, string const &, Group * const )

} /* namespace geos */
