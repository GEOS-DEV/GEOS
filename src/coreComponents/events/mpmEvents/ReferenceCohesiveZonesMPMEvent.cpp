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
  m_normalsAndPositionsMethod( mpm::NormalsAndPositionsMethodOption::LogisticRegression ),
  m_czSurfaceDisplacementUpdate( mpm::CohesiveSurfaceDisplacementUpdateOption::Nodal )
{
  registerWrapper( viewKeyStruct::regionNamesString(), &m_regionNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Region names for cohesive zones");

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

  registerWrapper( viewKeyStruct::czSurfaceDisplacementUpdateString(), &m_czSurfaceDisplacementUpdate ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_czSurfaceDisplacementUpdate ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Cohesive surface displacement update method. Nodal preserves the bugfix branch nodal update; TypeA uses the stored surface-position vector; TypeB uses the CPDI particle-face vector. Options are:\n* " +
                    EnumStrings< mpm::CohesiveSurfaceDisplacementUpdateOption >::concat( "\n* " ) );
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

  GEOS_ERROR_IF( !( m_czVolumeNormalization == 0 || m_czVolumeNormalization == 1 ),
                 "czVolumeNormalization can only be 0 or 1." );
}

REGISTER_CATALOG_ENTRY( MPMEventBase, ReferenceCohesiveZonesMPMEvent, string const &, Group * const )

} /* namespace geos */
