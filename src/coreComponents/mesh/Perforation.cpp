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

/*
 * @file Perforation.cpp
 */

#include "Perforation.hpp"
#include "dataRepository/InputFlags.hpp"

namespace geos
{

using namespace dataRepository;

Perforation::Perforation( string const & name, Group * const parent )
  : Group( name, parent ),
  m_distanceFromHead( 0 ),
  m_wellTransmissibility( 0 ),
  m_wellSkinFactor( 0 ),
  m_perfStatusTable( "PerfStatusTable", parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::distanceFromHeadString(), &m_distanceFromHead ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Linear distance from well head to the perforation" );

  registerWrapper( viewKeyStruct::wellTransmissibilityString(), &m_wellTransmissibility ).
    setApplyDefaultValue( -1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Perforation transmissibility" );

  registerWrapper( viewKeyStruct::wellSkinFactorString(), &m_wellSkinFactor ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Perforation skin factor" );

  registerWrapper( viewKeyStruct::targetRegionString(), &m_targetRegionName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Target region to connect the perforation" );
#if 0
  registerWrapper( "PerfStatusTable", &m_perfStatusTable ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The perforation table defining perforation state as a function of time. \n"
                    "If the status function evaluates to a positive value at the current time, the perforation will be open otherwise the perforation will be shut." );
#else
  registerWrapper( viewKeyStruct::perfStatusTableString(), &m_perfStatus ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "The perforation table defining perforation state as a function of time." );
#endif
}


void Perforation::postInputInitialization()
{
  GEOS_ERROR_IF( m_distanceFromHead <= 0,
                 getWrapperDataContext( viewKeyStruct::distanceFromHeadString() ) <<
                 ": distance from well head to perforation cannot be negative." );

  // 12) Create the time-dependent perforation status table
  if( m_perfStatusTable.numDimensions() == 0 )
  {
    real64 constantValue = 1.0;
    array1d< array1d< real64 > > timeCoord;
    timeCoord.resize( 1 );
    timeCoord[0].emplace_back( 0 );
    array1d< real64 > constantValueArray;
    constantValueArray.emplace_back( constantValue );
    m_perfStatusTable.setTableCoordinates( timeCoord, { units::Time } );
    m_perfStatusTable.setTableValues( constantValueArray );
    m_perfStatusTable.setInterpolationMethod( TableFunction::InterpolationType::Lower );
  }
}

Perforation::~Perforation()
{}

} //namespace geos
