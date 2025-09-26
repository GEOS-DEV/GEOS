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
#include "functions/FunctionManager.hpp"

namespace geos
{

using namespace dataRepository;

Perforation::Perforation( string const & name, Group * const parent )
  : Group( name, parent ),
  m_distanceFromHead( 0 ),
  m_wellTransmissibility( 0 ),
  m_wellSkinFactor( 0 )
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

  registerWrapper( viewKeyStruct::perfStatusTableNameString(), &m_perfStatusTableName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription(
    "Name of table function defining perforation status as a function of time:\n"
    "- The first row sets the perforation time coordinates,\n"
    "- The following row is the perforation status ( 0.0 for closed, 1.0 for opened) for each time coordinate.\n"
    "If none specified, a table function with be created internally and assigned the name formed by well and perforation names, separated by '_'.\n"
    "A table that opens the perforation at t=100, and closes it at t=200 could be `{{ 0, 100, 200 },{ 0.0, 1.0, 0.0 }}`" );
  registerWrapper( viewKeyStruct::perfStatusTableString(), &m_perfStatusTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Table defining perforation status as a function of time. The name must correspond to a TableFunction in the Functions section." );


}


void Perforation::postInputInitialization()
{
  GEOS_ERROR_IF( m_distanceFromHead <= 0,
                 getWrapperDataContext( viewKeyStruct::distanceFromHeadString() ) <<
                 ": distance from well head to perforation cannot be negative.",
                 getWrapperDataContext( viewKeyStruct::distanceFromHeadString() ) );

  // Setup perforation status function
  FunctionManager & functionManager = FunctionManager::getInstance();

  if( m_perfStatusTableName.empty() )
  {
    // No table function provided, defaults name to perforation name
    m_perfStatusTableName = getParent().getName()+"_"+getName();
  }
  else
  {
    // Table name provided as input, check that it exists
    GEOS_THROW_IF( !functionManager.hasGroup< TableFunction >( m_perfStatusTableName ),
                   GEOS_FMT( "{}: missing perforation status table `{}`",
                             getDataContext(), m_perfStatusTableName ),
                   InputError );
  }

  if( !functionManager.hasGroup< TableFunction >( m_perfStatusTableName ) )
  {
    //  Create the time-dependent perforation status table
    TableFunction * tableFunction = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", m_perfStatusTableName ) );

    array1d< array1d< real64 > > timeCoord;
    timeCoord.resize( 1 );
    array1d< real64 > values;

    if( m_perfStatusTable.size( 0 ) == 0 )
    {
      // No table supplied set all perfs to open
      real64 alwaysOpen = 1.0;
      timeCoord[0].emplace_back( 0 );
      values.emplace_back( alwaysOpen );
    }
    else
    {
      // User supplied table in Perforation section

      GEOS_THROW_IF( m_perfStatusTable[0].size() != m_perfStatusTable[1].size(),
                     GEOS_FMT( "{}: Perforation status table `{}` missing time or status.", m_perfStatusTableName, getDataContext() ),
                     InputError );

      for( std::ptrdiff_t i=0; i<m_perfStatusTable[0].size(); i++ )
      {
        timeCoord[0].emplace_back( m_perfStatusTable[0][i] );
        GEOS_THROW_IF( ( !isZero( m_perfStatusTable[1][i] ) && !isZero( 1 - m_perfStatusTable[1][i] ) ),
                       GEOS_FMT( "{}: Perforation status value must be 0 or 1.", getDataContext() ),
                       InputError );
        values.emplace_back( m_perfStatusTable[1][i] );
      }
    }
    tableFunction->setTableCoordinates( timeCoord, { units::Time } );
    tableFunction->setTableValues( values );
    tableFunction->setInterpolationMethod( TableFunction::InterpolationType::Lower );
  }
}

Perforation::~Perforation()
{}

} //namespace geos
