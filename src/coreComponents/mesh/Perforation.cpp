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

  registerWrapper( viewKeyStruct::perfStatusTableString(), &m_perfStatusTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Table defining perforation state as a function of time. If enterned in Functions section, the name must be the same as the perforation name" );


}


void Perforation::postInputInitialization()
{
  GEOS_ERROR_IF( m_distanceFromHead <= 0,
                 getWrapperDataContext( viewKeyStruct::distanceFromHeadString() ) <<
                 ": distance from well head to perforation cannot be negative." );

  // Setup perforation status function
  FunctionManager & functionManager = FunctionManager::getInstance();
 
  if( !functionManager.hasGroup< TableFunction >( getName() ) )
  {
    TableFunction * tableFunction = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", getName() ) );

    array1d< array1d< real64 > > timeCoord;
    timeCoord.resize( 1 );
    array1d< real64 > values;
    //  Create the time-dependent perforation status table

    if(  m_perfStatusTable[0].size() == 0 )
    {
      real64 alwaysOpen = 1.0;
    timeCoord[0].emplace_back( 0 );
      values.emplace_back( alwaysOpen );
    }
    else
    {
      // If a name is explicitly given, then check that it exists
      GEOS_THROW_IF( m_perfStatusTable[0].size() != m_perfStatusTable[1].size(),
            GEOS_FMT( "Perforation status table missing time or status : {}", getName() ),
              InputError );
      for (std::ptrdiff_t i=0;i<m_perfStatusTable[0].size();i++)
      {
        timeCoord[0].emplace_back( m_perfStatusTable[0][i]);
        values.emplace_back( m_perfStatusTable[1][i]);
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
