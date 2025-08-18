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
 * @file WellConstraintBase.cpp
 */

#include "LogLevelsInfo.hpp"
#include "WellConstraintsBase.hpp"
#include "WellConstants.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionManager.hpp"


namespace geos
{

using namespace dataRepository;

namespace
{


#if 0  
/// Utility function to create a one-value table internally when not provided by the user
TableFunction * createConstraintScheduleTable( string const & tableName,
                                               real64 const & constantValue )
{
  array1d< array1d< real64 > > timeCoord;
  timeCoord.resize( 1 );
  timeCoord[0].emplace_back( 0 );
  array1d< real64 > constantValueArray;
  constantValueArray.emplace_back( constantValue );

  FunctionManager & functionManager = FunctionManager::getInstance();
  TableFunction * table = dynamicCast< TableFunction * >( functionManager.createChild( TableFunction::catalogName(), tableName ));
  table->setTableCoordinates( timeCoord, { units::Time } );
  table->setTableValues( constantValueArray );
  table->setInterpolationMethod( TableFunction::InterpolationType::Lower );
  return table;
}
#endif


}

WellConstraintBase::WellConstraintBase( string const & name, Group * const parent )
  : Group( name, parent ),
  m_isConstraintActive( true ),
  m_useScheduleTable( false ),
  m_constraintValue( 0 ),
  m_constraintScheduleTable( nullptr ),
  m_rateSign( 1.0 ) // Default to positive rate sign for injection, set to -1.0 for production wells

{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::constraintScheduleTableNameString(), &m_constraintScheduleTableName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the well constraint schedule table when the constraint value  is a time dependent function. \n" );

}


WellConstraintBase::~WellConstraintBase()
{}


void WellConstraintBase::postInputInitialization()
{

  // check constraint value
  GEOS_THROW_IF( m_constraintValue < 0,
                 getWrapperDataContext( constraintViewStruct::constraintValueKey::constraintValueString() ) << ": Target value is negative",
                 InputError );


  GEOS_THROW_IF( ((m_constraintValue > 0.0 && !m_constraintScheduleTableName.empty())),
                 getConstraintKey() << " " << getDataContext() << ": You have provided redundant information for well constraint value ." <<
                 " A constraint value and table of constraint values cannot be specified together",
                 InputError );

  //  Create time-dependent constraint table
  if( !m_constraintScheduleTableName.empty() )
  {
    FunctionManager & functionManager = FunctionManager::getInstance();
    m_constraintScheduleTable = &(functionManager.getGroup< TableFunction const >( m_constraintScheduleTableName ));

    GEOS_THROW_IF( m_constraintScheduleTable->getInterpolationMethod() != TableFunction::InterpolationType::Lower,
                   getConstraintKey() << " " << getDataContext() << ": The interpolation method for the schedule table "
                                      << m_constraintScheduleTable->getName() << " should be TableFunction::InterpolationType::Lower",
                   InputError );
  }


  GEOS_THROW_IF  ((m_constraintValue <= 0.0 && m_constraintScheduleTableName.empty()),
                  getConstraintKey() << " " << getDataContext() << ": You need to specify a volume rate constraint. \n" <<
                  "The  rate constraint can be specified using " <<
                  "either " << constraintViewStruct::constraintValueKey::constraintValueString() <<
                  " or " << constraintViewStruct::constraintValueKey::constraintScheduleTableNameString(),
                  InputError );

}

void WellConstraintBase::setNextDtFromTables( real64 const currentTime, real64 & nextDt )
{
  setNextDtFromTable( m_constraintScheduleTable, currentTime, nextDt );
}

void WellConstraintBase::setNextDtFromTable( TableFunction const * table, real64 const currentTime, real64 & nextDt )
{
  if( table )
  {
    // small epsilon to make sure we land on the other side of table interval and pick up the right rate
    real64 const eps = 1e-6;
    real64 const dtLimit = (table->getCoord( &currentTime, 0, TableFunction::InterpolationType::Upper ) - currentTime) * ( 1.0 + eps );
    if( dtLimit > eps && dtLimit < nextDt )
    {
      nextDt = dtLimit;
    }
  }
}

} //namespace geos
