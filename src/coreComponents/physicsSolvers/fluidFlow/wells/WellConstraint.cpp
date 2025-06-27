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
 * @file WellConstraint.cpp
 */

#include "LogLevelsInfo.hpp"
#include "WellConstraint.hpp"
#include "WellConstants.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionManager.hpp"


namespace geos
{

using namespace dataRepository;

namespace
{

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
}

WellConstraint::WellConstraint( string const & name, Group * const parent )
  : Group( name, parent ),
  m_constraintScheduleTable( nullptr )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );


  registerWrapper( viewKeyStruct::constraintScheduleTableNameString(), &m_constraintScheduleTableName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the well constraint schedule table when the constraint value  is a time dependent function. \n" );

}


WellConstraint::~WellConstraint()
{}



void WellConstraint::postInputInitialization()
{}



void WellConstraint::setNextDtFromTables( real64 const currentTime, real64 & nextDt )
{
  setNextDtFromTable( m_constraintScheduleTable, currentTime, nextDt );
}

void WellConstraint::setNextDtFromTable( TableFunction const * table, real64 const currentTime, real64 & nextDt )
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


WellBHPConstraint::WellBHPConstraint( string const & name, Group * const parent )
  : WellConstraint( name, parent ),
  m_targetBHP( 0.0 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::targetBHPString(), &m_targetBHP ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Target bottom-hole pressure [Pa]" );

}


WellBHPConstraint::~WellBHPConstraint()
{}

void WellBHPConstraint::postInputInitialization()
{

  // check target BHP
  GEOS_THROW_IF( m_targetBHP < 0,
                 getWrapperDataContext( viewKeyStruct::targetBHPString() ) <<
                 ": Target bottom-hole pressure is negative",
                 InputError );


  if( m_targetBHP <= 0.0 && m_constraintScheduleTableName.empty() )
  {
    m_targetBHP = isProductionConstraint() ? WellConstants::defaultProducerBHP : WellConstants::defaultInjectorBHP;
    GEOS_LOG_LEVEL_RANK_0( logInfo::WellControl,
                           GEOS_FMT( "WellBHPConstraint {}: Setting {}  to default value {}", getDataContext(), viewKeyStruct::targetBHPString(), m_targetBHP ));

  }

  //  Create time-dependent BHP table
  if( m_constraintScheduleTableName.empty() )
  {
    m_constraintScheduleTableName = getName()+"_ConstantBHP_table";
    m_constraintScheduleTable = createConstraintScheduleTable( m_constraintScheduleTableName, m_targetBHP );
  }
  else
  {
    FunctionManager & functionManager = FunctionManager::getInstance();
    m_constraintScheduleTable = &(functionManager.getGroup< TableFunction const >( m_constraintScheduleTableName ));

    GEOS_THROW_IF( m_constraintScheduleTable->getInterpolationMethod() != TableFunction::InterpolationType::Lower,
                   "WellBHPConstraint " << getDataContext() << ": The interpolation method for the schedule table "
                                        << m_constraintScheduleTable->getName() << " should be TableFunction::InterpolationType::Lower",
                   InputError );
  }


}


} //namespace geos
