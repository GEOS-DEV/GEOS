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
 * @file WellLiquidRateConstraint.cpp
 */

#include "LogLevelsInfo.hpp"
#include "WellLiquidRateConstraint.hpp"
#include "WellConstants.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionManager.hpp"

namespace geos
{

using namespace dataRepository;


LiquidRateConstraint::LiquidRateConstraint( string const & name, Group * const parent )
  : WellConstraintBase( name, parent )
{
  this->registerWrapper( viewKeyStruct::liquidRateString(), &this->m_constraintValue ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Phase rate,  (if useSurfaceCondSitions: [surface m^3/s]; else [reservoir m^3/s]) " );

  this->registerWrapper( viewKeyStruct::phaseNamesString(), &m_phaseNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "List of fluid phase names defining the liquid" );
}

LiquidRateConstraint::~LiquidRateConstraint()
{}

void LiquidRateConstraint::postInputInitialization()
{
  // Validate table options
  WellConstraintBase::postInputInitialization();

  // check constraint value
  GEOS_THROW_IF( m_constraintValue < 0,
                 getWrapperDataContext( viewKeyStruct::liquidRateString() ) << ": Target value is negative",
                 InputError );

  GEOS_THROW_IF  ((m_constraintValue <= 0.0 && m_constraintScheduleTableName.empty()),
                  getName() << " " << getDataContext() << ": You need to specify a liquid rate constraint. \n" <<
                  "The  rate constraint can be specified using " <<
                  "either " <<  viewKeyStruct::liquidRateString() <<
                  " or " <<  WellConstraintBase::viewKeyStruct::constraintScheduleTableNameString(),
                  InputError );
}


bool LiquidRateConstraint::checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const
{
  real64 const currentValue = currentConstraint.liquidRate();
  real64 const constraintValue = this->getConstraintValue( currentTime );
  return ( LvArray::math::abs( currentValue ) <= LvArray::math::abs( constraintValue ) );
}

} //namespace geos
