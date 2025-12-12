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
 * @file WellMassRateConstraints.cpp
 */

#include "LogLevelsInfo.hpp"
#include "WellMassRateConstraint.hpp"
#include "WellConstants.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionManager.hpp"

namespace geos
{

using namespace dataRepository;

MassRateConstraint::MassRateConstraint( string const & name, Group * const parent )
  : WellConstraintBase( name, parent )
{
  this->setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  this->registerWrapper( viewKeyStruct::massRateString(), &this->m_constraintValue ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Maximum mass rate (kg/s)" );
}

MassRateConstraint::~MassRateConstraint()
{}


void MassRateConstraint::postInputInitialization()
{
  // Validate table options
  WellConstraintBase::postInputInitialization();

  // check constraint value
  GEOS_THROW_IF( m_constraintValue < 0,
                 getWrapperDataContext( viewKeyStruct::massRateString() ) << ": Target value is negative",
                 InputError );

  GEOS_THROW_IF  ((m_constraintValue <= 0.0 && m_constraintScheduleTableName.empty()),
                  getName() << " " << getDataContext() << ": You need to specify a mass rate constraint. \n" <<
                  "The  rate constraint can be specified using " <<
                  "either " <<  viewKeyStruct::massRateString() <<
                  " or " <<  WellConstraintBase::viewKeyStruct::constraintScheduleTableNameString(),
                  InputError );
}


bool MassRateConstraint::checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime )const
{
  // isViolated is defined as a static method on the specific WellConstraintType (Injection/Production)
  // Evaluate violation according to the sign set for injectors/producers
  real64 const currentValue = currentConstraint.massRate();
  real64 const constraintValue = this->getConstraintValue( currentTime );
  return ( LvArray::math::abs( currentValue ) > LvArray::math::abs( constraintValue ) );

}
} //namespace geos
