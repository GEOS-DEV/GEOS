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
#include "WellMassRateConstraints.hpp"
#include "WellConstants.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionManager.hpp"


namespace geos
{

using namespace dataRepository;


MassConstraint::MassConstraint( string const & name, Group * const parent )
  : WellConstraintBase( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
}

MassConstraint::~MassConstraint()
{}

void MassConstraint::postInputInitialization()
{
  // Validate value and table options
  WellConstraintBase::postInputInitialization();

}

MassProductionConstraint::MassProductionConstraint( string const & name, Group * const parent )
  : MassConstraint( name, parent )
{
  m_rateSign = -1.0;
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
  registerWrapper( constraintViewStruct::constraintValueKey::constraintValueString(), &m_constraintValue ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Maximum mass injection rate (if useSurfaceConditions: [surface m^3/s]; else [reservoir m^3/s])" );

}


MassProductionConstraint::~MassProductionConstraint()
{}

void MassProductionConstraint::postInputInitialization()
{
  // Validate value and table options
  MassConstraint::postInputInitialization();

}

bool MassProductionConstraint::checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const
{
  return currentConstraint.massRate() < getConstraintValue( currentTime );
}


MassInjectionConstraint::MassInjectionConstraint( string const & name, Group * const parent )
  : MassConstraint( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( constraintViewStruct::constraintValueKey::constraintValueString(), &m_constraintValue ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Maximum mass injection rate (if useSurfaceConditions: [surface m^3/s]; else [reservoir m^3/s])" );
}

MassInjectionConstraint::~MassInjectionConstraint()
{}

void MassInjectionConstraint::postInputInitialization()
{
  // Validate value and table options
  MassConstraint::postInputInitialization();

// Validate the injection stream and temperature
  validateInjectionStream( m_injectionStream, m_injectionTemperature, getConstraintKey(), *this );
}

bool MassInjectionConstraint::checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime )const
{
  return currentConstraint.massRate() >  getConstraintValue( currentTime );
}


} //namespace geos
