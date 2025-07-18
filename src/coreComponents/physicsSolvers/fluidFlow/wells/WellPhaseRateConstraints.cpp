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
#include "WellPhaseRateConstraints.hpp"
#include "WellConstants.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionManager.hpp"


namespace geos
{

using namespace dataRepository;

// *** Phase Constraint for Production Well  ***************************************************************
PhaseConstraint::PhaseConstraint( string const & name, Group * const parent )
  : WellConstraintBase( name, parent ),
  m_phaseIndex( -1 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::phaseRateString(), &m_constraintValue ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Maximum phase production rate,  (if useSurfaceConditions: [surface m^3/s]; else [reservoir m^3/s]) " );

  registerWrapper( viewKeyStruct::phaseNameString(), &m_phaseName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setDefaultValue( "" ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Name of the target phase" );

}

PhaseConstraint::~PhaseConstraint()
{}

void PhaseConstraint::postInputInitialization()
{
  // Validate value and table options
  WellConstraintBase::postInputInitialization();
}


// *** Phase Constraint for Production Well  ***************************************************************
PhaseProductionConstraint::PhaseProductionConstraint( string const & name, Group * const parent )
  : PhaseConstraint( name, parent )
{}

PhaseProductionConstraint::~PhaseProductionConstraint()
{}

void PhaseProductionConstraint::postInputInitialization()
{
  // Validate value and table options
  PhaseConstraint::postInputInitialization();

}

bool PhaseProductionConstraint::checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const
{
  return -1.0*currentConstraint.phaseVolumeRates()[m_phaseIndex] > getConstraintValue( currentTime );
}



// *** Phase Constraint for Injection Well  ***************************************************************
PhaseInjectionConstraint::PhaseInjectionConstraint( string const & name, Group * const parent )
  : PhaseConstraint( name, parent )
{
  registerWrapper( constraintViewStruct::injectionStreamKey::injectionStreamString(), &m_injectionStream ).
    setDefaultValue( -1 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Global component densities of the injection stream [moles/m^3 or kg/m^3]" );

  registerWrapper( constraintViewStruct::injectionStreamKey::injectionTemperatureString(), &m_injectionTemperature ).
    setDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Temperature of the injection stream [K]" );
}


PhaseInjectionConstraint::~PhaseInjectionConstraint()
{}

void PhaseInjectionConstraint::postInputInitialization()
{

  // Validate value and table options
  PhaseConstraint::postInputInitialization();

// Validate the injection stream and temperature
  validateInjectionStream( m_injectionStream, m_injectionTemperature, dataRepository::keys::phaseInjectionConstraint, *this );

}

bool PhaseInjectionConstraint::checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const
{
  return currentConstraint.phaseVolumeRates()[m_phaseIndex] > getConstraintValue( currentTime );
}



} //namespace geos
