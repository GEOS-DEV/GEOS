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
PhaseProductionConstraint::PhaseProductionConstraint( string const & name, Group * const parent )
  : WellConstraintBase( name, parent ),
  m_useSurfaceConditions( 0 ),
  m_surfacePres( 0.0 ),
  m_surfaceTemp( 0.0 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( constraintViewStruct::phaseConstraintKey::targetPhaseRateString(), &m_constraintValue ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Maximum phase production rate,  (if useSurfaceConditions: [surface m^3/s]; else [reservoir m^3/s]) " );

  registerWrapper( constraintViewStruct::phaseConstraintKey::targetPhaseNameString(), &m_targetPhaseName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setDefaultValue( "" ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Name of the target phase" );

  // Field registration
  registerSurfaceConditions( m_useSurfaceConditions, m_surfacePres, m_surfaceTemp, *this );

}

PhaseProductionConstraint::~PhaseProductionConstraint()
{}

void PhaseProductionConstraint::postInputInitialization()
{
  // Validate value and table options
  WellConstraintBase::postInputInitialization();

}

// *** Phase Constraint for Injection Well  ***************************************************************
PhaseInjectionConstraint::PhaseInjectionConstraint( string const & name, Group * const parent )
  : WellConstraintBase( name, parent ),
  m_useSurfaceConditions( 0 ),
  m_surfacePres( 0.0 ),
  m_surfaceTemp( 0.0 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( constraintViewStruct::phaseConstraintKey::targetPhaseRateString(), &m_constraintValue ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Maximum phase injection rate,  (if useSurfaceConditions: [surface m^3/s]; else [reservoir m^3/s]) " );

  registerWrapper( constraintViewStruct::phaseConstraintKey::targetPhaseNameString(), &m_targetPhaseName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setDefaultValue( "" ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Name of the target phase" );

  // Field registration
  registerSurfaceConditions( m_useSurfaceConditions, m_surfacePres, m_surfaceTemp, *this );
  registerInjectionStream( m_injectionStream, m_injectionTemperature, *this );


}


PhaseInjectionConstraint::~PhaseInjectionConstraint()
{}

void PhaseInjectionConstraint::postInputInitialization()
{

  // Validate value and table options
  WellConstraintBase::postInputInitialization();

// Validate the injection stream and temperature
  validateInjectionStream( m_injectionStream, m_injectionTemperature, dataRepository::keys::phaseInjectionConstraint, *this );


}

} //namespace geos
