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


MassProductionConstraint::MassProductionConstraint( string const & name, Group * const parent )
  : WellConstraintBase( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( constraintViewStruct::constraintValueKey::constraintValueString(), &m_constraintValue ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Maximum mass production rate (if useSurfaceConditions: [surface m^3/s]; else [reservoir m^3/s])" );

  // Field registration
  registerSurfaceConditions( m_useSurfaceConditions, m_surfacePres, m_surfaceTemp, *this );


}


MassProductionConstraint::~MassProductionConstraint()
{}

void MassProductionConstraint::postInputInitialization()
{
  // Validate value and table options
  WellConstraintBase::postInputInitialization();

  // Validate surface conditions
  validateSurfaceConditions( m_useSurfaceConditions, dataRepository::keys::MassInjectionConstraint, *this );
}

MassInjectionConstraint::MassInjectionConstraint( string const & name, Group * const parent )
  : WellConstraintBase( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( constraintViewStruct::constraintValueKey::constraintValueString(), &m_constraintValue ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Maximum mass injection rate (if useSurfaceConditions: [surface m^3/s]; else [reservoir m^3/s])" );

  // Field registration
  registerSurfaceConditions( m_useSurfaceConditions, m_surfacePres, m_surfaceTemp, *this );
  registerInjectionStream( m_injectionStream, m_injectionTemperature, *this );

}


MassInjectionConstraint::~MassInjectionConstraint()
{}

void MassInjectionConstraint::postInputInitialization()
{

  // Validate value and table options
  WellConstraintBase::postInputInitialization();

  // Validate surface conditions
  validateSurfaceConditions( m_useSurfaceConditions, dataRepository::keys::MassInjectionConstraint, *this );
// Validate the injection stream and temperature
  validateInjectionStream( m_injectionStream, m_injectionTemperature, dataRepository::keys::MassProductionConstraint, *this );


}

} //namespace geos
