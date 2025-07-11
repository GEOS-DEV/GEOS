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
 * @file WellVolumeRateConstraints.cpp
 */

#include "LogLevelsInfo.hpp"
#include "WellVolumeRateConstraints.hpp"
#include "WellConstants.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionManager.hpp"


namespace geos
{

using namespace dataRepository;


VolumeProductionConstraint::VolumeProductionConstraint( string const & name, Group * const parent )
  : WellConstraintBase( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( constraintViewStruct::constraintValueKey::constraintValueString(), &m_constraintValue ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Maximum volumetric production rate (if useSurfaceConditions: [surface m^3/s]; else [reservoir m^3/s])" );

  registerWrapper( constraintViewStruct::surfaceConditionsKey::useSurfaceConditionsString(), &m_useSurfaceConditions ).
    setDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to specify whether rates are checked at surface or reservoir conditions.\n"
                    "Equal to 1 for surface conditions, and to 0 for reservoir conditions" );

  registerWrapper( constraintViewStruct::surfaceConditionsKey::surfacePressureString(), &m_surfacePres ).
    setDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Surface pressure used to compute volumetric rates when surface conditions are used [Pa]" );

  registerWrapper( constraintViewStruct::surfaceConditionsKey::surfaceTemperatureString(), &m_surfaceTemp ).
    setDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Surface temperature used to compute volumetric rates when surface conditions are used [K]" );


}


VolumeProductionConstraint::~VolumeProductionConstraint()
{}

void VolumeProductionConstraint::postInputInitialization()
{
  WellConstraintBase::postInputInitialization();

}

VolumeInjectionConstraint::VolumeInjectionConstraint( string const & name, Group * const parent )
  : WellConstraintBase( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( constraintViewStruct::constraintValueKey::constraintValueString(), &m_constraintValue ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Maximum volumetric injection rate (if useSurfaceConditions: [surface m^3/s]; else [reservoir m^3/s])" );


  // Field registration
  registerSurfaceConditions( m_useSurfaceConditions, m_surfacePres, m_surfaceTemp, *this );
  registerInjectionStream( m_injectionStream, m_injectionTemperature, *this );
}


VolumeInjectionConstraint::~VolumeInjectionConstraint()
{}

void VolumeInjectionConstraint::postInputInitialization()
{

// Validate the injection stream and temperature
  validateInjectionStream( m_injectionStream, m_injectionTemperature, dataRepository::keys::volumeProductionConstraint, *this );


}

} //namespace geos
