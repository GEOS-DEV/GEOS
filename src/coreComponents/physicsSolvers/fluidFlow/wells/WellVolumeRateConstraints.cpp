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

VolumeConstraint::VolumeConstraint( string const & name, Group * const parent )
  : WellConstraintBase( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
  registerWrapper( viewKeyStruct::volumeRateString(), &m_constraintValue ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Volumetric rate (if useSurfaceConditions: [surface m^3/s]; else [reservoir m^3/s])" );

}


VolumeConstraint::~VolumeConstraint()
{}

void VolumeConstraint::postInputInitialization()
{
  WellConstraintBase::postInputInitialization();

}


VolumeProductionConstraint::VolumeProductionConstraint( string const & name, Group * const parent )
  : VolumeConstraint( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
  m_rateSign=-1.0;

}


VolumeProductionConstraint::~VolumeProductionConstraint()
{}

void VolumeProductionConstraint::postInputInitialization()
{
  VolumeConstraint::postInputInitialization();

}

bool VolumeProductionConstraint::checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime )const
{
  return currentConstraint.totalVolumeRate() <  getConstraintValue( currentTime );
}


VolumeInjectionConstraint::VolumeInjectionConstraint( string const & name, Group * const parent )
  : VolumeConstraint( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );


  // Field registration
  registerInjectionStream( m_injectionStream, m_injectionTemperature, *this );
}


VolumeInjectionConstraint::~VolumeInjectionConstraint()
{}

void VolumeInjectionConstraint::postInputInitialization()
{

  VolumeConstraint::postInputInitialization();
// Validate the injection stream and temperature
  validateInjectionStream( m_injectionStream, m_injectionTemperature, getConstraintKey(), *this );

}

bool VolumeInjectionConstraint::checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime )const
{
  return currentConstraint.totalVolumeRate() >  getConstraintValue( currentTime );
}

} //namespace geos
