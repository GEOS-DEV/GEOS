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

/**
 * @file SolidThermalExpansion.cpp
 */

#include "SolidThermalExpansion.hpp"

#include "codingUtilities/Utilities.hpp"
#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

SolidThermalExpansion::SolidThermalExpansion( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_referenceDrainedLinearTEC( 0.0 ),
  m_dDrainedLinearTEC_dTemperature( 0.0 ),
  m_referenceTemperature( 0.0 )
{
  registerWrapper( viewKeyStruct::referenceDrainedLinearTECString(), &m_referenceDrainedLinearTEC ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Drained linear thermal expansion coefficient at the reference temperature [1/temperature]. "
                    "Required unless drainedLinearTECTableName is given." );

  registerWrapper( viewKeyStruct::dDrainedLinearTEC_dTemperatureString(), &m_dDrainedLinearTEC_dTemperature ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Derivative of the thermal expansion coefficient with respect to temperature [1/temperature^2]. "
                    "The coefficient is referenceDrainedLinearTEC + dDrainedLinearTEC_dTemperature * ( T - referenceTemperature )." );

  registerWrapper( viewKeyStruct::referenceTemperatureString(), &m_referenceTemperature ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Temperature at which the thermal expansion coefficient equals referenceDrainedLinearTEC [temperature]. "
                    "Required when dDrainedLinearTEC_dTemperature is non-zero. "
                    "This is a property of the material, not the initial temperature of the simulation." );

  registerWrapper( viewKeyStruct::drainedLinearTECTableNameString(), &m_drainedLinearTECTableName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "" ).
    setDescription( "Name of a one-dimensional TableFunction of temperature, with linear interpolation, giving the thermal "
                    "expansion coefficient. Cannot be combined with referenceDrainedLinearTEC, dDrainedLinearTEC_dTemperature or referenceTemperature." );

  // Copies of the table data, filled during initialization. They are registered so that the
  // per-subregion copies of this model carry them as well; they are never written to restart files.
  registerWrapper( viewKeyStruct::tableKnotsString(), &m_tableKnots ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( 0 ).
    setDescription( "Table temperatures, copied from the TableFunction [temperature]" );

  registerWrapper( viewKeyStruct::tableValuesString(), &m_tableValues ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( 0 ).
    setDescription( "Thermal expansion coefficient at each table temperature, copied from the TableFunction [1/temperature]" );

  registerWrapper( viewKeyStruct::tableCumulativeIntegralString(), &m_tableCumulativeIntegral ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( 0 ).
    setDescription( "Integral of the thermal expansion coefficient from the first table temperature to each table temperature [dimensionless]" );
}

void SolidThermalExpansion::postInputInitialization()
{
  ConstitutiveBase::postInputInitialization();

  bool const hasReferenceDrainedLinearTEC = getWrapperBase( viewKeyStruct::referenceDrainedLinearTECString() ).getSuccessfulReadFromInput();
  bool const hasDerivative = getWrapperBase( viewKeyStruct::dDrainedLinearTEC_dTemperatureString() ).getSuccessfulReadFromInput();
  bool const hasReferenceTemperature = getWrapperBase( viewKeyStruct::referenceTemperatureString() ).getSuccessfulReadFromInput();

  if( !m_drainedLinearTECTableName.empty() )
  {
    GEOS_THROW_IF( hasReferenceDrainedLinearTEC || hasDerivative || hasReferenceTemperature,
                   GEOS_FMT( "{} cannot be combined with {}, {} or {}",
                             viewKeyStruct::drainedLinearTECTableNameString(),
                             viewKeyStruct::referenceDrainedLinearTECString(),
                             viewKeyStruct::dDrainedLinearTEC_dTemperatureString(),
                             viewKeyStruct::referenceTemperatureString() ),
                   InputError, getWrapperDataContext( viewKeyStruct::drainedLinearTECTableNameString() ) );
    return;
  }

  GEOS_THROW_IF( !hasReferenceDrainedLinearTEC,
                 GEOS_FMT( "Either {} or {} must be given",
                           viewKeyStruct::referenceDrainedLinearTECString(),
                           viewKeyStruct::drainedLinearTECTableNameString() ),
                 InputError, getDataContext() );

  GEOS_THROW_IF( !isZero( m_dDrainedLinearTEC_dTemperature, 0.0 ) && !hasReferenceTemperature,
                 GEOS_FMT( "{} must be given when {} is non-zero",
                           viewKeyStruct::referenceTemperatureString(),
                           viewKeyStruct::dDrainedLinearTEC_dTemperatureString() ),
                 InputError, getWrapperDataContext( viewKeyStruct::dDrainedLinearTEC_dTemperatureString() ) );
}

void SolidThermalExpansion::initializePreSubGroups()
{
  ConstitutiveBase::initializePreSubGroups();

  if( m_drainedLinearTECTableName.empty() )
  {
    return;
  }

  // TableFunctions are initialized when the input is read, before this step.
  FunctionManager const & functionManager = FunctionManager::getInstance();
  GEOS_THROW_IF( !functionManager.hasGroup< TableFunction >( m_drainedLinearTECTableName ),
                 GEOS_FMT( "The TableFunction '{}' could not be found", m_drainedLinearTECTableName ),
                 InputError, getWrapperDataContext( viewKeyStruct::drainedLinearTECTableNameString() ) );

  TableFunction const & table = functionManager.getGroup< TableFunction >( m_drainedLinearTECTableName );

  GEOS_THROW_IF( table.numDimensions() != 1,
                 GEOS_FMT( "The TableFunction '{}' must have one input (temperature), it has {}",
                           m_drainedLinearTECTableName, table.numDimensions() ),
                 InputError, getWrapperDataContext( viewKeyStruct::drainedLinearTECTableNameString() ) );

  GEOS_THROW_IF( table.getInterpolationMethod() != TableFunction::InterpolationType::Linear,
                 GEOS_FMT( "The TableFunction '{}' must use linear interpolation", m_drainedLinearTECTableName ),
                 InputError, getWrapperDataContext( viewKeyStruct::drainedLinearTECTableNameString() ) );

  ArrayOfArraysView< real64 const > const coordinates = table.getCoordinates();
  arrayView1d< real64 const > const values = table.getValues();
  localIndex const numKnots = coordinates.sizeOfArray( 0 );

  GEOS_THROW_IF( numKnots < 1 || values.size() != numKnots,
                 GEOS_FMT( "The TableFunction '{}' must have as many values as temperatures, and at least one",
                           m_drainedLinearTECTableName ),
                 InputError, getWrapperDataContext( viewKeyStruct::drainedLinearTECTableNameString() ) );

  m_tableKnots.resize( numKnots );
  m_tableValues.resize( numKnots );
  m_tableCumulativeIntegral.resize( numKnots );

  for( localIndex i = 0; i < numKnots; ++i )
  {
    m_tableKnots[i] = coordinates( 0, i );
    m_tableValues[i] = values[i];
  }

  // Exact integral of the piecewise linear coefficient from the first knot to each knot.
  m_tableCumulativeIntegral[0] = 0.0;
  for( localIndex i = 1; i < numKnots; ++i )
  {
    m_tableCumulativeIntegral[i] = m_tableCumulativeIntegral[i - 1]
                                   + 0.5 * ( m_tableValues[i - 1] + m_tableValues[i] ) * ( m_tableKnots[i] - m_tableKnots[i - 1] );
  }
}

bool SolidThermalExpansion::isTemperatureDependent() const
{
  return !m_drainedLinearTECTableName.empty() || !isZero( m_dDrainedLinearTEC_dTemperature, 0.0 );
}

SolidThermalExpansion::KernelWrapper
SolidThermalExpansion::createKernelUpdates() const
{
  if( m_drainedLinearTECTableName.empty() )
  {
    return KernelWrapper( m_referenceDrainedLinearTEC, m_dDrainedLinearTEC_dTemperature, m_referenceTemperature );
  }

  GEOS_ERROR_IF( m_tableKnots.empty(),
                 GEOS_FMT( "The data of the TableFunction '{}' has not been copied yet", m_drainedLinearTECTableName ),
                 getDataContext() );

  return KernelWrapper( m_tableKnots.toViewConst(),
                        m_tableValues.toViewConst(),
                        m_tableCumulativeIntegral.toViewConst() );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, SolidThermalExpansion, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
