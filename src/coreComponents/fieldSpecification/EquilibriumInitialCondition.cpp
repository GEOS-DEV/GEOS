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
 * @file EquilibriumInitialCondition.cpp
 */

#include "EquilibriumInitialCondition.hpp"

#include "functions/TableFunction.hpp"

namespace geos
{

using namespace dataRepository;

EquilibriumInitialCondition::EquilibriumInitialCondition( string const & name, Group * parent ):
  FieldSpecificationBase( name, parent )
{
  registerWrapper( viewKeyStruct::datumElevationString(), &m_datumElevation ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Datum elevation [m]" );

  registerWrapper( viewKeyStruct::datumPressureString(), &m_datumPressure ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Datum pressure [Pa]" );

  registerWrapper( viewKeyStruct::maxNumEquilibrationIterationsString(), &m_maxNumEquilibrationIterations ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 5 ).
    setDescription( "Maximum number of equilibration iterations" );

  registerWrapper( viewKeyStruct::equilibrationToleranceString(), &m_equilibrationTolerance ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1e-3 ).
    setDescription( "Tolerance in the fixed-point iteration scheme used for hydrostatic initialization" );

  registerWrapper( viewKeyStruct::elevationIncrementString(), &m_elevationIncrement ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.6096 ). // 2 feet
    setDescription( "Elevation increment [m] in the hydrostatic pressure table constructed internally" );

  registerWrapper( viewKeyStruct::initPhaseNameString(), &m_initPhaseName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the phase initially saturating the reservoir" );

  registerWrapper( viewKeyStruct::componentNamesString(), &m_componentNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Names of the fluid components" );

  registerWrapper( viewKeyStruct::componentFractionVsElevationTableNamesString(), &m_componentFractionVsElevationTableNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Names of the tables specifying the (component fraction vs elevation) relationship for each component" );

  registerWrapper( viewKeyStruct::temperatureVsElevationTableNameString(), &m_temperatureVsElevationTableName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the table specifying the (temperature [K] vs elevation) relationship" );

  registerWrapper( viewKeyStruct::phaseContactsString(), &m_phaseContacts ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Phase contacts' elevations [m]" );

  getWrapper< string >( FieldSpecificationBase::viewKeyStruct::fieldNameString() ).
    setInputFlag( InputFlags::FALSE );
  setFieldName( catalogName() );

  getWrapper< int >( FieldSpecificationBase::viewKeyStruct::componentString() ).
    setInputFlag( InputFlags::FALSE );

  getWrapper< int >( FieldSpecificationBase::viewKeyStruct::initialConditionString() ).
    setInputFlag( InputFlags::FALSE );
  initialCondition( false ); // to make sure this is not called by applyInitialConditions

  getWrapper< string_array >( FieldSpecificationBase::viewKeyStruct::setNamesString() ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::FALSE );
  addSetName( "all" );
}

void EquilibriumInitialCondition::postInputInitialization()
{
  FieldSpecificationBase::postInputInitialization();

  FunctionManager const & functionManager = FunctionManager::getInstance();

  if( !m_componentFractionVsElevationTableNames.empty() )
  {
    GEOS_THROW_IF( m_componentFractionVsElevationTableNames.size() <= 1,
                   GEOS_FMT( "At least two component names must be specified in {}",
                             viewKeyStruct::componentNamesString() ),
                   InputError, getDataContext() );
    GEOS_THROW_IF( m_componentFractionVsElevationTableNames.size() != m_componentNames.size(),
                   GEOS_FMT( "Mismatch between the size of {} and {}",
                             viewKeyStruct::componentNamesString(),
                             viewKeyStruct::componentFractionVsElevationTableNamesString() ),
                   InputError, getDataContext() );

    integer const numberOfComponents = static_cast< integer >(m_componentNames.size());

    if( 1 < numberOfComponents )
    {
      integer const numberOfContacts = static_cast< integer >(m_phaseContacts.size());
      GEOS_THROW_IF( m_initPhaseName.empty() && numberOfContacts == 0,
                     GEOS_FMT( ": for a multiphase simulation either the initial phase name must be provided using {} "
                               "or the phase contact elevations number be provided using {}",
                               viewKeyStruct::initPhaseNameString(),
                               viewKeyStruct::phaseContactsString() ),
                     InputError, getDataContext() );

      if( !m_initPhaseName.empty() && 0 < numberOfContacts )
      {
        GEOS_WARNING( GEOS_FMT( "both {} and {} have been specified. The phase contacts will be ignored and "
                                "single phase initialisation performed",
                                viewKeyStruct::initPhaseNameString(),
                                viewKeyStruct::phaseContactsString() ),
                      getDataContext() );
      }

      // Contacts if provided must be non-decreasing
      if( 1 < numberOfContacts )
      {
        for( integer i = 1; i < numberOfContacts; i++ )
        {
          GEOS_THROW_IF( m_phaseContacts[i] - m_phaseContacts[i-1] < -LvArray::NumericLimits< real64 >::epsilon,
                         "The phase contacts must be increasing",
                         InputError, getDataContext() );
        }
      }
    }

    array1d< localIndex > tableSizes( numberOfComponents );
    for( integer ic = 0; ic < numberOfComponents; ++ic )
    {
      GEOS_THROW_IF( m_componentFractionVsElevationTableNames[ic].empty(),
                     GEOS_FMT( "The component fraction vs elevation table name is missing for component {}", ic ),
                     InputError, getDataContext() );

      GEOS_THROW_IF( !m_componentFractionVsElevationTableNames[ic].empty() &&
                     !functionManager.hasGroup( m_componentFractionVsElevationTableNames[ic] ),
                     GEOS_FMT( "The component fraction vs elevation table {} could not be found for component {}",
                               m_componentFractionVsElevationTableNames[ic],
                               ic ),
                     InputError, getDataContext() );

      TableFunction const & compFracTable = functionManager.getGroup< TableFunction >( m_componentFractionVsElevationTableNames[ic] );
      GEOS_THROW_IF( compFracTable.getInterpolationMethod() != TableFunction::InterpolationType::Linear,
                     GEOS_FMT( "The interpolation method for the component fraction vs elevation table {} "
                               "should be TableFunction::InterpolationType::Linear",
                               compFracTable.getName() ),
                     InputError, getDataContext() );

    }
  }

  if( !m_temperatureVsElevationTableName.empty() )
  {

    GEOS_THROW_IF( !functionManager.hasGroup( m_temperatureVsElevationTableName ),
                   GEOS_FMT( "The temperature vs elevation table {} could not be found",
                             m_temperatureVsElevationTableName ),
                   InputError, getDataContext() );

    TableFunction const & tempTable = functionManager.getGroup< TableFunction >( m_temperatureVsElevationTableName );
    GEOS_THROW_IF( tempTable.getInterpolationMethod() != TableFunction::InterpolationType::Linear,
                   GEOS_FMT( "The interpolation method for the temperature vs elevation table {} "
                             "should be TableFunction::InterpolationType::Linear",
                             tempTable.getName() ),
                   InputError, getDataContext() );
  }
}

void EquilibriumInitialCondition::initializePreSubGroups()
{

  if( !m_componentFractionVsElevationTableNames.empty() )
  {

    FunctionManager const & functionManager = FunctionManager::getInstance();

    array1d< localIndex > tableSizes( m_componentNames.size() );
    for( size_t ic = 0; ic < m_componentNames.size(); ++ic )
    {
      TableFunction const & compFracTable = functionManager.getGroup< TableFunction >( m_componentFractionVsElevationTableNames[ic] );
      arrayView1d< real64 const > compFracValues = compFracTable.getValues();
      GEOS_THROW_IF( compFracValues.size() <= 1,
                     GEOS_FMT( "The component fraction vs elevation table {} must contain at least two values",
                               compFracTable.getName() ),
                     InputError, getDataContext() );

      tableSizes[ic] = compFracValues.size();
      if( ic >= 1 )
      {
        GEOS_THROW_IF( tableSizes[ic] != tableSizes[ic-1],
                       "All the component fraction vs elevation tables must contain the same number of values",
                       InputError, getDataContext() );
      }
    }

    array2d< real64 > elevation( m_componentNames.size(), tableSizes[0] );
    array1d< real64 > sumCompFrac( tableSizes[0] );
    for( size_t ic = 0; ic < m_componentNames.size(); ++ic )
    {
      TableFunction const & compFracTable = functionManager.getGroup< TableFunction >( m_componentFractionVsElevationTableNames[ic] );

      ArrayOfArraysView< real64 const > elevationValues = compFracTable.getCoordinates();
      arrayView1d< real64 const > compFracValues = compFracTable.getValues();
      for( localIndex i = 0; i < compFracValues.size(); ++i )
      {
        elevation[ic][i] = elevationValues[0][i];
        sumCompFrac[i] += compFracValues[i];

        if( ic >= 1 )
        {
          GEOS_THROW_IF( !isZero( elevation[ic][i] - elevation[ic-1][i] ),
                         "The elevation values must be the same in all the component vs elevation tables",
                         InputError, getDataContext() );
        }

        if( ic == m_componentNames.size() - 1 )
        {
          GEOS_THROW_IF( !isZero( sumCompFrac[i] - 1 ),
                         "At a given elevation, the component fraction sum must be equal to one",
                         InputError, getDataContext() );
        }
      }
    }
  }
}

REGISTER_CATALOG_ENTRY( FieldSpecificationBase, EquilibriumInitialCondition, string const &, Group * const )


} /* namespace geos */
