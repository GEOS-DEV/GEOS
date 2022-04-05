/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
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

namespace geosx
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
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the phase initially saturating the reservoir" );

  registerWrapper( viewKeyStruct::componentNamesString(), &m_componentNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Names of the fluid components" );

  registerWrapper( viewKeyStruct::componentFractionVsElevationTableNamesString(), &m_componentFractionVsElevationTableNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Names of the tables specifying the (component fraction vs elevation) relationship for each component" );

  registerWrapper( viewKeyStruct::temperatureVsElevationTableNameString(), &m_temperatureVsElevationTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the table specifying the (temperature [K] vs elevation) relationship" );

  getWrapper< string >( FieldSpecificationBase::viewKeyStruct::fieldNameString() ).
    setInputFlag( InputFlags::FALSE );
  setFieldName( catalogName() );

  getWrapper< int >( FieldSpecificationBase::viewKeyStruct::componentString() ).
    setInputFlag( InputFlags::FALSE );

  getWrapper< int >( FieldSpecificationBase::viewKeyStruct::initialConditionString() ).
    setInputFlag( InputFlags::FALSE );
  initialCondition( false ); // to make sure this is not called by applyInitialConditions

  getWrapper< string_array >( FieldSpecificationBase::viewKeyStruct::setNamesString() ).
    setInputFlag( InputFlags::FALSE );
  addSetName( "all" );
}

void EquilibriumInitialCondition::postProcessInput()
{

  GEOSX_THROW_IF( ( m_temperatureVsElevationTableName.empty() != m_componentFractionVsElevationTableNames.empty() ),
                  getCatalogName() << " " << getName()
                                   << ": both " <<
                  viewKeyStruct::componentFractionVsElevationTableNamesString() << " and " <<
                  viewKeyStruct::temperatureVsElevationTableNameString() << " must be provided for a multiphase simulation",
                  InputError );

  FunctionManager const & functionManager = FunctionManager::getInstance();

  if( !m_componentFractionVsElevationTableNames.empty() )
  {
    GEOSX_THROW_IF( m_componentFractionVsElevationTableNames.size() <= 1,
                    getCatalogName() << " " << getName() << ": at least two component names must be specified in "
                                     << viewKeyStruct::componentNamesString(),
                    InputError );
    GEOSX_THROW_IF( m_componentFractionVsElevationTableNames.size() != m_componentNames.size(),
                    getCatalogName() << " " << getName() << ": mismatch between the size of "
                                     << viewKeyStruct::componentNamesString() << " and " << viewKeyStruct::componentFractionVsElevationTableNamesString(),
                    InputError );
    GEOSX_THROW_IF( m_componentNames.size() >= 2 && m_initPhaseName.empty(),
                    getCatalogName() << " " << getName() << ": for now, the keyword: " << viewKeyStruct::initPhaseNameString()
                                     << " must be filled for a multiphase simulation",
                    InputError );

    array1d< localIndex > tableSizes( m_componentNames.size() );
    for( localIndex ic = 0; ic < m_componentNames.size(); ++ic )
    {
      GEOSX_THROW_IF( m_componentFractionVsElevationTableNames[ic].empty(),
                      getCatalogName() << " " << getName()
                                       << ": the component fraction vs elevation table name is missing for component " << ic,
                      InputError );

      GEOSX_THROW_IF( !m_componentFractionVsElevationTableNames[ic].empty() &&
                      !functionManager.hasGroup( m_componentFractionVsElevationTableNames[ic] ),
                      getCatalogName() << " " << getName()
                                       << ": the component fraction vs elevation table " << m_componentFractionVsElevationTableNames[ic]
                                       << " could not be found"  << " for component " << ic,
                      InputError );

      TableFunction const & compFracTable = functionManager.getGroup< TableFunction >( m_componentFractionVsElevationTableNames[ic] );
      GEOSX_THROW_IF( compFracTable.getInterpolationMethod() != TableFunction::InterpolationType::Linear,
                      getCatalogName() << " " << getName() << ": the interpolation method for the component fraction vs elevation table "
                                       << compFracTable.getName() << " should be TableFunction::InterpolationType::Linear",
                      InputError );

    }
  }

  if( !m_temperatureVsElevationTableName.empty() )
  {

    GEOSX_THROW_IF( !functionManager.hasGroup( m_temperatureVsElevationTableName ),
                    getCatalogName() << " " << getName() << ": the temperature vs elevation table " << m_temperatureVsElevationTableName << " could not be found",
                    InputError );

    TableFunction const & tempTable = functionManager.getGroup< TableFunction >( m_temperatureVsElevationTableName );
    GEOSX_THROW_IF( tempTable.getInterpolationMethod() != TableFunction::InterpolationType::Linear,
                    getCatalogName() << " " << getName() << ": The interpolation method for the temperature vs elevation table "
                                     << tempTable.getName() << " should be TableFunction::InterpolationType::Linear",
                    InputError );
  }
}

void EquilibriumInitialCondition::initializePreSubGroups()
{

  if( !m_componentFractionVsElevationTableNames.empty() )
  {

    FunctionManager const & functionManager = FunctionManager::getInstance();

    array1d< localIndex > tableSizes( m_componentNames.size() );
    for( localIndex ic = 0; ic < m_componentNames.size(); ++ic )
    {
      TableFunction const & compFracTable = functionManager.getGroup< TableFunction >( m_componentFractionVsElevationTableNames[ic] );
      arrayView1d< real64 const > compFracValues = compFracTable.getValues();
      GEOSX_THROW_IF( compFracValues.size() <= 1,
                      getCatalogName() << " " << getName() << ": the component fraction vs elevation table "
                                       << compFracTable.getName() << " must contain at least two values",
                      InputError );

      tableSizes[ic] = compFracValues.size();
      if( ic >= 1 )
      {
        GEOSX_THROW_IF( tableSizes[ic] != tableSizes[ic-1],
                        getCatalogName() << " " << getName()
                                         << ": all the component fraction vs elevation tables must contain the same number of values",
                        InputError );
      }
    }

    array2d< real64 > elevation( m_componentNames.size(), tableSizes[0] );
    array1d< real64 > sumCompFrac( tableSizes[0] );
    for( localIndex ic = 0; ic < m_componentNames.size(); ++ic )
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
          GEOSX_THROW_IF( !isZero( elevation[ic][i] - elevation[ic-1][i] ),
                          getCatalogName() << " " << getName()
                                           << ": the elevation values must be the same in all the component vs elevation tables",
                          InputError );
        }

        if( ic == m_componentNames.size() - 1 )
        {
          GEOSX_THROW_IF( !isZero( sumCompFrac[i] - 1 ),
                          getCatalogName() << " " << getName()
                                           << ": at a given elevation, the component fraction sum must be equal to one",
                          InputError );
        }
      }
    }
  }
}

REGISTER_CATALOG_ENTRY( FieldSpecificationBase, EquilibriumInitialCondition, string const &, Group * const )


} /* namespace geosx */
