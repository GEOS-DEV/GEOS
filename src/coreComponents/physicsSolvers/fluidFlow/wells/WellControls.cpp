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

/*
 * @file WellControls.cpp
 */

#include "WellControls.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionManager.hpp"

namespace geosx
{

using namespace dataRepository;

WellControls::WellControls( string const & name, Group * const parent )
  : Group( name, parent ),
  m_type( Type::PRODUCER ),
  m_refElevation( 0.0 ),
  m_refGravCoef( 0.0 ),
  m_currentControl( Control::BHP ),
  m_targetBHP( 0.0 ),
  m_targetTotalRate( 0.0 ),
  m_targetPhaseRate( 0.0 ),
  m_useSurfaceConditions( 0 ),
  m_surfacePres( 0.0 ),
  m_surfaceTemp( 0.0 ),
  m_disableInjectorCrossflow( 0 ),
  m_targetTotalRateTable( nullptr ),
  m_targetPhaseRateTable( nullptr ),
  m_targetBHPTable( nullptr )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::typeString(), &m_type ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Well type. Valid options:\n* " + EnumStrings< Type >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::controlString(), &m_currentControl ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Well control. Valid options:\n* " + EnumStrings< Control >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::targetBHPString(), &m_targetBHP ).
    setDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Target bottom-hole pressure" );

  registerWrapper( viewKeyStruct::targetTotalRateString(), &m_targetTotalRate ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Target total volumetric rate" );

  registerWrapper( viewKeyStruct::targetPhaseRateString(), &m_targetPhaseRate ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Target phase volumetric rate" );

  registerWrapper( viewKeyStruct::targetPhaseNameString(), &m_targetPhaseName ).
    setDefaultValue( "" ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the target phase" );

  registerWrapper( viewKeyStruct::refElevString(), &m_refElevation ).
    setDefaultValue( -1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference elevation where BHP control is enforced" );

  registerWrapper( viewKeyStruct::injectionStreamString(), &m_injectionStream ).
    setDefaultValue( -1 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Global component densities of the injection stream" );

  registerWrapper( viewKeyStruct::injectionTemperatureString(), &m_injectionTemperature ).
    setDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Temperature of the injection stream" );

  registerWrapper( viewKeyStruct::useSurfaceConditionsString(), &m_useSurfaceConditions ).
    setDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to specify whether rates are checked at surface or reservoir conditions.\n"
                    "Equal to 1 for surface conditions, and to 0 for reservoir conditions" );

  registerWrapper( viewKeyStruct::surfacePressureString(), &m_surfacePres ).
    setDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Surface pressure used to compute volumetric rates when surface conditions are used" );

  registerWrapper( viewKeyStruct::surfaceTemperatureString(), &m_surfaceTemp ).
    setDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Surface temperature used to compute volumetric rates when surface conditions are used" );

  registerWrapper( viewKeyStruct::disableInjectorCrossflowString(), &m_disableInjectorCrossflow ).
    setDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to disable crossflow in injecting wells: \n"
                    " - If the flag is set to 0, both reservoir-to-well flow and well-to-reservoir flow are allowed at the perforations. \n"
                    " - If the flag is set to 1, we only allow well-to-reservoir flow at the perforations. \n"
                    "This flag is only available for injectors." );

  registerWrapper( viewKeyStruct::targetBHPTableNameString(), &m_targetBHPTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the BHP table when the rate is a time dependent function" );

  registerWrapper( viewKeyStruct::targetTotalRateTableNameString(), &m_targetTotalRateTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the total rate table when the rate is a time dependent function" );

  registerWrapper( viewKeyStruct::targetPhaseRateTableNameString(), &m_targetPhaseRateTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the phase rate table when the rate is a time dependent function" );
}


WellControls::~WellControls()
{}

void WellControls::switchToBHPControl( real64 const & val )
{
  m_currentControl = Control::BHP;
  m_targetBHP = val;
}

void WellControls::switchToTotalRateControl( real64 const & val )
{
  m_currentControl = Control::TOTALVOLRATE;
  m_targetTotalRate = val;
}

void WellControls::switchToPhaseRateControl( real64 const & val )
{
  m_currentControl = Control::PHASEVOLRATE;
  m_targetPhaseRate = val;
}

namespace
{

/// Utility function to create a one-value table internally when not provided by the user
TableFunction * createWellTable( string const & tableName,
                                 real64 const & constantValue )
{
  array1d< array1d< real64 > > timeCoord;
  timeCoord.resize( 1 );
  timeCoord[0].emplace_back( 0 );
  array1d< real64 > constantValueArray;
  constantValueArray.emplace_back( constantValue );

  FunctionManager & functionManager = FunctionManager::getInstance();
  TableFunction * table = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", tableName ));
  table->setTableCoordinates( timeCoord );
  table->setTableValues( constantValueArray );
  table->setInterpolationMethod( TableFunction::InterpolationType::Lower );
  return table;
}

}

void WellControls::postProcessInput()
{
  // 3.a) check target BHP
  GEOSX_THROW_IF( m_targetBHP < 0,
                  "WellControls '" << getName() << "': Target bottom-hole pressure is negative",
                  InputError );

  // 3.b) check target rates
  GEOSX_THROW_IF( m_targetTotalRate < 0,
                  "WellControls '" << getName() << "': Target rate is negative",
                  InputError );

  GEOSX_THROW_IF( m_targetPhaseRate < 0,
                  "WellControls '" << getName() << "': Target oil rate is negative",
                  InputError );

  GEOSX_THROW_IF( (m_injectionStream.empty()  && m_injectionTemperature >= 0) ||
                  (!m_injectionStream.empty() && m_injectionTemperature < 0),
                  "WellControls '" << getName() << "': Both "
                                   << viewKeyStruct::injectionStreamString() << " and " << viewKeyStruct::injectionTemperatureString()
                                   << " must be specified for multiphase simulations",
                  InputError );

  // 4) check injection stream
  if( !m_injectionStream.empty())
  {
    real64 sum = 0.0;
    for( localIndex ic = 0; ic < m_injectionStream.size(); ++ic )
    {
      GEOSX_ERROR_IF( m_injectionStream[ic] < 0.0 || m_injectionStream[ic] > 1.0,
                      "WellControls '" << getName() << "': Invalid injection stream" );
      sum += m_injectionStream[ic];
    }
    GEOSX_THROW_IF( LvArray::math::abs( 1.0 - sum ) > std::numeric_limits< real64 >::epsilon(),
                    "WellControls '" << getName() << "': Invalid injection stream",
                    InputError );
  }

  // 5) check the flag for surface / reservoir conditions
  GEOSX_THROW_IF( m_useSurfaceConditions != 0 && m_useSurfaceConditions != 1,
                  "WellControls '" << getName() << "': The flag to select surface/reservoir conditions must be equal to 0 or 1",
                  InputError );

  // 6) check the flag for surface / reservoir conditions
  GEOSX_THROW_IF( m_useSurfaceConditions == 1 && m_surfacePres <= 0,
                  "WellControls '" << getName() << "': When useSurfaceConditions == 1, the surface pressure must be defined",
                  InputError );

  // 7) check that at least one rate constraint has been defined
  GEOSX_THROW_IF( ((m_targetPhaseRate <= 0.0 && m_targetPhaseRateTableName.empty()) &&
                   (m_targetTotalRate <= 0.0 && m_targetTotalRateTableName.empty())),
                  "WellControls '" << getName() << "': You need to specify a phase rate constraint or a total rate constraint for injectors",
                  InputError );

  // 8) check whether redundant information has been provided
  GEOSX_THROW_IF( ((m_targetPhaseRate > 0.0 && !m_targetPhaseRateTableName.empty())),
                  "WellControls '" << getName() << "': You have provided redundant information for well phase rate." <<
                  " The keywords " << viewKeyStruct::targetPhaseRateString() << " and " << viewKeyStruct::targetPhaseRateTableNameString() << " cannot be specified together",
                  InputError );

  GEOSX_THROW_IF( ((m_targetTotalRate > 0.0 && !m_targetTotalRateTableName.empty())),
                  "WellControls '" << getName() << "': You have provided redundant information for well total rate." <<
                  " The keywords " << viewKeyStruct::targetTotalRateString() << " and " << viewKeyStruct::targetTotalRateTableNameString() << " cannot be specified together",
                  InputError );

  GEOSX_THROW_IF( ((m_targetBHP > 0.0 && !m_targetBHPTableName.empty())),
                  "WellControls '" << getName() << "': You have provided redundant information for well BHP." <<
                  " The keywords " << viewKeyStruct::targetBHPString() << " and " << viewKeyStruct::targetBHPTableNameString() << " cannot be specified together",
                  InputError );

  GEOSX_THROW_IF( ((m_targetBHP <= 0.0 && !m_targetBHPTableName.empty())),
                  "WellControls '" << getName() << "': You have to provide well BHP by specifying either " << viewKeyStruct::targetBHPString() << " or " << viewKeyStruct::targetBHPTableNameString(),
                  InputError );

  //  9) Make sure that the flag disabling crossflow is not used for producers
  GEOSX_THROW_IF( getType() == Type::PRODUCER && m_disableInjectorCrossflow != 0,
                  "WellControls '" << getName() << "': The option '" << viewKeyStruct::disableInjectorCrossflowString() << "' is not available for producers",
                  InputError );

  //  10) Create time-dependent BHP table
  if( m_targetBHPTableName.empty() )
  {
    m_targetBHPTable = createWellTable( getName()+"_ConstantBHP_table", m_targetBHP );
  }
  else
  {
    FunctionManager & functionManager = FunctionManager::getInstance();
    m_targetBHPTable = &(functionManager.getGroup< TableFunction >( m_targetBHPTableName ));

    GEOSX_THROW_IF( m_targetBHPTable->getInterpolationMethod() != TableFunction::InterpolationType::Lower,
                    "WellControls '" << getName() << "': The interpolation method for the time-dependent rate table "
                                     << m_targetBHPTable->getName() << " should be TableFunction::InterpolationType::Lower",
                    InputError );
  }

  //  11) Create time-dependent total rate table
  if( m_targetTotalRateTableName.empty() )
  {
    m_targetTotalRateTable = createWellTable( getName()+"_ConstantTotalRate_table", m_targetTotalRate );
  }
  else
  {
    FunctionManager & functionManager = FunctionManager::getInstance();
    m_targetTotalRateTable = &(functionManager.getGroup< TableFunction >( m_targetTotalRateTableName ));

    GEOSX_THROW_IF( m_targetTotalRateTable->getInterpolationMethod() != TableFunction::InterpolationType::Lower,
                    "WellControls '" << getName() << "': The interpolation method for the time-dependent rate table "
                                     << m_targetTotalRateTable->getName() << " should be TableFunction::InterpolationType::Lower",
                    InputError );
  }

  //  12) Create time-dependent phase rate table
  if( m_targetPhaseRateTableName.empty() )
  {
    m_targetPhaseRateTable = createWellTable( getName()+"_ConstantPhaseRate_table", m_targetPhaseRate );
  }
  else
  {
    FunctionManager & functionManager = FunctionManager::getInstance();
    m_targetPhaseRateTable = &(functionManager.getGroup< TableFunction >( m_targetPhaseRateTableName ));

    GEOSX_THROW_IF( m_targetPhaseRateTable->getInterpolationMethod() != TableFunction::InterpolationType::Lower,
                    "WellControls '" << getName() << "': The interpolation method for the time-dependent rate table "
                                     << m_targetPhaseRateTable->getName() << " should be TableFunction::InterpolationType::Lower",
                    InputError );
  }

  // For producer, the user has specified a positive rate
  // Therefore, we multiply it by -1 before the simulation starts to have the correct sign in the equations
  if( getType() == Type::PRODUCER )
  {
    array1d< real64 > & phaseRatetableValues = m_targetPhaseRateTable->getValues();
    for( localIndex i = 0; i < phaseRatetableValues.size(); ++i )
    {
      phaseRatetableValues( i ) *= -1;
    }

    array1d< real64 > & totalRatetableValues = m_targetTotalRateTable->getValues();
    for( localIndex i = 0; i < totalRatetableValues.size(); ++i )
    {
      totalRatetableValues( i ) *= -1;
    }
  }
}

bool WellControls::wellIsOpen( real64 const & currentTime ) const
{
  bool isOpen = true;
  if( ( m_currentControl == Control::TOTALVOLRATE && isZero( getTargetTotalRate( currentTime ) ) ) ||
      ( m_currentControl == Control::PHASEVOLRATE && isZero( getTargetPhaseRate( currentTime ) ) ) )
  {
    isOpen = false;
  }
  return isOpen;
}


} //namespace geosx
