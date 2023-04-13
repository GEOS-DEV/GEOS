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
#include "WellConstants.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionManager.hpp"

namespace geos
{

using namespace dataRepository;

WellControls::WellControls( string const & name, Group * const parent )
  : Group( name, parent ),
  m_type( Type::PRODUCER ),
  m_refElevation( 0.0 ),
  m_refGravCoef( 0.0 ),
  m_inputControl( Control::UNINITIALIZED ),
  m_currentControl( Control::UNINITIALIZED ),
  m_targetBHP( 0.0 ),
  m_targetTotalRate( 0.0 ),
  m_targetPhaseRate( 0.0 ),
  m_useSurfaceConditions( 0 ),
  m_surfacePres( 0.0 ),
  m_surfaceTemp( 0.0 ),
  m_isCrossflowEnabled( 1 ),
  m_initialPressureCoefficient( 0.1 ),
  m_rateSign( -1.0 ),
  m_targetTotalRateTable( nullptr ),
  m_targetPhaseRateTable( nullptr ),
  m_targetBHPTable( nullptr ),
  m_statusTable( nullptr )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  enableLogLevelInput();

  registerWrapper( viewKeyStruct::typeString(), &m_type ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Well type. Valid options:\n* " + EnumStrings< Type >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::inputControlString(), &m_inputControl ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Well control. Valid options:\n* " + EnumStrings< Control >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::currentControlString(), &m_currentControl ).
    setDefaultValue( Control::UNINITIALIZED ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Current well control" );

  registerWrapper( viewKeyStruct::targetBHPString(), &m_targetBHP ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Target bottom-hole pressure [Pa]" );

  registerWrapper( viewKeyStruct::targetTotalRateString(), &m_targetTotalRate ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Target total volumetric rate (if useSurfaceConditions: [surface m^3/s]; else [reservoir m^3/s])" );

  registerWrapper( viewKeyStruct::targetPhaseRateString(), &m_targetPhaseRate ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Target phase volumetric rate (if useSurfaceConditions: [surface m^3/s]; else [reservoir m^3/s])" );

  registerWrapper( viewKeyStruct::targetPhaseNameString(), &m_targetPhaseName ).
    setDefaultValue( "" ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the target phase" );

  registerWrapper( viewKeyStruct::refElevString(), &m_refElevation ).
    setDefaultValue( -1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference elevation where BHP control is enforced [m]" );

  registerWrapper( viewKeyStruct::injectionStreamString(), &m_injectionStream ).
    setDefaultValue( -1 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Global component densities of the injection stream [moles/m^3 or kg/m^3]" );

  registerWrapper( viewKeyStruct::injectionTemperatureString(), &m_injectionTemperature ).
    setDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Temperature of the injection stream [K]" );

  registerWrapper( viewKeyStruct::useSurfaceConditionsString(), &m_useSurfaceConditions ).
    setDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to specify whether rates are checked at surface or reservoir conditions.\n"
                    "Equal to 1 for surface conditions, and to 0 for reservoir conditions" );

  registerWrapper( viewKeyStruct::surfacePressureString(), &m_surfacePres ).
    setDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Surface pressure used to compute volumetric rates when surface conditions are used [Pa]" );

  registerWrapper( viewKeyStruct::surfaceTemperatureString(), &m_surfaceTemp ).
    setDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Surface temperature used to compute volumetric rates when surface conditions are used [K]" );

  registerWrapper( viewKeyStruct::enableCrossflowString(), &m_isCrossflowEnabled ).
    setDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to enable crossflow. Currently only supported for injectors: \n"
                    " - If the flag is set to 1, both reservoir-to-well flow and well-to-reservoir flow are allowed at the perforations. \n"
                    " - If the flag is set to 0, we only allow well-to-reservoir flow at the perforations." );

  registerWrapper( viewKeyStruct::initialPressureCoefficientString(), &m_initialPressureCoefficient ).
    setDefaultValue( 0.1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Tuning coefficient for the initial well pressure of rate-controlled wells: \n"
                    " - Injector pressure at reference depth initialized as: (1+initialPressureCoefficient)*reservoirPressureAtClosestPerforation + density*g*( zRef - zPerf ) \n"
                    " - Producer pressure at reference depth initialized as: (1-initialPressureCoefficient)*reservoirPressureAtClosestPerforation + density*g*( zRef - zPerf ) " );

  registerWrapper( viewKeyStruct::targetBHPTableNameString(), &m_targetBHPTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the BHP table when the rate is a time dependent function" );

  registerWrapper( viewKeyStruct::targetTotalRateTableNameString(), &m_targetTotalRateTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the total rate table when the rate is a time dependent function" );

  registerWrapper( viewKeyStruct::targetPhaseRateTableNameString(), &m_targetPhaseRateTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the phase rate table when the rate is a time dependent function" );

  registerWrapper( viewKeyStruct::statusTableNameString(), &m_statusTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the well status table when the status of the well is a time dependent function. \n"
                    "If the status function evaluates to a positive value at the current time, the well will be open otherwise the well will be shut." );
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
  TableFunction * table = dynamicCast< TableFunction * >( functionManager.createChild( TableFunction::catalogName(), tableName ));
  table->setTableCoordinates( timeCoord );
  table->setTableValues( constantValueArray );
  table->setInterpolationMethod( TableFunction::InterpolationType::Lower );
  return table;
}

}

void WellControls::postProcessInput()
{
  // 0) Assign the value of the current well control
  // When the simulation starts from a restart file, we don't want to use the inputControl,
  // because the control may have switched in the simulation that generated the restart
  GEOSX_THROW_IF_IF( m_inputControl == Control::UNINITIALIZED,
                  "WellControls '" << getName() << "': Input well control cannot be uninitialized",
                  InputError );

  if( m_currentControl == Control::UNINITIALIZED )
  {
    m_currentControl = m_inputControl;
  }

  // 1.a) check target BHP
  GEOSX_THROW_IF_IF( m_targetBHP < 0,
                  "WellControls '" << getName() << "': Target bottom-hole pressure is negative",
                  InputError );

  // 1.b) check target rates
  GEOSX_THROW_IF_IF( m_targetTotalRate < 0,
                  "WellControls '" << getName() << "': Target rate is negative",
                  InputError );

  GEOSX_THROW_IF_IF( m_targetPhaseRate < 0,
                  "WellControls '" << getName() << "': Target oil rate is negative",
                  InputError );

  GEOSX_THROW_IF_IF( (m_injectionStream.empty()  && m_injectionTemperature >= 0) ||
                  (!m_injectionStream.empty() && m_injectionTemperature < 0),
                  "WellControls '" << getName() << "': Both "
                                   << viewKeyStruct::injectionStreamString() << " and " << viewKeyStruct::injectionTemperatureString()
                                   << " must be specified for multiphase simulations",
                  InputError );

  // 1.c) Set the multiplier for the rates
  if( isProducer() )
  {
    m_rateSign = -1.0;
  }
  else
  {
    m_rateSign = 1.0;
  }

  // 2) check injection stream
  if( !m_injectionStream.empty())
  {
    real64 sum = 0.0;
    for( localIndex ic = 0; ic < m_injectionStream.size(); ++ic )
    {
      GEOS_ERROR_IF( m_injectionStream[ic] < 0.0 || m_injectionStream[ic] > 1.0,
                      "WellControls '" << getName() << "': Invalid injection stream" );
      sum += m_injectionStream[ic];
    }
    GEOSX_THROW_IF_IF( LvArray::math::abs( 1.0 - sum ) > std::numeric_limits< real64 >::epsilon(),
                    "WellControls '" << getName() << "': Invalid injection stream",
                    InputError );
  }

  // 3) check the flag for surface / reservoir conditions
  GEOSX_THROW_IF_IF( m_useSurfaceConditions != 0 && m_useSurfaceConditions != 1,
                  "WellControls '" << getName() << "': The flag to select surface/reservoir conditions must be equal to 0 or 1",
                  InputError );

  // 4) check the flag for surface / reservoir conditions
  GEOSX_THROW_IF_IF( m_useSurfaceConditions == 1 && m_surfacePres <= 0,
                  "WellControls '" << getName() << "': When " << viewKeyStruct::useSurfaceConditionsString() << " == 1, the surface pressure must be defined",
                  InputError );

  // 5) check that at least one rate constraint has been defined
  GEOSX_THROW_IF_IF( ((m_targetPhaseRate <= 0.0 && m_targetPhaseRateTableName.empty()) &&
                   (m_targetTotalRate <= 0.0 && m_targetTotalRateTableName.empty())),
                  "WellControls '" << getName() << "': You need to specify a phase rate constraint or a total rate constraint. \n" <<
                  "The phase rate constraint can be specified using " <<
                  "either " << viewKeyStruct::targetPhaseRateString() <<
                  " or " << viewKeyStruct::targetPhaseRateTableNameString() << ".\n" <<
                  "The total rate constraint can be specified using " <<
                  "either " << viewKeyStruct::targetTotalRateString() <<
                  " or " << viewKeyStruct::targetTotalRateTableNameString(),
                  InputError );

  // 6) check whether redundant information has been provided
  GEOSX_THROW_IF_IF( ((m_targetPhaseRate > 0.0 && !m_targetPhaseRateTableName.empty())),
                  "WellControls '" << getName() << "': You have provided redundant information for well phase rate." <<
                  " The keywords " << viewKeyStruct::targetPhaseRateString() << " and " << viewKeyStruct::targetPhaseRateTableNameString() << " cannot be specified together",
                  InputError );

  GEOSX_THROW_IF_IF( ((m_targetTotalRate > 0.0 && !m_targetTotalRateTableName.empty())),
                  "WellControls '" << getName() << "': You have provided redundant information for well total rate." <<
                  " The keywords " << viewKeyStruct::targetTotalRateString() << " and " << viewKeyStruct::targetTotalRateTableNameString() << " cannot be specified together",
                  InputError );

  GEOSX_THROW_IF_IF( ((m_targetBHP > 0.0 && !m_targetBHPTableName.empty())),
                  "WellControls '" << getName() << "': You have provided redundant information for well BHP." <<
                  " The keywords " << viewKeyStruct::targetBHPString() << " and " << viewKeyStruct::targetBHPTableNameString() << " cannot be specified together",
                  InputError );

  // 6.1) If the well is under BHP control then the BHP must be specified.
  //      Otherwise the BHP will be set to a default value.
  if( m_currentControl == Control::BHP )
  {
    GEOSX_THROW_IF_IF( ((m_targetBHP <= 0.0 && m_targetBHPTableName.empty())),
                    "WellControls '" << getName() << "': You have to provide well BHP by specifying either "
                                     << viewKeyStruct::targetBHPString() << " or " << viewKeyStruct::targetBHPTableNameString(),
                    InputError );
  }
  else if( m_targetBHP <= 0.0 && m_targetBHPTableName.empty() )
  {
    m_targetBHP = isProducer() ? WellConstants::defaultProducerBHP : WellConstants::defaultInjectorBHP;
    GEOS_LOG_LEVEL_RANK_0( 1, "WellControls '" << getName() << "': Setting " << viewKeyStruct::targetBHPString() << " to default value "
                                                << m_targetBHP << "." );
  }

  // 6.2) Check incoherent information
  // A producer must be controlled by PhaseVolRate
  GEOSX_THROW_IF_IF( (isProducer() && (m_inputControl == Control::TOTALVOLRATE)),
                  "WellControls '" << getName() << "': You have to control a producer with "
                                   << EnumStrings< Control >::toString( Control::PHASEVOLRATE ),
                  InputError );

  // An injector must be controlled by TotalVolRate
  GEOSX_THROW_IF_IF( (isInjector() && (m_inputControl == Control::PHASEVOLRATE)),
                  "WellControls '" << getName() << "': You have to control an injector with "
                                   << EnumStrings< Control >::toString( Control::TOTALVOLRATE ),
                  InputError );

  // 7) Make sure that the flag disabling crossflow is not used for producers
  GEOSX_THROW_IF_IF( isProducer() && m_isCrossflowEnabled == 0,
                  "WellControls '" << getName() << "': The option '" << viewKeyStruct::enableCrossflowString() << "' cannot be set to '0' for producers",
                  InputError );

  // 8) Make sure that the initial pressure coefficient is positive
  GEOSX_THROW_IF_IF( m_initialPressureCoefficient < 0,
                  "WellControls '" << getName() << "': The tuning coefficient " << viewKeyStruct::initialPressureCoefficientString() << " is negative",
                  InputError );


  // 9) Create time-dependent BHP table
  if( m_targetBHPTableName.empty() )
  {
    m_targetBHPTableName = getName()+"_ConstantBHP_table";
    m_targetBHPTable = createWellTable( m_targetBHPTableName, m_targetBHP );
  }
  else
  {
    FunctionManager & functionManager = FunctionManager::getInstance();
    m_targetBHPTable = &(functionManager.getGroup< TableFunction const >( m_targetBHPTableName ));

    GEOSX_THROW_IF_IF( m_targetBHPTable->getInterpolationMethod() != TableFunction::InterpolationType::Lower,
                    "WellControls '" << getName() << "': The interpolation method for the time-dependent rate table "
                                     << m_targetBHPTable->getName() << " should be TableFunction::InterpolationType::Lower",
                    InputError );
  }

  // 10) Create time-dependent total rate table
  if( m_targetTotalRateTableName.empty() )
  {
    m_targetTotalRateTableName = getName()+"_ConstantTotalRate_table";
    m_targetTotalRateTable = createWellTable( m_targetTotalRateTableName, m_targetTotalRate );
  }
  else
  {
    FunctionManager & functionManager = FunctionManager::getInstance();
    m_targetTotalRateTable = &(functionManager.getGroup< TableFunction const >( m_targetTotalRateTableName ));

    GEOSX_THROW_IF_IF( m_targetTotalRateTable->getInterpolationMethod() != TableFunction::InterpolationType::Lower,
                    "WellControls '" << getName() << "': The interpolation method for the time-dependent rate table "
                                     << m_targetTotalRateTable->getName() << " should be TableFunction::InterpolationType::Lower",
                    InputError );
  }

  // 11) Create time-dependent phase rate table
  if( m_targetPhaseRateTableName.empty() )
  {
    m_targetPhaseRateTableName = getName()+"_ConstantPhaseRate_table";
    m_targetPhaseRateTable = createWellTable( m_targetPhaseRateTableName, m_targetPhaseRate );
  }
  else
  {
    FunctionManager & functionManager = FunctionManager::getInstance();
    m_targetPhaseRateTable = &(functionManager.getGroup< TableFunction const >( m_targetPhaseRateTableName ));

    GEOSX_THROW_IF_IF( m_targetPhaseRateTable->getInterpolationMethod() != TableFunction::InterpolationType::Lower,
                    "WellControls '" << getName() << "': The interpolation method for the time-dependent rate table "
                                     << m_targetPhaseRateTable->getName() << " should be TableFunction::InterpolationType::Lower",
                    InputError );
  }

  // 12) Create the time-dependent well status table
  if( m_statusTableName.empty())
  {
    // All well controls without a specified status function will use the same "Open" status function.
    m_statusTableName = GEOSX_FMT( "{0}_OpenStatus_table", dataRepository::keys::wellControls );
    FunctionManager & functionManager = FunctionManager::getInstance();
    m_statusTable = functionManager.getGroupPointer< TableFunction const >( m_statusTableName );
    if( m_statusTable==nullptr )
    {
      m_statusTable = createWellTable( m_statusTableName, 1.0 );
    }
  }
  else
  {
    FunctionManager & functionManager = FunctionManager::getInstance();
    m_statusTable = &(functionManager.getGroup< TableFunction const >( m_statusTableName ));

    GEOSX_THROW_IF_IF( m_statusTable->getInterpolationMethod() != TableFunction::InterpolationType::Lower,
                    "WellControls '" << getName() << "': The interpolation method for the time-dependent rate table "
                                     << m_targetPhaseRateTable->getName() << " should be TableFunction::InterpolationType::Lower",
                    InputError );
  }
}

bool WellControls::isWellOpen( real64 const & currentTime ) const
{
  bool isOpen = true;
  if( isZero( getTargetTotalRate( currentTime ) ) && isZero( getTargetPhaseRate( currentTime ) ) )
  {
    isOpen = false;
  }
  if( m_statusTable->evaluate( &currentTime ) < LvArray::NumericLimits< real64 >::epsilon )
  {
    isOpen = false;
  }
  return isOpen;
}

} //namespace geos
