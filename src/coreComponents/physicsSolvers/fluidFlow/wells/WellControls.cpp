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
 * @file WellControls.cpp
 */

#include "LogLevelsInfo.hpp"
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
  m_currentControl( Control::UNINITIALIZED ), // tjb remove
  m_useSurfaceConditions( 0 ),
  m_surfacePres( -1.0 ),
  m_surfaceTemp( -1.0 ),
  m_isCrossflowEnabled( 1 ),
  m_initialPressureCoefficient( 0.1 ),
  m_rateSign( -1.0 ),
  m_statusTable( nullptr ),
  m_wellOpen( false ),
  m_estimateSolution( 0 ),
  m_constraintSwitch( true ),
  m_currentConstraint( nullptr ),
  m_wellStatus( WellControls::Status::OPEN ),
  m_regionAveragePressure( -1 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::typeString(), &m_type ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Well type. Valid options:\n* " + EnumStrings< Type >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::currentControlString(), &m_currentControl ).
    setDefaultValue( Control::UNINITIALIZED ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Current well control" );



  registerWrapper( viewKeyStruct::useSurfaceConditionsString(), &m_useSurfaceConditions ).
    setDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to specify whether rates are checked at surface or reservoir conditions.\n"
                    "Equal to 1 for surface conditions, and to 0 for reservoir conditions.\n"
                    "See note on referenceReservoirRegion for reservoir condition options" );

  registerWrapper( viewKeyStruct::referenceReservoirRegionString(), &m_referenceReservoirRegion ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setDefaultValue( "" ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of reservoir region used for obtaining average region pressure used in volume rate constraint calculations.\n"
                    "Frequency of pressure update is set in Single/CompositionalMultiPhaseStatistics definition.\n"
                    "Setting cycleFrequency='1' will update the pressure every timestep, note that is a lagged property in constraint properties" );

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


  addLogLevel< logInfo::WellControl >();
}


WellControls::~WellControls()
{}

Group * WellControls::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( GEOS_FMT( "{}: adding {} {}", getName(), childKey, childName ) );
  ////const auto childTypes = { viewKeyStruct::perforationString() };
  //GEOS_ERROR_IF( childKey != viewKeyStruct::perforationString(),
  //               CatalogInterface::unknownTypeError( childKey, getDataContext(), childTypes ) );

  Group * constraint = nullptr;
  if( childKey == viewKeyStruct::minimumBHPConstraintString() )
  {
    MinimumBHPConstraint & bhpConstraint = registerGroup< MinimumBHPConstraint >( childName );
    m_minBHPConstraint =   &bhpConstraint;
    constraint = &bhpConstraint;
  }
  else if( childKey == viewKeyStruct::maximumBHPConstraintString() )
  {
    MaximumBHPConstraint & bhpConstraint = registerGroup< MaximumBHPConstraint >( childName );
    m_maxBHPConstraint =  &bhpConstraint;
    constraint = &bhpConstraint;
  }
  else if( childKey == viewKeyStruct::phaseProductionConstraintString() )
  {
    //PhaseProductionConstraint & phaseConstraint = registerGroup< PhaseProductionConstraint >( childName );
    PhaseProductionConstraint1 & phaseConstraint = registerGroup< PhaseProductionConstraint1 >( childName );
    m_productionRateConstraintList.emplace_back( &phaseConstraint );
    constraint = &phaseConstraint;
  }
  else if( childKey == viewKeyStruct::phaseInjectionConstraintString() )
  {
    //PhaseInjectionConstraint & phaseConstraint = registerGroup< PhaseInjectionConstraint >( childName );
    PhaseInjectionConstraint1 & phaseConstraint = registerGroup< PhaseInjectionConstraint1 >( childName );
    m_injectionRateConstraintList.emplace_back( &phaseConstraint );
    constraint = &phaseConstraint;
  }
  else if( childKey == viewKeyStruct::totalVolProductionConstraintString() )
  {
    TotalVolProductionConstraint & volConstraint = registerGroup< TotalVolProductionConstraint >( childName );
    m_productionRateConstraintList.emplace_back( &volConstraint );
    constraint = &volConstraint;
  }
  else if( childKey == viewKeyStruct::totalVolInjectionConstraintString() )
  {
    TotalVolInjectionConstraint & volConstraint = registerGroup< TotalVolInjectionConstraint >( childName );
    m_injectionRateConstraintList.emplace_back( &volConstraint );
    constraint = &volConstraint;
  }
  else if( childKey == viewKeyStruct::massProductionConstraintString() )
  {
    MassProductionConstraint & massConstraint = registerGroup< MassProductionConstraint >( childName );
    m_productionRateConstraintList.emplace_back( &massConstraint );
    constraint = &massConstraint;

  }
  else if( childKey == viewKeyStruct::massInjectionConstraintString() )
  {
    MassInjectionConstraint & massConstraint = registerGroup< MassInjectionConstraint >( childName );
    m_injectionRateConstraintList.emplace_back( &massConstraint );
    constraint = &massConstraint;
  }
  else if( childKey == viewKeyStruct::liquidProductionConstraintString() )
  {
    LiquidProductionConstraint & liquidConstraint = registerGroup< LiquidProductionConstraint >( childName );
    m_productionRateConstraintList.emplace_back( &liquidConstraint );
    constraint = &liquidConstraint;
  }
  return constraint;
}

void WellControls::expandObjectCatalogs()
{
  //createChild( keys::wellControls, keys::wellControls );
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
  table->setTableCoordinates( timeCoord, { units::Time } );
  table->setTableValues( constantValueArray );
  table->setInterpolationMethod( TableFunction::InterpolationType::Lower );
  return table;
}

}

void WellControls::postInputInitialization()
{

  // 1.c) Set the multiplier for the rates
  if( isProducer() )
  {
    m_rateSign = -1.0;
  }
  else
  {
    m_rateSign = 1.0;
  }

  // 3) check the flag for surface / reservoir conditions
  GEOS_THROW_IF( m_useSurfaceConditions != 0 && m_useSurfaceConditions != 1,
                 getWrapperDataContext( viewKeyStruct::useSurfaceConditionsString() ) << ": The flag to select surface/reservoir conditions must be equal to 0 or 1",
                 InputError );


  //GEOS_THROW_IF( ((m_targetMassRate > 0.0 &&  m_useSurfaceConditions==0)),
  //               "WellControls " << getDataContext() << ": Option only valid if useSurfaceConditions set to 1",
  //               InputError );


  // 8) Make sure that the initial pressure coefficient is positive
  GEOS_THROW_IF( m_initialPressureCoefficient < 0,
                 getWrapperDataContext( viewKeyStruct::initialPressureCoefficientString() ) <<
                 ": This tuning coefficient is negative",
                 InputError );



  // 12) Create the time-dependent well status table
  if( m_statusTableName.empty())
  {
    // All well controls without a specified status function will use the same "Open" status function.
    m_statusTableName = GEOS_FMT( "{0}_OpenStatus_table", dataRepository::keys::wellControls );
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

    GEOS_THROW_IF( m_statusTable->getInterpolationMethod() != TableFunction::InterpolationType::Lower,
                   "WellControls " << getDataContext() << ": The interpolation method for the time-dependent status table "
                                   << m_statusTable->getName() << " should be TableFunction::InterpolationType::Lower",
                   InputError );
  }

}

void WellControls::setWellStatus( real64 const & currentTime, WellControls::Status status )
{
  m_wellStatus = status;
  if( m_wellStatus == WellControls::Status::OPEN )
  {

    if( isZero( getTargetTotalRate( currentTime ) ) && isZero( getTargetPhaseRate( currentTime ) )
        && isZero( getTargetMassRate( currentTime ) ) )
    {
      m_wellStatus =  WellControls::Status::CLOSED;
    }
    if( m_statusTable->evaluate( &currentTime ) < LvArray::NumericLimits< real64 >::epsilon )
    {
      m_wellStatus =  WellControls::Status::CLOSED;
    }
  }
}

bool WellControls::isWellOpen() const
{
  return getWellStatus() == WellControls::Status::OPEN;
}

void WellControls::setWellState( bool open )
{
  m_wellOpen = open;
}

bool WellControls::getWellState() const
{
  return m_wellOpen;
}

void WellControls::setConstraintSwitch( bool constraintSwitch )
{
  m_constraintSwitch = constraintSwitch;
}

bool WellControls::getConstraintSwitch() const
{
  return m_constraintSwitch;
}


void WellControls::setNextDtFromTables( real64 const & currentTime, real64 & nextDt )
{
  if( isProducer() )
  {
    getMinBHPConstraint()->setNextDtFromTables( currentTime, nextDt );
    for( auto const & constraint : m_productionRateConstraintList )
    {
      constraint->setNextDtFromTables( currentTime, nextDt );
    }
  }
  else
  {
    getMaxBHPConstraint()->setNextDtFromTables( currentTime, nextDt );
    for( auto const & constraint : m_injectionRateConstraintList )
    {
      constraint->setNextDtFromTables( currentTime, nextDt );
    }
  }

  WellControls::setNextDtFromTable( m_statusTable, currentTime, nextDt );
}

void WellControls::setNextDtFromTable( TableFunction const * table, real64 const currentTime, real64 & nextDt )
{
  if( table )
  {
    // small epsilon to make sure we land on the other side of table interval and pick up the right rate
    real64 const eps = 1e-6;
    real64 const dtLimit = (table->getCoord( &currentTime, 0, TableFunction::InterpolationType::Upper ) - currentTime) * ( 1.0 + eps );
    if( dtLimit > eps && dtLimit < nextDt )
    {
      nextDt = dtLimit;
    }
  }
}

} //namespace geos
