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
  m_inputControl( Control::UNINITIALIZED ),
  m_currentControl( Control::UNINITIALIZED ), // tjb remove
  m_useSurfaceConditions( 0 ),
  m_surfacePres( -1.0 ),
  m_surfaceTemp( -1.0 ),
  m_isCrossflowEnabled( 1 ),
  m_initialPressureCoefficient( 0.1 ),
  m_rateSign( -1.0 ),
  m_statusTable( nullptr ),
  m_wellOpen( false ),
  m_estimateSolution( false ),
  m_currentConstraint( nullptr ),
  m_wellStatus( WellControls::Status::OPEN ),
  m_regionAveragePressure( -1 ),
  m_enableIsoThermalEstimator( 0 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::typeString(), &m_type ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Well type. Valid options:\n* " + EnumStrings< Type >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::currentControlString(), &m_currentControl ).
    setDefaultValue( Control::UNINITIALIZED ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Current well control" );

  registerWrapper( viewKeyStruct::inputControlString(), &m_inputControl ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Well control. Valid options:\n* " + EnumStrings< Control >::concat( "\n* " ) );

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
                    "Frequency of pressure update is set in SinglePhase/CompositionalMultiphaseStatistics definition.\n"
                    "Setting cycleFrequency='1' will update the pressure every timestep, note that is a lagged property in constraint properties"
                    "Note the event associated with the statists task must be entered before the solver event.\n" );

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

  this->registerWrapper( viewKeyStruct::estimateWellSolutionString(), &m_estimateSolution ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to esitmate well solution prior to coupled reservoir and well solve." );

  registerWrapper( viewKeyStruct::enableIsoThermalEstimatorString(), &m_enableIsoThermalEstimator ).
    setDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Estimator configuration option to disable thermal effects on initial well constraint solve and then converge solution with thermal effects enabled: \n"
                    " - If the flag is set to 1, thermal effects are enabled during the initial constraint solve. \n"
                    " - If the flag is set to 0, thermal effects are disabled during the initial constraint solve." );


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
  else if( childKey == viewKeyStruct::productionPhaseVolumeRateConstraintString() )
  {
    ProductionConstraint< PhaseVolumeRateConstraint > & phaseConstraint = registerGroup< ProductionConstraint< PhaseVolumeRateConstraint > >( childName );
    m_productionRateConstraintList.emplace_back( &phaseConstraint );
    constraint = &phaseConstraint;
  }
  else if( childKey == viewKeyStruct::injectionPhaseVolumeRateConstraint() )
  {

    InjectionConstraint< PhaseVolumeRateConstraint > & phaseConstraint = registerGroup< InjectionConstraint< PhaseVolumeRateConstraint > >( childName );
    m_injectionRateConstraintList.emplace_back( &phaseConstraint );
    constraint = &phaseConstraint;
  }
  else if( childKey == viewKeyStruct::productionVolumeRateConstraint() )
  {
    ProductionConstraint< VolumeRateConstraint > & volConstraint = registerGroup< ProductionConstraint< VolumeRateConstraint > >( childName );
    m_productionRateConstraintList.emplace_back( &volConstraint );
    constraint = &volConstraint;
  }
  else if( childKey == viewKeyStruct::injectionVolumeRateConstraint() )
  {
    InjectionConstraint< VolumeRateConstraint > & volConstraint = registerGroup< InjectionConstraint< VolumeRateConstraint > >( childName );
    m_injectionRateConstraintList.emplace_back( &volConstraint );
    constraint = &volConstraint;
  }
  else if( childKey == viewKeyStruct::productionMassRateConstraint() )
  {
    ProductionConstraint< MassRateConstraint > & massConstraint = registerGroup< ProductionConstraint< MassRateConstraint > >( childName );
    m_productionRateConstraintList.emplace_back( &massConstraint );
    constraint = &massConstraint;

  }
  else if( childKey == viewKeyStruct::injectionMassRateConstraint() )
  {
    InjectionConstraint< MassRateConstraint > & massConstraint = registerGroup< InjectionConstraint< MassRateConstraint > >( childName );
    m_injectionRateConstraintList.emplace_back( &massConstraint );
    constraint = &massConstraint;
  }
  else if( childKey == viewKeyStruct::productionLiquidRateConstraint() )
  {
    ProductionConstraint< LiquidRateConstraint > & liquidConstraint = registerGroup< ProductionConstraint< LiquidRateConstraint > >( childName );
    m_productionRateConstraintList.emplace_back( &liquidConstraint );
    constraint = &liquidConstraint;
  }
  return constraint;
}

void WellControls::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from ConstitutiveBase here
  for( auto & catalogIter: WellConstraintBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
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

  // 0) Assign the value of the current well control
  // When the simulation starts from a restart file, we don't want to use the inputControl,
  // because the control may have switched in the simulation that generated the restart
  GEOS_THROW_IF( m_inputControl == Control::UNINITIALIZED,
                 getWrapperDataContext( viewKeyStruct::inputControlString() ) <<
                 ": Input well control cannot be uninitialized",
                 InputError, getWrapperDataContext( viewKeyStruct::inputControlString() ) );

  if( m_currentControl == Control::UNINITIALIZED )
  {
    m_currentControl = m_inputControl;
  }


  // 3) check the flag for surface / reservoir conditions
  GEOS_THROW_IF( m_useSurfaceConditions != 0 && m_useSurfaceConditions != 1,
                 getWrapperDataContext( viewKeyStruct::useSurfaceConditionsString() ) << ": The flag to select surface/reservoir conditions must be equal to 0 or 1",
                 InputError, getWrapperDataContext( viewKeyStruct::useSurfaceConditionsString() ) );

  // tjb add more constraint validation
  // 1) liquid rate - phase names consistent with fluild model
  // 2) at least one bhp and one rate constraint defined
  // 3) constraint type and well type compatibility

  //GEOS_THROW_IF( ((m_targetMassRate > 0.0 &&  m_useSurfaceConditions==0)),
  //               "WellControls " << getDataContext() << ": Option only valid if useSurfaceConditions set to 1",
  //               InputError );


  // 8) Make sure that the initial pressure coefficient is positive
  GEOS_THROW_IF( m_initialPressureCoefficient < 0,
                 getWrapperDataContext( viewKeyStruct::initialPressureCoefficientString() ) <<
                 ": This tuning coefficient is negative",
                 InputError );

  // 6.2) Check incoherent information

  // An injector must be controlled by TotalVolRate
  GEOS_THROW_IF( (isInjector() && (m_inputControl == Control::PHASEVOLRATE)),
                 "WellControls " << getDataContext() << ": You have to control an injector with "
                                 << EnumStrings< Control >::toString( Control::TOTALVOLRATE ),
                 InputError, getDataContext() );

  // An injector must be controlled by TotalVolRate
  GEOS_THROW_IF( (isProducer() && (m_inputControl == Control::MASSRATE)),
                 "WellControls " << getDataContext() << ": You have to control an injector with "
                                 << EnumStrings< Control >::toString( Control::MASSRATE ),
                 InputError, getDataContext() );

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
                   InputError, getDataContext() );
  }

}

void WellControls::setWellStatus( real64 const & currentTime, WellControls::Status status )
{
  m_wellStatus = status;
  if( m_wellStatus == WellControls::Status::OPEN )
  {
    if( isProducer())
    {
      std::vector< WellConstraintBase * >  const constraints =  getProdRateConstraints();
      for( auto const & constraint : constraints )
      {
        if( isZero( constraint->getConstraintValue( currentTime ) ) )
        {
          m_wellStatus =  WellControls::Status::CLOSED;
          break;
        }
      }
    }
    else
    {
      std::vector< WellConstraintBase * >  const constraints =  getInjRateConstraints();
      for( auto const & constraint : constraints )
      {
        if( isZero( constraint->getConstraintValue( currentTime ) ) )
        {
          m_wellStatus =  WellControls::Status::CLOSED;
          break;
        }
      }
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

void WellControls::setNextDtFromTables( real64 const & currentTime, real64 & nextDt )
{
  if( isProducer() )
  {
    if( getMinBHPConstraint() != nullptr )
    {
      getMinBHPConstraint()->setNextDtFromTables( currentTime, nextDt );
    }
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

real64 WellControls::getTargetBHP( real64 const & targetTime ) const
{
  if( isProducer())
  {
    return m_minBHPConstraint->getConstraintValue( targetTime );
  }
  return m_maxBHPConstraint->getConstraintValue( targetTime );
}


real64 WellControls::getInjectionTemperature() const
{
  real64 injectionTemperature = 0.0;
  this->forInjectionConstraints< InjectionConstraint< PhaseVolumeRateConstraint >, InjectionConstraint< VolumeRateConstraint >, InjectionConstraint< MassRateConstraint > >( [&] ( auto & constraint )
  {
    if( constraint.isConstraintActive())
    {
      injectionTemperature =  constraint.getInjectionTemperature();
      return;
    }
  } );
  return injectionTemperature;
}


arrayView1d< real64 const > WellControls::getInjectionStream() const
{
  arrayView1d< real64 const > injectionStream;
  forInjectionConstraints< InjectionConstraint< PhaseVolumeRateConstraint >, InjectionConstraint< VolumeRateConstraint >, InjectionConstraint< MassRateConstraint > >( [&] ( auto & constraint )
  {
    if( constraint.isConstraintActive() )
    {
      injectionStream = constraint.getInjectionStream();
      return;
    }
  } );

  return injectionStream;
}

integer WellControls::getConstraintPhaseIndex() const
{
  integer phaseIndex = -1;

  if( isProducer() )
  {
    forProductionConstraints< ProductionConstraint< PhaseVolumeRateConstraint > >( [&] ( auto & constraint )
    {
      if( constraint.isConstraintActive() )
      {
        phaseIndex = constraint.getPhaseIndex();
      }
    } );
  }
  else
  {
    forInjectionConstraints< InjectionConstraint< PhaseVolumeRateConstraint > >( [&] ( auto & constraint )
    {
      if( constraint.isConstraintActive() )
      {
        phaseIndex = constraint.getPhaseIndex();
      }
    } );
  }

  return phaseIndex;
}

real64 WellControls::getReferenceElevation() const
{
  if( isProducer () )
  {
    return getMinBHPConstraint()->getReferenceElevation();
  }
  return getMaxBHPConstraint()->getReferenceElevation();
}
} //namespace geos
