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

#include "WellControls.hpp"

#include "physicsSolvers/fluidFlow/wells/WellInjectionConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellProductionConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellBHPConstraints.hpp"
#include "physicsSolvers/fluidFlow/wells/WellVolumeRateConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellPhaseVolumeRateConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellMassRateConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellLiquidRateConstraint.hpp"

#include "LogLevelsInfo.hpp"
#include "WellConstants.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionManager.hpp"
#include "mesh/PerforationFields.hpp"
#include "fileIO/Outputs/OutputBase.hpp"
#include "physicsSolvers/fluidFlow/wells/WellFields.hpp"

#include "functions/FunctionManager.hpp"

namespace geos
{

using namespace dataRepository;

WellControls::WellControls( string const & name, Group * const parent )
  : Group( name, parent ),
  m_type( Type::PRODUCER ),
  m_numPhases( 0 ),
  m_numComponents( 0 ),
  m_numDofPerWellElement( 0 ),
  m_numDofPerResElement( 0 ),
  m_isThermal( 0 ),
  m_keepVariablesConstantDuringInitStep( false ),
  m_ratesOutputDir( joinPath( OutputBase::getOutputDirectory(), parent->getName() + "_rates" ) ),
  m_inputControl( Control::UNINITIALIZED ),
  m_currentControl( Control::UNINITIALIZED ),
  m_useSurfaceConditions( 0 ),
  m_surfacePres( -1.0 ),
  m_surfaceTemp( -1.0 ),
  m_isCrossflowEnabled( 1 ),
  m_initialPressureCoefficient( 0.5 ),
  m_currentConstraint( nullptr ),
  m_wellStatus( WellControls::Status::OPEN ),
  m_wellOpen( false ),
  m_statusTable( nullptr ),
  m_regionAveragePressure( -1 ),
  m_estimateSolution( 0 ),
  m_enableIsoThermalEstimator( 0 ),
  /// Nonlinear solver parameters
  m_wellNewtonSolver( groupKeyStruct::wellNewtonSolverString(), this ),
  m_estimatorDoFManager( name ),
  m_dofManagerInitialized( false )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::typeString(), &m_type ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Well type. Valid options:\n* " + EnumStrings< Type >::concat( "\n* " ) );

  this->registerWrapper( viewKeyStruct::writeCSVFlagString(), &m_writeCSV ).
    setApplyDefaultValue( 1 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "When set to 1, write the rates into a CSV file." );

  this->registerWrapper( viewKeyStruct::timeStepFromTablesFlagString(), &m_timeStepFromTables ).
    setApplyDefaultValue( 0 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Choose time step to honor rates/bhp tables time intervals" );

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

  this->registerWrapper( viewKeyStruct::enableIsoThermalEstimatorString(), &m_enableIsoThermalEstimator ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to enable isothermal estimator prior to coupled reservoir and well solve." );

  registerWrapper( viewKeyStruct::statusTableNameString(), &m_statusTableName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the well status table when the status of the well is a time dependent function. \n"
                    "If the status function evaluates to a positive value at the current time, the well will be open otherwise the well will be shut." );

  registerGroup( groupKeyStruct::wellNewtonSolverString(), &m_wellNewtonSolver );

  addLogLevel< logInfo::WellControl >();
}


WellControls::~WellControls()
{}

Group * WellControls::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( GEOS_FMT( "{}: adding {} {}", getName(), childKey, childName ) );
  std::unique_ptr< WellConstraintBase > constraint =
    WellConstraintBase::CatalogInterface::factory( childKey, getDataContext(), childName, this );
  return &registerGroup< WellConstraintBase >( childName, std::move( constraint ) );
}

void WellControls::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from WellConstraintBase here
  for( auto & catalogIter : WellConstraintBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
  // tjbcreateChild( WellNewtonSolver::catalogName(), WellNewtonSolver::catalogName() );
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

void WellControls::registerWellDataOnMesh( WellElementSubRegion & subRegion )
{
  std::string const & regionName = subRegion.getName();
  std::string addrWithMask( regionName );
  std::size_t pos = addrWithMask.find( "UniqueSubRegion" );
  std::string addr = addrWithMask.substr( 0, pos );
  m_targetRegionNames.push_back( addr );

  registerWrapper< real64 >( viewKeyStruct::currentBHPString() );
  registerWrapper< real64 >( viewKeyStruct::currentVolRateString() );

  registerWrapper< array1d< real64 > >( viewKeyStruct::currentPhaseVolRateString() ).
    setSizedFromParent( 0 ).
    reference().resizeDimension< 0 >( m_numPhases );
  registerWrapper< real64 >( viewKeyStruct::massDensityString() );

  registerWrapper< real64 >( viewKeyStruct::currentTotalVolRateString() );
  registerWrapper< real64 >( viewKeyStruct::currentMassRateString() );

  // If estimator is used including thermal effects set during constraint evaluation
  // otherwise they are always included
  if( isThermal() )
  {
    if( m_estimateSolution == 0 )
    {
      enableThermalEffects( true );
    }
  }
}

void WellControls::postInputInitialization()
{
  Group::postInputInitialization();
  // 0) Assign the value of the current well control
  // When the simulation starts from a restart file, we don't want to use the inputControl,
  // because the control may have switched in the simulation that generated the restart
  GEOS_THROW_IF( m_inputControl == Control::UNINITIALIZED,
                 "Input well control cannot be uninitialized",
                 InputError, getWrapperDataContext( viewKeyStruct::inputControlString() ) );

  if( m_currentControl == Control::UNINITIALIZED )
  {
    m_currentControl = m_inputControl;
  }


  // 3) check the flag for surface / reservoir conditions
  GEOS_THROW_IF( m_useSurfaceConditions != 0 && m_useSurfaceConditions != 1,
                 "The flag to select surface/reservoir conditions must be equal to 0 or 1",
                 InputError, getWrapperDataContext( viewKeyStruct::useSurfaceConditionsString() ) );



  // 6.2) Check incoherent information

  // An injector must be controlled by TotalVolRate
  GEOS_THROW_IF( (isInjector() && (m_inputControl == Control::PHASEVOLRATE)),
                 GEOS_FMT( "You have to control an injector with {}",
                           EnumStrings< Control >::toString( Control::TOTALVOLRATE ) ),
                 InputError, getDataContext() );

  // An injector must be controlled by TotalVolRate
  GEOS_THROW_IF( (isProducer() && (m_inputControl == Control::MASSRATE)),
                 GEOS_FMT( "You have to control an injector with {}",
                           EnumStrings< Control >::toString( Control::MASSRATE ) ),
                 InputError, getDataContext() );

  // 8) Make sure that the initial pressure coefficient is positive
  GEOS_THROW_IF( m_initialPressureCoefficient < 0,
                 GEOS_FMT( "{}This tuning coefficient is negative",
                           viewKeyStruct::initialPressureCoefficientString() ),
                 InputError, getWrapperDataContext( viewKeyStruct::initialPressureCoefficientString() ) );



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
                   GEOS_FMT( "The interpolation method for the time-dependent status table {} "
                             "should be TableFunction::InterpolationType::Lower",
                             m_statusTable->getName() ),
                   InputError, getDataContext() );
  }

//##  // 13) Validate constraints
//##  bool const producer = isProducer();
//##  ((void)producer);
//##  mapBase <string,integer,std::false_type> typeCount;
//##  for( auto & catalog : WellConstraintBase::getCatalog())
//##  {
//##    typeCount[catalog.first] = 0;
//##  }
//##  forSubGroups<
//##  >([&](auto const & group)
//##  {
//##    string const catalogName = group.getCatalogName();
//##    ((void)catalogName);
//##//if (typeCount. in )
//##  });
}

void WellControls::postRestartInitialization( )
{}

void WellControls::logConstraint( WellConstraintBase const * constraint, WellElementSubRegion const & region, real64 time, bool isLimiting ) const
{
  bool const needsLog = (constraint != nullptr) && (getLogLevel() > 4) && region.isLocallyOwned();
  if( isLimiting )
  {
    GEOS_LOG_RANK_IF ( needsLog,
                       GEOS_FMT( " Well {}: Limiting Constraint {} - BHP {}, Volume rates {}, Total rate {}, Mass rate {}",
                                 region.getName(), constraint->getName(), constraint->bottomHolePressure(),
                                 constraint->phaseVolumeRates(), constraint->totalVolumeRate(), constraint->massRate()) );
  }
  else
  {
    GEOS_LOG_RANK_IF ( needsLog,
                       GEOS_FMT( " Well {}: Constraint {} - active {}, value {}",
                                 region.getName(), constraint->getName(), constraint->isConstraintActive(),
                                 constraint->getConstraintValue( time ) ) );
  }
}

void WellControls::setWellStatus( real64 const & currentTime, WellControls::Status status )
{
  m_wellStatus = status;
  if( m_wellStatus == WellControls::Status::OPEN )
  {
    bool hasZeroRate = false;

    if( m_statusTable->evaluate( &currentTime ) < LvArray::NumericLimits< real64 >::epsilon )
    {
      hasZeroRate = true;
    }

    if( !hasZeroRate )
    {
      for( auto const * constraint : getRateConstraints() )
      {
        if( isZero( constraint->getConstraintValue( currentTime ) ) )
        {
          hasZeroRate = true;
        }
      }
    }

    if( hasZeroRate )
    {
      m_wellStatus = WellControls::Status::CLOSED;
      m_currentConstraint = nullptr;
    }
  }
}

real64 WellControls::setNextDt( real64 const & currentTime,
                                real64 const & currentDt,
                                WellElementSubRegion & subRegion )
{
  real64 nextDt = currentDt;
  real64 nextDt_perf=nextDt;

  // Find min dt from perf status tables
  PerforationData & perforationData = *subRegion.getPerforationData();
  string_array const & perfStatusTableName = perforationData.getPerfStatusTableName();
  FunctionManager & functionManager = FunctionManager::getInstance();
  // Get dt for local perforations
  for( integer i=0; i<perforationData.size(); i++ )
  {
    TableFunction * tableFunction =  functionManager.getGroupPointer< TableFunction >( perfStatusTableName[i] );
    setNextDtFromTable( tableFunction, currentTime, nextDt_perf );
  }
  nextDt = MpiWrapper::min< real64 >( nextDt_perf );
  // Find min dt including rate and status tables
  real64 const nextDt_orig = nextDt;
  setNextDtFromTables( currentTime, nextDt );
  //if( m_nonlinearSolverParameters.getLogLevel() > 0 && nextDt < nextDt_orig )
  if( getLogLevel() > 0 && nextDt < nextDt_orig )
    GEOS_LOG_RANK_0( GEOS_FMT( "{}: next time step based on tables coordinates = {}", getName(), nextDt ));
  return nextDt;
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
  for( auto const * constraint : getAllConstraints() )
  {
    constraint->setNextDtFromTables( currentTime, nextDt );
  }

  setNextDtFromTable( m_statusTable, currentTime, nextDt );
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
  return getBHPConstraint()->getConstraintValue( targetTime );
}

real64 WellControls::getInjectionTemperature() const
{
  real64 injectionTemperature = 0.0;
  localIndex firstIndex = -1;  // Used to "capture" the first one
  forSubGroupsIndex< InjectionConstraint< MassRateConstraint >,
                     InjectionConstraint< VolumeRateConstraint >,
                     InjectionConstraint< PhaseVolumeRateConstraint >,
                     InjectionConstraint< LiquidRateConstraint > >( [&] ( localIndex index, auto const & constraint )
  {
    if( firstIndex < 0 )
    {
      if( constraint.isConstraintActive() )
      {
        injectionTemperature = constraint.getInjectionTemperature();
        firstIndex = index;
      }
    }
  } );
  return injectionTemperature;
}

arrayView1d< real64 const > WellControls::getInjectionStream() const
{
  arrayView1d< real64 const > injectionStream;
  localIndex firstIndex = -1;  // Used to "capture" the first one
  forSubGroupsIndex< InjectionConstraint< MassRateConstraint >,
                     InjectionConstraint< VolumeRateConstraint >,
                     InjectionConstraint< PhaseVolumeRateConstraint >,
                     InjectionConstraint< LiquidRateConstraint > >( [&] ( localIndex index, auto const & constraint )
  {
    if( firstIndex < 0 )
    {
      if( constraint.isConstraintActive())
      {
        injectionStream = constraint.getInjectionStream();
        firstIndex = index;
      }
    }
  } );
  return injectionStream;
}

integer WellControls::getConstraintPhaseIndex() const
{
  integer phaseIndex = -1;
  // Validation should make sure we are not mixing constraints.
  // Here we assume that we have zero or one or the other but not both
  forSubGroups< ProductionConstraint< PhaseVolumeRateConstraint >,
                InjectionConstraint< PhaseVolumeRateConstraint > >( [&] ( auto & constraint )
  {
    if( constraint.isConstraintActive() )
    {
      phaseIndex = constraint.getPhaseIndex();
    }
  } );
  return phaseIndex;
}

real64 WellControls::getReferenceElevation() const
{
  real64 referenceElevation = 0.0;
  // Validation should make sure we are not mixing constraints.
  // Here we assume that we have zero or one or the other but not both
  forSubGroups< MinimumBHPConstraint,
                MaximumBHPConstraint >( [&] ( auto & constraint )
  {
    referenceElevation = constraint.getReferenceElevation();
  } );
  return referenceElevation;
}

void WellControls::implicitStepSetup( real64 const & time_n,
                                      real64 const & GEOS_UNUSED_PARAM( dt ),
                                      ElementRegionManager & elemManager,
                                      WellElementSubRegion & subRegion )
{

  GEOS_UNUSED_VAR( elemManager );
  // Set perforation status
  setPerforationStatus( time_n, subRegion );
}

void WellControls::setPerforationStatus( real64 const & time_n, WellElementSubRegion & subRegion )
{
  FunctionManager & functionManager = FunctionManager::getInstance();

  // Set perforation status

  PerforationData & perforationData = *subRegion.getPerforationData();
  string_array const & perfStatusTableName = perforationData.getPerfStatusTableName();
  arrayView1d< integer > perfStatus = perforationData.getLocalPerfStatus();
  // for now set to open
  for( integer i=0; i<perforationData.size(); i++ )
  {
    TableFunction * tableFunction =  functionManager.getGroupPointer< TableFunction >( perfStatusTableName[i] );
    perfStatus[i]=PerforationData::PerforationStatus::OPEN;
    if( tableFunction->evaluate( &time_n ) < LvArray::NumericLimits< real64 >::epsilon )
    {
      perfStatus[i]=PerforationData::PerforationStatus::CLOSED;
    }
  }

  array1d< localIndex > const perfWellElemIndex = perforationData.getField< fields::perforation::wellElementIndex >();
  // global index local elements (size == subregion.size)
  arrayView1d< globalIndex const > globalWellElementIndex = subRegion.getGlobalWellElementIndex();

  arrayView1d< integer const > const elemGhostRank  = subRegion.ghostRank();
  array1d< integer > & currentStatus = subRegion.getWellElementStatus();
  // Local elements
  array1d< integer > & localElemStatus = subRegion.getWellLocalElementStatus();

  integer numLocalElements = subRegion.getNumLocalElements();
  array1d< integer > segStatus( numLocalElements );

  // Local perforations
  for( integer j = 0; j < perforationData.size(); j++ )
  {
    localIndex const iwelem = perfWellElemIndex[j];
    if( elemGhostRank[iwelem] < 0 )
    {
      if( perfStatus[j] )
      {
        segStatus[iwelem] +=1;
      }
    }
  }
  // Broadcast segment status so all cores have same well status
  subRegion.setElementStatus( segStatus );
  integer numOpenElements = 0;
  array1d< integer > const & updatedStatus = subRegion.getWellElementStatus();
  for( integer i=0; i<currentStatus.size(); i++ )
  {
    numOpenElements += updatedStatus[i];
  }
  numOpenElements>0 ?  setWellStatus( time_n, WellControls::Status::OPEN ) :  setWellStatus( time_n, WellControls::Status::CLOSED );


  // Set local well element status array
  for( integer i=0; i<subRegion.size(); i++ )
  {
    integer gi = globalWellElementIndex[i];
    localElemStatus[i] = currentStatus[gi];
  }

}

void WellControls::setGravCoef( WellElementSubRegion & subRegion, R1Tensor const & gravVector )
{
  PerforationData & perforationData = *subRegion.getPerforationData();

  real64 const refElev =  getReferenceElevation();

  arrayView2d< real64 const > const wellElemLocation = subRegion.getElementCenter();
  arrayView1d< real64 > const wellElemGravCoef = subRegion.getField< fields::well::gravityCoefficient >();

  arrayView2d< real64 const > const perfLocation = perforationData.getField< fields::perforation::location >();
  arrayView1d< real64 > const perfGravCoef = perforationData.getField< fields::well::gravityCoefficient >();

  forAll< serialPolicy >( perforationData.size(), [=]( localIndex const iperf )
  {
    // precompute the depth of the perforations
    perfGravCoef[iperf] = LvArray::tensorOps::AiBi< 3 >( perfLocation[iperf], gravVector );
  } );

  forAll< serialPolicy >( subRegion.size(), [=]( localIndex const iwelem )
  {
    // precompute the depth of the well elements
    wellElemGravCoef[iwelem] = LvArray::tensorOps::AiBi< 3 >( wellElemLocation[iwelem], gravVector );
  } );

  forSubGroups< MinimumBHPConstraint, MaximumBHPConstraint >( [&]( auto & constraint )
  {
    // set the reference well element where the BHP control is applied
    real64 const refElev1 = constraint.getReferenceElevation();
    constraint.setReferenceGravityCoef( refElev1 * gravVector[2] );
  } );
  // set the reference well element where the BHP control is applied
  setReferenceGravityCoef( refElev * gravVector[2] );       // tjb remove
}

void WellControls::selectWellConstraint( real64 const & time_n,
                                         real64 const & dt,
                                         integer const cycleNumber,
                                         integer const coupledIterationNumber,
                                         DomainPartition & domain,
                                         MeshLevel & meshLevel,
                                         ElementRegionManager & elemManager,
                                         WellElementSubRegion & subRegion,
                                         DofManager const & dofManager )

{

  // Well state estimated from reservoir conditions
  if( isWellOpen() )
  {
    if( !getWellState() )
    {
      setWellState( 1 );

      initializeWell( domain, meshLevel, subRegion, time_n );
    }
  }
  else
  {
    setWellState( 0 );
  }

  bool useEstimator = coupledIterationNumber < estimateSolution();

  if( getWellState())
  {
    if( useEstimator )
    {
      // Estimate well solution prior to coupled solve
      evaluateConstraints( time_n,
                           dt,
                           cycleNumber,
                           coupledIterationNumber,
                           domain,
                           meshLevel,
                           elemManager,
                           subRegion,
                           dofManager );
    }
    else
    {
      // Evaluate well constraints based on current solution
      evaluateConstraints( time_n,
                           subRegion );
    }

    // If a well is opened and then timestep is cut resulting in the well being shut, if the well is opened
    // the well initialization code requires control type to by synced
    integer owner = -1;
    // Only subregion owner evaluates well control and control changes need to be broadcast to all ranks
    if( subRegion.isLocallyOwned() )
    {
      owner = MpiWrapper::commRank( MPI_COMM_GEOS );
    }
    owner = MpiWrapper::max( owner );
    WellControls::Control wellControl = getControl();
    MpiWrapper::broadcast( wellControl, owner );
    setControl( wellControl );
  }
}

bool WellControls::evaluateConstraints( real64 const & time_n,
                                        WellElementSubRegion & subRegion )
{

  // create list of all constraints to process
  stdVector< WellConstraintBase * > constraintList = getAllConstraints();

  // Get current constraint
  WellConstraintBase * limitingConstraint = nullptr;
  for( auto & constraint : constraintList )
  {
    if( constraint->getName() == getCurrentConstraint()->getName())
    {
      limitingConstraint = constraint;
      logConstraint( limitingConstraint, subRegion, time_n, true );
    }
  }
  constraintList.erase( std::find( constraintList.begin(), constraintList.end(), limitingConstraint ) );

  // Check current against other constraints
  for( auto & constraint : constraintList )
  {

    if( limitingConstraint->getName() != constraint->getName())
    {
      if( constraint->checkViolation( *limitingConstraint, time_n ) )
      {
        limitingConstraint = constraint;
        setControl( static_cast< WellControls::Control >(constraint->getControl()) );      // tjb old
        setCurrentConstraint( constraint );
        GEOS_LOG_RANK_IF ( subRegion.isLocallyOwned(),
                           " Well " << subRegion.getName() << " Control switch " << constraint->getName() << " "  << constraint->getConstraintValue( time_n )  );
      }
    }
  }
  logConstraint( limitingConstraint, subRegion, time_n, true );

  return true;
}

bool WellControls::evaluateConstraints( real64 const & time_n,
                                        real64 const & dt,
                                        integer const cycleNumber,
                                        integer const coupledIterationNumber,
                                        DomainPartition & domain,
                                        MeshLevel & mesh,
                                        ElementRegionManager & elemManager,
                                        WellElementSubRegion & subRegion,
                                        DofManager const & dofManager )
{
  // create list of all constraints to solve
  // note that initializeWells sets the initial constraint
  // tjb reactive control schema to allow use to set if needed
  stdVector< WellConstraintBase * > constraintList = getRateConstraints();
  WellConstraintBase * limitingConstraint = getCurrentConstraint();
  if( isProducer() )
  {
    if( limitingConstraint->getControl() != ConstraintTypeId::BHP )
    {
      {   // remove from list and add BHP constraint
        auto it = std::find( constraintList.begin(), constraintList.end(), limitingConstraint );
        if( it != constraintList.end() )
        {
          constraintList.erase( it );
        }
        forSubGroups< MinimumBHPConstraint >( [&] ( auto & constraint )
        {
          constraintList.emplace_back( &constraint );
        } );
        constraintList.insert( constraintList.begin(), limitingConstraint );
      }
    }
    // Solve minimum bhp constraint first
    if( limitingConstraint == nullptr )
    {
      limitingConstraint = constraintList[0];
    }
  }
  else
  {
    if( limitingConstraint->getControl() != ConstraintTypeId::BHP )
    {
      constraintList.emplace_back( getBHPConstraint() );
    }
  }
  if( isoThermalEstimatorEnabled() )
  {
    enableThermalEffects( false );
    solveConstraint ( limitingConstraint, time_n,
                      dt,
                      cycleNumber,
                      coupledIterationNumber,
                      domain,
                      mesh,
                      elemManager,
                      subRegion,
                      dofManager );
    enableThermalEffects( true );
  }
  solveConstraint ( limitingConstraint, time_n,
                    dt,
                    cycleNumber,
                    coupledIterationNumber,
                    domain,
                    mesh,
                    elemManager,
                    subRegion,
                    dofManager );

  for( auto const & constraint : constraintList )
  {
    GEOS_LOG_RANK_IF ( getLogLevel() > 4 && subRegion.isLocallyOwned(),
                       " Well " << subRegion.getName() << " Constraint " << constraint->getName() << " active " << constraint->isConstraintActive() <<
                       " value " << constraint->getConstraintValue( time_n ) );
    if( constraint->isConstraintActive() && constraint->checkViolation( *limitingConstraint, time_n ))
    {
      limitingConstraint=constraint;
      setControl( static_cast< WellControls::Control >(constraint->getControl()) );                       // tjb old
      setCurrentConstraint( limitingConstraint );
      GEOS_LOG_RANK_IF ( getLogLevel() > 4 && subRegion.isLocallyOwned(),
                         " Well " << subRegion.getName() << " New Limiting Constraint " << constraint->getName() << " active " << constraint->isConstraintActive() <<
                         " value " << constraint->getConstraintValue( time_n ) );
      solveConstraint ( constraint, time_n,
                        dt,
                        cycleNumber,
                        coupledIterationNumber,
                        domain,
                        mesh,
                        elemManager,
                        subRegion,
                        dofManager );
    }
  }
  solveConstraint ( limitingConstraint, time_n,
                    dt,
                    cycleNumber,
                    coupledIterationNumber,
                    domain,
                    mesh,
                    elemManager,
                    subRegion,
                    dofManager );

  logConstraint( limitingConstraint, subRegion, time_n, true );

  return true;
}

void WellControls::setupWellDofs( DomainPartition & domain, WellElementRegion & wellElementRegion,
                                  string const & meshBodyName, MeshLevel const & meshLevel )
{
  if( !m_dofManagerInitialized )
  {
    m_dofManagerInitialized=true;
    m_wellNewtonSolver.setupSystem( *this, domain, meshBodyName, meshLevel, wellElementRegion );

  }
}


void WellControls::solveConstraint( WellConstraintBase *constraint,
                                    real64 const & time_n,
                                    real64 const & dt,
                                    integer const cycleNumber,
                                    integer const coupledIterationNumber,
                                    DomainPartition & domain,
                                    MeshLevel & mesh,
                                    ElementRegionManager & elemManager,
                                    WellElementSubRegion & subRegion,
                                    DofManager const & dofManager )
{
  GEOS_UNUSED_VAR( dt );
  GEOS_UNUSED_VAR( cycleNumber );
  GEOS_UNUSED_VAR( domain );
  GEOS_UNUSED_VAR( mesh );
  GEOS_UNUSED_VAR( elemManager );
  GEOS_UNUSED_VAR( dofManager );

  bool useEstimator =   coupledIterationNumber < estimateSolution();
  if( useEstimator )
  {
    if( getLogLevel() > 4 )
    {
      GEOS_LOG_RANK_0( "Well " <<getName() << " Evaluating constraint " << constraint->getName() << " value " << constraint->getConstraintValue( time_n ) << " active " <<
                       constraint->isConstraintActive() );
    }
    if( constraint->isConstraintActive() )
    {
      setControl( static_cast< WellControls::Control >(constraint->getControl()) );    // tjb old
      setCurrentConstraint( constraint );
      // If a well is opened and then timestep is cut resulting in the well being shut, if the well is opened
// the well initialization code requires control type to by synced
      integer owner = -1;
// Only subregion owner evaluates well control and control changes need to be broadcast to all ranks
      if( subRegion.isLocallyOwned() )
      {
        owner = MpiWrapper::commRank( MPI_COMM_GEOS );
      }
      owner = MpiWrapper::max( owner );
      WellControls::Control wellControl = getControl();
      MpiWrapper::broadcast( wellControl, owner );
      setControl( wellControl );

      m_wellNewtonSolver.solveNonlinearSystem( *this, time_n,
                                               dt,
                                               cycleNumber,
                                               domain,
                                               mesh,
                                               elemManager,
                                               subRegion );

      if( getLogLevel() > 4 )
      {
        GEOS_LOG_RANK_0( "Well " <<getName() << " aft solve Constraint rates " << constraint->getName() << " bhp " << constraint->bottomHolePressure() << " phaseVolRate " <<
                         constraint->phaseVolumeRates() << " totalVolRate " << constraint->totalVolumeRate() << " massRate " << constraint->massRate());
      }
    }


  }

}

void
WellControls::assembleSystem( real64 const & time_n,
                              real64 const & dt,
                              integer const cycleNumber,
                              ElementRegionManager & elemManager,
                              WellElementSubRegion & subRegion,
                              DofManager const & dofManager,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )
{
  GEOS_UNUSED_VAR( cycleNumber );
  assembleWellAccumulationTerms( time_n, dt, subRegion, dofManager, localMatrix.toViewConstSizes(), localRhs );

  assembleWellConstraintTerms( time_n, dt, subRegion, dofManager, localMatrix.toViewConstSizes(), localRhs );

  assembleWellPressureRelations( time_n, dt, subRegion, dofManager, localMatrix.toViewConstSizes(), localRhs );
  computeWellPerforationRates( time_n, dt, elemManager, subRegion );
  assembleWellFluxTerms( time_n, dt, subRegion, dofManager, localMatrix.toViewConstSizes(), localRhs );

}
} //namespace geos
