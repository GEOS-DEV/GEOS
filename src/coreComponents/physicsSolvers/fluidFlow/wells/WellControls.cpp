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
#include "ProdPipeFlowTableFunction.hpp"
#include "InjPipeFlowTableFunction.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionManager.hpp"
#include "mesh/PerforationFields.hpp"
#include "fileIO/Outputs/OutputBase.hpp"
#include "physicsSolvers/fluidFlow/wells/WellFields.hpp"

#include "functions/FunctionManager.hpp"
namespace geos
{
void printlmat( std::string loc, CRSMatrixView< real64, globalIndex const > const & localMatrix, arrayView1d< real64 > const & localRhs )
{
  std::cout << "Local matrix at " << loc << ":" << std::endl;
  // Print matrix information using proper CRS matrix access methods
  std::cout << "LocalMatrix info:" << std::endl;
  std::cout << "Number of rows: " << localMatrix.numRows() << std::endl;
  std::cout << "Number of columns: " << localMatrix.numColumns() << std::endl;
  // Print first few rows of the matrix
  for( localIndex row = 0; row < std::min( localIndex( 5 ), localMatrix.numRows() ); ++row )
  {
    std::cout << "Row " << row << " (nnz=" << localMatrix.numNonZeros( row ) << "): ";
    if( localMatrix.numNonZeros( row ) > 0 )
    {
      auto entries = localMatrix.getEntries( row );
      auto columns = localMatrix.getColumns( row );
      for( localIndex j = 0; j < std::min( localIndex( 5 ), localMatrix.numNonZeros( row ) ); ++j )
      {
        std::cout << "[" << columns[j] << "]=" << entries[j] << " ";
      }
      if( localMatrix.numNonZeros( row ) > 5 )
      {
        std::cout << " .. ";
      }
    }

    std::cout << localRhs[row] << std::endl;
  }
}
using namespace dataRepository;

WellControls::WellControls( string const & name, Group * const parent )
  : Group( name, parent ),
  m_type( Type::PRODUCER ),
  m_useMass( 0 ),
  m_numPhases( 0 ),
  m_numComponents( 0 ),
  m_numDofPerWellElement( 0 ),
  m_numDofPerResElement( 0 ),
  m_isThermal( 0 ),
  m_keepVariablesConstantDuringInitStep( false ),
  //m_ratesOutputDir( joinPath( OutputBase::getOutputDirectory(), "_rates" ) ),
  m_ratesOutputDir( joinPath( OutputBase::getOutputDirectory(), parent->getName() + "_rates" ) ),
  m_inputControl( Control::UNINITIALIZED ),
  m_currentControl( Control::UNINITIALIZED ),
  m_useSurfaceConditions( 0 ),
  m_surfacePres( -1.0 ),
  m_surfaceTemp( -1.0 ),
  m_isCrossflowEnabled( 1 ),
  m_initialPressureCoefficient( 0.5 ),
  m_currentConstraint( nullptr ),
  m_minBHPConstraint( nullptr ),
  m_maxBHPConstraint( nullptr ),
  m_minWHPConstraint( nullptr ),
  m_maxWHPConstraint( nullptr ),
  m_wellStatus( WellControls::Status::OPEN ),
  m_wellOpen( false ),
  m_statusTable( nullptr ),
  m_regionAveragePressure( -1 ),
  m_estimateSolution( 0 ),
  m_enableIsoThermalEstimator( 0 ),
  /// Nonlinear solver parameters
  m_wellNewtonSolver( viewKeyStruct::wellNewtonSolverString(), this ),
  m_estimatorDoFManager( name ),
  m_dofManagerInitialized( false ),
  m_writeSegDebug( 0 ),

  m_numTimesteps( 0 ),
  m_wellDebugInit( false )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
  m_minWHPConstraint=nullptr;
  registerWrapper( viewKeyStruct::typeString(), &m_type ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Well type. Valid options:\n* " + EnumStrings< Type >::concat( "\n* " ) );

  this->registerWrapper( viewKeyStruct::writeCSVFlagString(), &m_writeCSV ).
    setApplyDefaultValue( 1 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "When set to 1, write the rates into a CSV file." );

  this->registerWrapper( viewKeyStruct::writeSegDebugFlagString(), &m_writeSegDebug ).
    setApplyDefaultValue( 0 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "When set to 1, write the segment debug information into a CSV file." );

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


  registerGroup( viewKeyStruct::wellNewtonSolverString(), &m_wellNewtonSolver );

  addLogLevel< logInfo::WellControl >();
}


WellControls::~WellControls()
{}

Group * WellControls::createChild( string const & childKey, string const & childName )
{

  Group * child = nullptr;
  if( childKey == viewKeyStruct::minimumBHPConstraintString() )
  {
    MinimumBHPConstraint & bhpConstraint = registerGroup< MinimumBHPConstraint >( childName );
    m_minBHPConstraint =   &bhpConstraint;
    child = &bhpConstraint;
  }
  else if( childKey == viewKeyStruct::maximumBHPConstraintString() )
  {
    MaximumBHPConstraint & bhpConstraint = registerGroup< MaximumBHPConstraint >( childName );
    m_maxBHPConstraint =  &bhpConstraint;
    child = &bhpConstraint;
  }
  else if( childKey == viewKeyStruct::minimumWHPConstraintString() )
  {
    MinimumWHPConstraint & whpConstraint = registerGroup< MinimumWHPConstraint >( childName );
    m_minWHPConstraint =  &whpConstraint;
    child = &whpConstraint;
  }
  else if( childKey == viewKeyStruct::maximumWHPConstraintString() )
  {
    MaximumWHPConstraint & whpConstraint = registerGroup< MaximumWHPConstraint >( childName );
    m_maxWHPConstraint =  &whpConstraint;
    child = &whpConstraint;
  }
  else if( childKey == viewKeyStruct::productionPhaseVolumeRateConstraintString() )
  {
    ProductionConstraint< PhaseVolumeRateConstraint > & phaseConstraint = registerGroup< ProductionConstraint< PhaseVolumeRateConstraint > >( childName );
    m_productionRateConstraintList.emplace_back( &phaseConstraint );
    child = &phaseConstraint;
  }
  else if( childKey == viewKeyStruct::injectionPhaseVolumeRateConstraint() )
  {
    InjectionConstraint< PhaseVolumeRateConstraint > & phaseConstraint = registerGroup< InjectionConstraint< PhaseVolumeRateConstraint > >( childName );
    m_injectionRateConstraintList.emplace_back( &phaseConstraint );
    child = &phaseConstraint;
  }
  else if( childKey == viewKeyStruct::productionVolumeRateConstraint() )
  {
    ProductionConstraint< VolumeRateConstraint > & volConstraint = registerGroup< ProductionConstraint< VolumeRateConstraint > >( childName );
    m_productionRateConstraintList.emplace_back( &volConstraint );
    child = &volConstraint;
  }
  else if( childKey == viewKeyStruct::injectionVolumeRateConstraint() )
  {
    InjectionConstraint< VolumeRateConstraint > & volConstraint = registerGroup< InjectionConstraint< VolumeRateConstraint > >( childName );
    m_injectionRateConstraintList.emplace_back( &volConstraint );
    child = &volConstraint;
  }
  else if( childKey == viewKeyStruct::productionMassRateConstraint() )
  {
    ProductionConstraint< MassRateConstraint > & massConstraint = registerGroup< ProductionConstraint< MassRateConstraint > >( childName );
    m_productionRateConstraintList.emplace_back( &massConstraint );
    child = &massConstraint;

  }
  else if( childKey == viewKeyStruct::injectionMassRateConstraint() )
  {
    InjectionConstraint< MassRateConstraint > & massConstraint = registerGroup< InjectionConstraint< MassRateConstraint > >( childName );
    m_injectionRateConstraintList.emplace_back( &massConstraint );
    child = &massConstraint;
  }
  else if( childKey == viewKeyStruct::productionLiquidRateConstraint() )
  {
    ProductionConstraint< LiquidRateConstraint > & liquidConstraint = registerGroup< ProductionConstraint< LiquidRateConstraint > >( childName );
    m_productionRateConstraintList.emplace_back( &liquidConstraint );
    child = &liquidConstraint;
  }

  return child;
}

void WellControls::createMinBHPConstraintForWHP()
{
  // Create constraint and set local pointer
  MinimumBHPConstraint & bhpConstraint = registerGroup< MinimumBHPConstraint >( m_minWHPConstraint->getName()+std::string( "MinimumBHPConstraint" ));
  m_minBHPConstraintForWHP =    &bhpConstraint;
  // Set properties from the original minBHP constraint
  m_minBHPConstraintForWHP->setReferenceElevation( m_minBHPConstraint->getReferenceElevation() );
  m_minBHPConstraintForWHP->setReferenceGravityCoef ( m_minBHPConstraint->getReferenceGravityCoef() );
  // Set to inactive. WHP estimator solve will set status
  m_minBHPConstraintForWHP->setConstraintActive( false );
}
void WellControls::createMaxBHPConstraintForWHP()
{
  // Create constraint and set local pointer
  MaximumBHPConstraint & bhpConstraint = registerGroup< MaximumBHPConstraint >( m_maxWHPConstraint->getName()+std::string( "MaximumBHPConstraint" ) );
  m_maxBHPConstraintForWHP =    &bhpConstraint;
  // Set properties from the original maxBHP constraint
  m_maxBHPConstraintForWHP->setReferenceElevation( m_maxBHPConstraint->getReferenceElevation() );
  m_maxBHPConstraintForWHP->setReferenceGravityCoef ( m_maxBHPConstraint->getReferenceGravityCoef() );
  // Set to inactive. WHP estimator solve will set status
  m_maxBHPConstraintForWHP->setConstraintActive( false );
}
void WellControls::createMaxLiquidConstraintForWHP()
{
  // Create constraint and set local pointer
  ProductionConstraint< LiquidRateConstraint > & liquidConstraint =
    registerGroup< ProductionConstraint< LiquidRateConstraint > >( m_minWHPConstraint->getName()+std::string( "LiquidProductionConstraint" ));
  m_maxLiquidConstraintForWHP =  &liquidConstraint;
  // Set properties from VFP table
  FunctionManager & functionManager = FunctionManager::getInstance();
  const ProdPipeFlowTableFunction & m_flowTable =  functionManager.getGroup< ProdPipeFlowTableFunction const >( m_minWHPConstraint->getFlowTableName());
  string_array ratePhases = m_flowTable.getRatePhases();
  m_maxLiquidConstraintForWHP->setPhaseNames( ratePhases );
  m_maxLiquidConstraintForWHP->validateLiquidType( getMultiFluidSeparator());
  // WHP estimator solve will set status
  m_maxLiquidConstraintForWHP->setConstraintActive( false );
}

void WellControls::createMaxVolumeInjConstraintForWHP()
{

  // Set properties from VFP table
  FunctionManager & functionManager = FunctionManager::getInstance();
  const InjPipeFlowTableFunction & m_flowTable =  functionManager.getGroup< InjPipeFlowTableFunction const >( m_maxWHPConstraint->getFlowTableName());
  std::string const & rateType = m_flowTable.getRateType();
  // Create constraint and set local pointer
  InjectionConstraint< PhaseVolumeRateConstraint > & volumeConstraint =
    registerGroup< InjectionConstraint< PhaseVolumeRateConstraint > >( m_maxWHPConstraint->getName()+std::string( "VolumeInjectionConstraint" ) );
  m_maxPhaseVolumeConstraintForWHP =  &volumeConstraint;

  bool foundMatchingConstraint = false;
  forSubGroups< InjectionConstraint< PhaseVolumeRateConstraint >
                >( [&]( auto & constraint )
  {
    std::cout << "check phase for whp constraint " << constraint.getName() << " rate type " << constraint.getPhaseName() << std::endl;
    if( constraint.isConstraintActive() && constraint.getPhaseName() ==  rateType )
    {
      foundMatchingConstraint = true;
      m_maxPhaseVolumeConstraintForWHP->setPhaseName( rateType );
      m_maxPhaseVolumeConstraintForWHP->setPhaseIndex( constraint.getPhaseIndex() );
      m_maxPhaseVolumeConstraintForWHP->setInjectionStream( constraint.getInjectionStream() );
      m_maxPhaseVolumeConstraintForWHP->setInjectionTemperature( constraint.getInjectionTemperature() );
    }
  } );
  GEOS_THROW_IF( !foundMatchingConstraint, "No active injection phase volume constraint with matching phase found for max WHP constraint " << getMaxWHPConstraint()->getName(),
                 InputError, getDataContext() );

  m_maxPhaseVolumeConstraintForWHP->setConstraintActive( false );
}

void WellControls::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from ConstitutiveBase here
  for( auto & catalogIter: WellConstraintBase::getCatalog())
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
  // TJB if( hasMinimumWHPConstraint() )
  registerWrapper< real64 >( viewKeyStruct::currentWHPString() );
  registerWrapper< real64 >( viewKeyStruct::currentVolRateString() );

  registerWrapper< array1d< real64 > >( viewKeyStruct::currentPhaseVolRateString() ).
    setSizedFromParent( 0 ).
    reference().resizeDimension< 0 >( m_numPhases );
  registerWrapper< real64 >( viewKeyStruct::massDensityString() );

  registerWrapper< real64 >( viewKeyStruct::currentTotalVolRateString() );
  registerWrapper< real64 >( viewKeyStruct::currentMassRateString() );
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
  //GEOS_THROW_IF( (isInjector() && (m_inputControl == Control::PHASEVOLRATE)),
  //               "You have to control an injector with "
  //               << EnumStrings< Control >::toString( Control::TOTALVOLRATE ),
  //               InputError, getDataContext() );

  // An injector must be controlled by TotalVolRate
  GEOS_THROW_IF( (isProducer() && (m_inputControl == Control::MASSRATE)),
                 "You have to control an injector with "
                 << EnumStrings< Control >::toString( Control::MASSRATE ),
                 InputError, getDataContext() );

  // 8) Make sure that the initial pressure coefficient is positive
  GEOS_THROW_IF( m_initialPressureCoefficient < 0,
                 getWrapperDataContext( viewKeyStruct::initialPressureCoefficientString() ) <<
                 ": This tuning coefficient is negative",
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
                   "The interpolation method for the time-dependent status table "
                   << m_statusTable->getName() << " should be TableFunction::InterpolationType::Lower",
                   InputError, getDataContext() );
  }

}
void WellControls::postRestartInitialization( )
{}
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
          m_currentConstraint=nullptr;
          break;
        }
      }
    }
    else
    {
      std::vector< WellConstraintBase * >  const constraints =  getInjRateConstraints();
      for( auto const & constraint : constraints )
      {
        std::cout << "Checking injection constraint " << constraint->getName() << " with value " << constraint->getConstraintValue( currentTime ) << std::endl;
        if( isZero( constraint->getConstraintValue( currentTime ) ) )
        {
          m_wellStatus =  WellControls::Status::CLOSED;
          m_currentConstraint=nullptr;
          break;
        }
      }
    }

    if( m_statusTable->evaluate( &currentTime ) < LvArray::NumericLimits< real64 >::epsilon )
    {
      m_wellStatus =  WellControls::Status::CLOSED;
      m_currentConstraint=nullptr;
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

  forSubGroups< BHPConstraint >( [&]( auto & constraint )
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
  if( isWellOpen()  )
  {
    if( !getWellState()   )
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
      std::cout << "Estimating well solution for well " << subRegion.getName() << " at time " << time_n << std::endl;
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
      std::cout << "Evaluating well constraints for well " << subRegion.getName() << " at time " << time_n << std::endl;
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
  std::vector< WellConstraintBase * > constraintList;
  if( isProducer() )
  {
    constraintList = getProdRateConstraints();
    //if constraints arent updated with estimator and WHP is binding dont check allow constraint to be switch during remainder of timestep
    if( hasMinimumWHPConstraint()  )
    {
      MinimumBHPConstraint * minBHPForWHP =  getMinimumBHPConstraintForWHP();
      if( minBHPForWHP != nullptr && minBHPForWHP->isConstraintActive())
      {
        std::cout << "we not active " << subRegion.getName() << " Constraint " << minBHPForWHP->getName() << " active " << minBHPForWHP->isConstraintActive() <<
          " value " << minBHPForWHP->getConstraintValue( time_n ) << std::endl;
        constraintList.insert( constraintList.begin(), minBHPForWHP );

      }
      else
      {
        ProductionConstraint< LiquidRateConstraint > * maxLiqForWHP =  getMaxLiquidConstraintForWHP();
        if( maxLiqForWHP != nullptr && maxLiqForWHP->isConstraintActive())
        {
          std::cout << "we  not active " << subRegion.getName() << " Constraint " << maxLiqForWHP->getName() << " active " << maxLiqForWHP->isConstraintActive() <<
            " value " << maxLiqForWHP->getConstraintValue( time_n ) << std::endl;
          constraintList.insert( constraintList.begin(), maxLiqForWHP );

        }
        else
        {
          // Solve minimum bhp constraint first
          if( getMinBHPConstraint()->isConstraintActive() )
          {
            std::cout << "we  not active " << subRegion.getName() << " Constraint add minbp " << std::endl;
            constraintList.insert( constraintList.begin(), getMinBHPConstraint() );
          }
        }
      }
    }
    else
    {
      std::cout << "we  not active " << subRegion.getName() << " Constraint " << getMinBHPConstraint()->getName() << " active " << getMinBHPConstraint()->isConstraintActive() <<
        " value " << getMinBHPConstraint()->getConstraintValue( time_n ) << std::endl;
      constraintList.insert( constraintList.begin(), getMinBHPConstraint() );
    }
  }
  else
  {
    constraintList = getInjRateConstraints();
    if( hasMaximumWHPConstraint()  )
    {
      MaximumBHPConstraint * maxBHPForWHP =  getMaximumBHPConstraintForWHP();
      if( maxBHPForWHP != nullptr && maxBHPForWHP->isConstraintActive())
      {
        std::cout << "we not active " << subRegion.getName() << " Constraint " << maxBHPForWHP->getName() << " active " << maxBHPForWHP->isConstraintActive() <<
          " value " << maxBHPForWHP->getConstraintValue( time_n ) << std::endl;
        constraintList.insert( constraintList.begin(), maxBHPForWHP );

      }
      else
      {
        InjectionConstraint< PhaseVolumeRateConstraint > * maxVolForWHP =  getMaxPhaseVolumeConstraintForWHP();
        if( maxVolForWHP != nullptr && maxVolForWHP->isConstraintActive())
        {
          std::cout << "we  not active " << subRegion.getName() << " Constraint " << maxVolForWHP->getName() << " active " << maxVolForWHP->isConstraintActive() <<
            " value " << maxVolForWHP->getConstraintValue( time_n ) << std::endl;
          constraintList.insert( constraintList.begin(), maxVolForWHP );

        }
        else
        {
          // Solve maximum bhp constraint first
          if( getMaxBHPConstraint()->isConstraintActive() )
          {
            std::cout << "we  not active " << subRegion.getName() << " Constraint add maxbp " << std::endl;
            constraintList.insert( constraintList.begin(), getMaxBHPConstraint() );
          }
        }
      }
    }
    else
    {
      // Solve maximum bhp constraint first;
      constraintList.insert( constraintList.begin(), getMaxBHPConstraint() );
    }
  }


// Get current constraint
  WellConstraintBase *  limitingConstraint = nullptr;
  for( auto & constraint : constraintList )
  {
    if( constraint->getName() == getCurrentConstraint()->getName())
    {
      limitingConstraint =  constraint;
      GEOS_LOG_RANK_IF ( getLogLevel() > 4 && subRegion.isLocallyOwned(),
                         " Well " << subRegion.getName() << " Limiting Constraint " << limitingConstraint->getName() << " "  << limitingConstraint->bottomHolePressure() << " " <<
                         limitingConstraint->phaseVolumeRates() << " " <<
                         limitingConstraint->totalVolumeRate() << " " << limitingConstraint->massRate());
    }
  }
// Check current against other constraints
  std::cout << "Current constraint for well " << subRegion.getName() << " is " << limitingConstraint->getName() << std::endl;
  constraintList.erase( std::remove( constraintList.begin(), constraintList.end(), limitingConstraint ), constraintList.end());
  std::vector< int > constraintChecked( constraintList.size(), 0 );
  for( int i = 0; i < static_cast< int >(constraintList.size()); ++i )
  {
    auto & constraint = constraintList[i];
    if( limitingConstraint->getName() != constraint->getName())
    {
      if( !constraintChecked[i] && constraint->checkViolation( *limitingConstraint, time_n ) )
      {
        limitingConstraint = constraint;
        setControl( static_cast< WellControls::Control >(constraint->getControl()) );      // tjb old
        setCurrentConstraint( constraint );
        GEOS_LOG_RANK_IF ( subRegion.isLocallyOwned(),
                           " Well " << subRegion.getName() << " Control switch " << constraint->getName() << " "  << constraint->getConstraintValue( time_n )  );
      }
    }
    constraintChecked[i]=1;
  }
  GEOS_LOG_RANK_IF ( getLogLevel() > 4 && subRegion.isLocallyOwned(),
                     " Well " << subRegion.getName() << " Limiting Constraint " << limitingConstraint->getName() << " "  << limitingConstraint->bottomHolePressure() << " " <<
                     limitingConstraint->phaseVolumeRates() << " " <<
                     limitingConstraint->totalVolumeRate() << " " << limitingConstraint->massRate());

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
  std::vector< WellConstraintBase * >  constraintList;
  WellConstraintBase *  limitingConstraint = getCurrentConstraint();
  if( isProducer() )
  {
    constraintList = getProdRateConstraints();

    if( limitingConstraint->getControl() != ConstraintTypeId::BHP )
    {
      {   // set BHP constraint to be first constraint evaluated
        if( getMinBHPConstraint()->isConstraintActive() )
        {
          constraintList.insert( constraintList.begin(), getMinBHPConstraint() );
        }
      }
    }
    // Solve minimum bhp constraint first
    if( false && getMinBHPConstraint()->isConstraintActive() )
    {
      // this is related to WHP option which introduces a new BHP constraint
      limitingConstraint = getMinBHPConstraint();
    }
    else if( limitingConstraint == nullptr )
    {
      limitingConstraint = constraintList[0];
    }
  }
  else
  {
    constraintList = getInjRateConstraints();
    if( limitingConstraint->getControl() != ConstraintTypeId::BHP )
    {
      {   // set BHP constraint to be first constraint evaluated
        if( getMaxBHPConstraint()->isConstraintActive() )
        {
          constraintList.insert( constraintList.begin(), getMaxBHPConstraint() );
        }
      }
    }
    // Solve minimum bhp constraint first
    if( false && getMinBHPConstraint()->isConstraintActive() )
    {
      // this is related to WHP option which introduces a new BHP constraint
      limitingConstraint = getMinBHPConstraint();
    }
    else if( limitingConstraint == nullptr )
    {
      limitingConstraint = constraintList[0];
    }

  }
  if( isoThermalEstimatorEnabled() )
  {
    std::cout << "Solving limiting constraint " << limitingConstraint->getName() << " without thermal effects for well " << subRegion.getName() << std::endl;

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

  constraintList.erase( std::remove( constraintList.begin(), constraintList.end(), limitingConstraint ), constraintList.end());


  std::vector< int > constraintChecked( constraintList.size(), 0 );
  for( int i = 0; i < static_cast< int >(constraintList.size()); ++i )
  {
    auto & constraint = constraintList[i];
    GEOS_LOG_RANK_IF ( getLogLevel() > 4 && subRegion.isLocallyOwned(),
                       " Well " << subRegion.getName() << " Constraint " << constraint->getName() << " active " << constraint->isConstraintActive() <<
                       " value " << constraint->getConstraintValue( time_n ) );
    if( limitingConstraint->getName() != constraint->getName())
    {
      if( !constraintChecked[i] &&constraint->isConstraintActive()  &&  constraint->checkViolation( *limitingConstraint, time_n ))
      {
        limitingConstraint=constraint;
        setControl( static_cast< WellControls::Control >(constraint->getControl()) );                     // tjb old
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
        // tjb. this is likely not needed. set in update state
        constraint->setBHP ( getReference< real64 >( viewKeyStruct::currentBHPString() ));
        if( isCompositional())
        {
          constraint->setPhaseVolumeRates ( getReference< array1d< real64 > >(
                                              viewKeyStruct::currentPhaseVolRateString() ) );
          constraint->setTotalVolumeRate ( getReference< real64 >(
                                             viewKeyStruct::currentTotalVolRateString() ));
          constraint->setMassRate( getReference< real64 >( viewKeyStruct::currentMassRateString() ));
        }
        else
        {
          constraint->setBHP ( getReference< real64 >( viewKeyStruct::currentBHPString() ));
          constraint->setTotalVolumeRate ( getReference< real64 >(
                                             viewKeyStruct::currentVolRateString() ));
        }

      }
    }
    constraintChecked[i]=1;
  }
#if 1
  solveConstraint ( limitingConstraint, time_n,
                    dt,
                    cycleNumber,
                    coupledIterationNumber,
                    domain,
                    mesh,
                    elemManager,
                    subRegion,
                    dofManager );
#endif
  GEOS_LOG_RANK_IF ( getLogLevel() > 4 && subRegion.isLocallyOwned(),
                     " Well " << subRegion.getName() << " Limiting Constraint " << limitingConstraint->getName() << " "  << limitingConstraint->bottomHolePressure() << " " << limitingConstraint->phaseVolumeRates() << " " <<
                     limitingConstraint->totalVolumeRate() << " " << limitingConstraint->massRate());

  if( hasMinimumWHPConstraint()     )
  {
    bool whpLimiting= solveMinWHPConstraint ( time_n,
                                              dt,
                                              cycleNumber,
                                              coupledIterationNumber,
                                              domain,
                                              mesh,
                                              elemManager,
                                              subRegion );

    if( whpLimiting )
    {
      // WHP option can use different constraint as limiting constraint, so need to update the rates for the current limiting constraint
      limitingConstraint= getCurrentConstraint();

      GEOS_LOG_RANK_IF ( getLogLevel() > 4 && subRegion.isLocallyOwned(),
                         " Well " << subRegion.getName() << " Limiting Constraint " << limitingConstraint->getName() << " "  << limitingConstraint->bottomHolePressure() << " " << limitingConstraint->phaseVolumeRates() << " " <<
                         limitingConstraint->totalVolumeRate() << " " << limitingConstraint->massRate());

    }
  }
  else if( hasMaximumWHPConstraint() )
  {
    bool whpLimiting= solveMaxWHPConstraint ( time_n,
                                              dt,
                                              cycleNumber,
                                              coupledIterationNumber,
                                              domain,
                                              mesh,
                                              elemManager,
                                              subRegion );

    if( whpLimiting )
    {
      // WHP option can use different constraint as limiting constraint, so need to update the rates for the current limiting constraint
      limitingConstraint= getCurrentConstraint();

      GEOS_LOG_RANK_IF ( getLogLevel() > 4 && subRegion.isLocallyOwned(),
                         " Well " << subRegion.getName() << " Limiting Constraint " << limitingConstraint->getName() << " "  << limitingConstraint->bottomHolePressure() << " " << limitingConstraint->phaseVolumeRates() << " " <<
                         limitingConstraint->totalVolumeRate() << " " << limitingConstraint->massRate());

    }
  }
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
  //printlmat( "Well subregion matrix after accumulation", localMatrix ,localRhs);
  assembleWellConstraintTerms( time_n, dt, subRegion, dofManager, localMatrix.toViewConstSizes(), localRhs );
  //printlmat( "Well subregion matrix after constraint terms", localMatrix ,localRhs);
  assembleWellPressureRelations( time_n, dt, subRegion, dofManager, localMatrix.toViewConstSizes(), localRhs );
  //printlmat( "Well subregion matrix after pressure relations", localMatrix ,localRhs);
  computeWellPerforationRates( time_n, dt, elemManager, subRegion );
  assembleWellFluxTerms( time_n, dt, subRegion, dofManager, localMatrix.toViewConstSizes(), localRhs );
  //printlmat( "Well subregion matrix after flux terms", localMatrix ,localRhs);

}
} //namespace geos
