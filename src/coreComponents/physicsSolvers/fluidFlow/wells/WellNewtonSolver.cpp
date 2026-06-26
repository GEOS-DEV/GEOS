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

#include "WellNewtonSolver.hpp"

#include "common/MpiWrapper.hpp"
#include "codingUtilities/RTTypes.hpp"
#include "common/format/EnumStrings.hpp"
#include "dataRepository/Group.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "common/format/LogPart.hpp"
#include "common/TimingMacros.hpp"
#include "linearAlgebra/solvers/KrylovSolver.hpp"
#include "mesh/DomainPartition.hpp"
#include "math/interpolation/Interpolation.hpp"
#include "common/Timer.hpp"
#include "common/Units.hpp"


namespace geos
{

using namespace dataRepository;

WellNewtonSolver::WellNewtonSolver( string const & name,
                                    Group * const parent )
  :
  dataRepository::Group( name, parent ),

  m_dofManager( name ),
  m_usePhysicsScaling( 1 ),
  m_linearSolverParameters( groupKeyStruct::linearSolverParametersString(), this ),
  m_nonlinearSolverParameters( groupKeyStruct::nonlinearSolverParametersString(), this ),
  m_solverStatistics( groupKeyStruct::solverStatisticsString(), this ),
  m_systemSetupTimestamp( 0 ),
  m_activeCoupledIterations( 1 ),
  m_enableIsoThermalEstimator( 0 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  this->registerWrapper( viewKeyStruct::activeCoupledIterationsString(), &m_activeCoupledIterations ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Number of coupled iterations to activate estimator solve." );

  registerWrapper( viewKeyStruct::enableIsoThermalEstimatorString(), &m_enableIsoThermalEstimator ).
    setDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Estimator configuration option to disable thermal effects on initial well constraint solve and then converge solution with thermal effects enabled: \n"
                    " - If the flag is set to 1, thermal effects are enabled during the initial constraint solve. \n"
                    " - If the flag is set to 0, thermal effects are disabled during the initial constraint solve." );



  registerWrapper( viewKeyStruct::writeLinearSystemString(), &m_writeLinearSystem ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Write matrix, rhs, solution to screen ( = 1) or file ( = 2)." );

  registerWrapper( viewKeyStruct::allowNonConvergedLinearSolverSolutionString(), &m_allowNonConvergedLinearSolverSolution ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Cut time step if linear solution fail without going until max nonlinear iterations." );

  registerWrapper( viewKeyStruct::usePhysicsScalingString(), &m_usePhysicsScaling ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Enable physics-based scaling of the linear system. Default: true." );


  registerWrapper( viewKeyStruct::writeStatisticsCSVString(), &m_writeStatisticsCSV ).
    setApplyDefaultValue( StatsOutputType::none ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( GEOS_FMT( "When set to `{}`, output iterations information to a csv\n"
                              "When set to `{}`, output convergence information to a csv\n"
                              "When set to `{}` output both convergence & iteration information to a csv.",
                              EnumStrings< StatsOutputType >::toString( StatsOutputType::iteration ),
                              EnumStrings< StatsOutputType >::toString( StatsOutputType::convergence ),
                              EnumStrings< StatsOutputType >::toString( StatsOutputType::all ) ));

  addLogLevel< logInfo::Convergence >();
  addLogLevel< logInfo::Fields >();
  addLogLevel< logInfo::LinearSolver >();
  addLogLevel< logInfo::ResidualNorm >();
  addLogLevel< logInfo::Solution >();
  addLogLevel< logInfo::TimeStep >();
  addLogLevel< logInfo::Timers >();

  registerGroup( groupKeyStruct::linearSolverParametersString(), &m_linearSolverParameters );
  registerGroup( groupKeyStruct::nonlinearSolverParametersString(), &m_nonlinearSolverParameters );
  registerGroup( groupKeyStruct::solverStatisticsString(), &m_solverStatistics );

  m_localMatrix.setName( this->getName() + "/localMatrix" );
  m_matrix.setDofManager( &m_dofManager );
}

void WellNewtonSolver::postInputInitialization()
{
  m_solverStatistics.setOutputFilesName( getName() );

  m_solverStatistics.makeDir( m_writeStatisticsCSV !=  StatsOutputType::none );

  getIterationStats().setTableName( getName() );
  getIterationStats().setLogOutputRequest( true );
  getIterationStats().setCSVOutputRequest( m_writeStatisticsCSV == StatsOutputType::iteration ||
                                           m_writeStatisticsCSV == StatsOutputType::all );
  getConvergenceStats().setCSVOutputRequest( m_writeStatisticsCSV == StatsOutputType::convergence ||
                                             m_writeStatisticsCSV == StatsOutputType::all );
}

WellNewtonSolver::~WellNewtonSolver() = default;

#if 0
void WellNewtonSolver::initialize_postMeshGeneration()
{
  //ExecutableGroup::initialize_postMeshGeneration();
  //DomainPartition const & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  //generateMeshTargetsFromTargetRegions( domain.getMeshBodies());
}
#endif
void WellNewtonSolver::generateMeshTargetsFromTargetRegions( Group const & meshBodies )
{
  for( auto const & target : m_targetRegionNames )
  {

    stdVector< string > targetTokens = stringutilities::tokenize( target, "/" );

    if( targetTokens.size()==1 ) // no MeshBody or MeshLevel specified
    {
      GEOS_ERROR_IF( meshBodies.numSubGroups() != 1,
                     getDataContext() << ": No MeshBody information is specified in" <<
                     " WellNewtonSolver::meshTargets, but there are multiple MeshBody objects",
                     getDataContext() );
      MeshBody const & meshBody = meshBodies.getGroup< MeshBody >( 0 );
      string const meshBodyName = meshBody.getName();

      string const meshLevelName = ""; //tjbm_discretizationName;

      string const regionName = target;
      auto const key = std::make_pair( meshBodyName, meshLevelName );
      m_meshTargets[key].emplace_back( regionName );
    }
    else if( targetTokens.size()==2 )
    {
      string const meshBodyName = targetTokens[0];
      GEOS_ERROR_IF( !meshBodies.hasGroup( meshBodyName ),
                     getWrapperDataContext( viewKeyStruct::targetRegionsString() ) << ": MeshBody (" <<
                     meshBodyName << ") is specified in targetRegions, but does not exist.",
                     getWrapperDataContext( viewKeyStruct::targetRegionsString() ) );

      string const meshLevelName = "";//tjbm_discretizationName;

      string const regionName = targetTokens[1];


      auto const key = std::make_pair( meshBodyName, meshLevelName );
      m_meshTargets[key].emplace_back( regionName );
    }
    else
    {
      GEOS_ERROR( getDataContext() << ": Invalid specification of targetRegions" );
    }
  }
}



Group * WellNewtonSolver::createChild( string const & GEOS_UNUSED_PARAM( childKey ), string const & GEOS_UNUSED_PARAM( childName ) )
{
  // Unused as all children are created within the constructor
  return nullptr;
}

WellNewtonSolver::CatalogInterface::CatalogType & WellNewtonSolver::getCatalog()
{
  static WellNewtonSolver::CatalogInterface::CatalogType catalog;
  return catalog;
}


bool WellNewtonSolver::registerCallback( void * func, const std::type_info & funcType )
{
  if( std::type_index( funcType ) == std::type_index( typeid( std::function< void( CRSMatrix< real64, globalIndex >, array1d< real64 > ) > ) ) )
  {
    m_assemblyCallback = *reinterpret_cast< std::function< void( CRSMatrix< real64, globalIndex >, array1d< real64 > ) > * >( func );
    return true;
  }

  return false;
}



void WellNewtonSolver::logEndOfCycleInformation( integer const cycleNumber,
                                                 integer const numOfSubSteps,
                                                 stdVector< real64 > const & subStepDts ) const
{
  LogPart logpart( "TIMESTEP", MpiWrapper::commRank() == 0 );
  logpart.addEndDescription( "- Cycle ", cycleNumber );
  logpart.addEndDescription( "- N substeps ", numOfSubSteps );

  std::stringstream logMessage;
  for( integer i = 0; i < numOfSubSteps; ++i )
  {
    if( i > 0 )
    {
      logMessage << ", ";
    }
    logMessage << subStepDts[i] << " " << units::getSymbol( units::Unit::Time );
  }

  if( logMessage.rdbuf()->in_avail() == 0 )
    logMessage << "/";

  logpart.addEndDescription( "- substep dts ", logMessage.str() );
  logpart.end();

  if( isLogLevelActive< logInfo::SolverExecutionDetails >( getLogLevel()))
    getIterationStats().outputStatistics();
}

#if 0
// tjb keep this history
void WellNewtonSolver::implicitStepSetup( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                          real64 const & GEOS_UNUSED_PARAM( dt ),
                                          DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  // clean the solution history
  while( m_solutionHistory.size() > 0 )
  {
    m_solutionHistory.eraseArray( 0 );
  }
}
#endif
void WellNewtonSolver::setupDofs( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                  DofManager & GEOS_UNUSED_PARAM( dofManager ) ) const
{
  GEOS_ERROR( "WellNewtonSolver::setupDofs called!. Should be overridden." );
}

void WellNewtonSolver::setSparsityPattern( DomainPartition & GEOS_UNUSED_PARAM( domain ),
                                           DofManager & dofManager,
                                           CRSMatrix< real64, globalIndex > & GEOS_UNUSED_PARAM( localMatrix ),
                                           SparsityPattern< globalIndex > & pattern )
{
  dofManager.setSparsityPattern( pattern );
}

void WellNewtonSolver::setSystemSetupTimestamp( Timestamp timestamp )
{
  m_systemSetupTimestamp = timestamp;

  std::ostringstream oss;
  m_dofManager.printFieldInfo( oss );
  GEOS_LOG_LEVEL( logInfo::Fields, oss.str());
}

std::unique_ptr< PreconditionerBase< LAInterface > >
WellNewtonSolver::createPreconditioner( DomainPartition & GEOS_UNUSED_PARAM( domain ) ) const
{
  // By default, do not create a preconditioner, one will be created internally inside LA backend
  return {};

  // TODO: refactor interfaces to always create preconditioner externally and pass to backends
  // return LAInterface::createPreconditioner( m_linearSolverParameters.get() );
}


namespace
{

/**
 * @brief Helper for debug output of linear algebra objects (matrices and vectors)
 * @tparam T type of LA object (must have stream insertion and .write() implemented)
 * @param obj                the object to output
 * @param cycleNumber        event cycle number
 * @param nonlinearIteration nonlinear iteration number
 * @param filePrefix          short filename prefix (e.g. "mat")
 * @param screenName           long name for screen output (e.g. "System matrix")
 * @param toScreen           whether to print on screen
 * @param toFile             whether to write to file
 */
template< typename T >
void debugOutputLAObject( T const & obj,
                          real64 const & GEOS_UNUSED_PARAM( time ),
                          integer const cycleNumber,
                          integer const nonlinearIteration,
                          string const & filePrefix,
                          string const & screenName,
                          bool const toScreen,
                          bool const toFile )
{
  if( toScreen )
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "{2:=>{1}}\n{0}:\n{2:=>{1}}", screenName, screenName.size() + 1, "" ) );
    GEOS_LOG( obj );
  }

  if( toFile )
  {
    string const filename = GEOS_FMT( "{}_{:06}_{:02}.mtx", filePrefix.c_str(), cycleNumber, nonlinearIteration );
    obj.write( filename, LAIOutputFormat::MATRIX_MARKET );
    GEOS_LOG_RANK_0( GEOS_FMT( "{} written to {}", screenName, filename ) );
  }
}

}

void WellNewtonSolver::debugOutputSystem( real64 const & time,
                                          integer const cycleNumber,
                                          integer const nonlinearIteration,
                                          ParallelMatrix const & matrix,
                                          ParallelVector const & rhs ) const
{
  // special case when flag value > 2
  if( m_writeLinearSystem > 2 && cycleNumber < m_writeLinearSystem )
    return;

  debugOutputLAObject( matrix,
                       time,
                       cycleNumber,
                       nonlinearIteration,
                       getName() + "_mat",
                       "System matrix",
                       m_writeLinearSystem == 1,
                       m_writeLinearSystem >= 2 );

  debugOutputLAObject( rhs,
                       time,
                       cycleNumber,
                       nonlinearIteration,
                       getName() + "_rhs",
                       "System right-hand side",
                       m_writeLinearSystem == 1,
                       m_writeLinearSystem >= 2 );
}

void WellNewtonSolver::debugOutputSolution( real64 const & time,
                                            integer const cycleNumber,
                                            integer const nonlinearIteration,
                                            ParallelVector const & solution ) const
{
  // special case when flag value > 2
  if( m_writeLinearSystem > 2 && cycleNumber < m_writeLinearSystem )
    return;

  debugOutputLAObject( solution,
                       time,
                       cycleNumber,
                       nonlinearIteration,
                       getName() + "_sol",
                       "System solution",
                       m_writeLinearSystem == 1,
                       m_writeLinearSystem >= 2 );
}

void WellNewtonSolver::updateAndWriteConvergenceStep( real64 const & time_n, real64 const & dt,
                                                      integer const cycleNumber, integer const iteration )
{
  getConvergenceStats().updateSolverStep( time_n, dt, cycleNumber, iteration );
  getConvergenceStats().writeConvergenceStatsToTable();
}


void WellNewtonSolver::solveLinearSystem( DofManager const & dofManager,
                                          ParallelMatrix & matrix,
                                          ParallelVector & rhs,
                                          ParallelVector & solution )
{
  GEOS_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  LinearSolverParameters const & params = m_linearSolverParameters.get();
  const bool isDirectSolver = (params.solverType == LinearSolverParameters::SolverType::direct);
  const bool isSetupNeeded = !(isDirectSolver && params.direct.reuseFactorization);

  matrix.setDofManager( &dofManager );

  GEOS_WARNING_IF( isDirectSolver && dofManager.numGlobalDofs() > 100000,
                   "Direct solver used for large system ( > 100,000 DOFs ). "
                   "This may lead to high memory consumption and long computation times. "
                   "Consider using an iterative solver for better performance." );

  // Apply physics-based scaling to the linear system if enabled
  if( m_usePhysicsScaling )
  {
    Timer timer_setup( m_timers.get_inserted( "linear solver scaling" ) );

    matrix.computeScalingVector( m_scaling );
    matrix.leftRightScale( m_scaling, m_scaling );
    rhs.pointwiseProduct( m_scaling );
    // Assume the solution is zeroed out, thus no need to scale it
  }

  if( isDirectSolver || !m_precond )
  {
    if( !m_linearSolver )
    {
      m_linearSolver = LAInterface::createSolver( params );
    }

    if( isSetupNeeded )
    {
      Timer timer_setup( m_timers.get_inserted( "linear solver setup" ) );
      m_linearSolver->setup( matrix );
    }

    {
      Timer timer_setup( m_timers.get_inserted( "linear solver solve" ) );
      m_linearSolver->solve( rhs, solution );
    }

    m_linearSolverResult = m_linearSolver->result();
  }
  else
  {
    {
      Timer timer_setup( m_timers.get_inserted( "linear solver setup" ) );
      m_precond->setup( matrix );
    }
    std::unique_ptr< KrylovSolver< ParallelVector > > solver = KrylovSolver< ParallelVector >::create( params, matrix, *m_precond );
    {
      Timer timer_setup( m_timers.get_inserted( "linear solver solve" ) );
      solver->solve( rhs, solution );
    }
    m_linearSolverResult = solver->result();
  }

  getIterationStats().accumulateSolverLinearTime( m_linearSolverResult.setupTime, m_linearSolverResult.solveTime );

  GEOS_LOG_LEVEL_RANK_0( logInfo::LinearSolver,
                         GEOS_FMT( "        Linear solve: ( iter, res ) = ( {:3}, {:4.2e} )",
                                   m_linearSolverResult.numIterations,
                                   m_linearSolverResult.residualReduction ));

  if( params.stopIfError )
  {
    GEOS_ERROR_IF( m_linearSolverResult.breakdown(),
                   getDataContext() << ": Linear solution breakdown -> simulation STOP",
                   getDataContext() );
  }
  else
  {
    GEOS_WARNING_IF( !m_linearSolverResult.success(),
                     getDataContext() << ": Linear solution failed",
                     getDataContext() );
  }

  // Unscale the solution vector if physics-based scaling was applied
  if( m_usePhysicsScaling )
  {
    Timer timer_setup( m_timers.get_inserted( "linear solver scaling" ) );

    solution.pointwiseProduct( m_scaling );
  }
}



// Detect oscillations for all dofs in the solution history
bool WellNewtonSolver::detectOscillations() const
{
  // grab the parameters
  integer const oscillationCheckDepth = m_nonlinearSolverParameters.m_oscillationCheckDepth;
  real64 const oscillationTolerance = m_nonlinearSolverParameters.m_oscillationTolerance;
  real64 const oscillationFraction = m_nonlinearSolverParameters.m_oscillationFraction;

  if( m_solutionHistory.size() < oscillationCheckDepth )
    return false; // not enough history to check oscillations

  RAJA::ReduceSum< parallelDeviceReduce, localIndex > oscillationCount( 0 );

  auto const solutionHistory = m_solutionHistory.toViewConst();
  localIndex const numDofs = m_solutionHistory[0].size();
  localIndex const historySize = m_solutionHistory.size();

  RAJA::forall< parallelDevicePolicy<> >( RAJA::TypedRangeSegment< localIndex >( 0, numDofs ),
                                          [=] GEOS_HOST_DEVICE ( localIndex const dof )
  {
    bool oscillationDetected = true;
    for( localIndex i = historySize - 1; i > historySize - oscillationCheckDepth; --i )
    {
      real64 dxCur = solutionHistory[i][dof];
      real64 dxPrev = solutionHistory[i-1][dof];

      if( LvArray::math::abs( dxCur ) < oscillationTolerance || LvArray::math::abs( dxPrev ) < oscillationTolerance )
      {
        oscillationDetected = false;
        break;   // solution changes are too small
      }

      real64 maxAbs = LvArray::math::max( LvArray::math::abs( dxCur ), LvArray::math::abs( dxPrev ) );
      if( LvArray::math::abs( dxCur + dxPrev ) / maxAbs > oscillationTolerance )
      {
        oscillationDetected = false;
        break;   // solution changes are not oscillating
      }

      if( dxCur * dxPrev > 0 )
      {
        oscillationDetected = false;
        break;   // sign is not oscillating
      }
    }

    if( oscillationDetected )
    {
      oscillationCount += 1;
    }
  } );

  real64 const f = static_cast< real64 >( MpiWrapper::sum( oscillationCount.get() ) ) / MpiWrapper::sum( numDofs );

  return f > oscillationFraction;
}


} // namespace geos
