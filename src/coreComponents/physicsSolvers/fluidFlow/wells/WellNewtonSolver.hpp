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
 * @file WellNewtonSolver.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WELLNEWTONSOLVER_HPP_
#define GEOS_PHYSICSSOLVERS_WELLNEWTONSOLVER_HPP_

#include "codingUtilities/traits.hpp"
#include "common/DataTypes.hpp"
#include "common/format/LogPart.hpp"

#include "dataRepository/RestartFlags.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/utilities/LinearSolverResult.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "physicsSolvers/NonlinearSolverParameters.hpp"
#include "physicsSolvers/LinearSolverParameters.hpp"
#include "physicsSolvers/SolverStatistics.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"

#include "physicsSolvers/fluidFlow/SolutionCheckHelpers.hpp"
#include "common/Timer.hpp"
#include <limits>

namespace geos
{

class DomainPartition;

/**
 * @class WellNewtonSolver
 * @brief Base class for all physics solvers
 *
 * This class provides the base interface for all physics solvers. It provides the basic
 * functionality for setting up and solving a linear system, as well as the interface for
 * performing a timestep.
 */
class WellNewtonSolver : public dataRepository::Group
{
public:

  /**
   * @brief Type of the stat output
   */
  enum class StatsOutputType : integer
  {
    none, iteration, convergence, all
  };

  /**
   * @brief Constructor for WellNewtonSolver
   * @param name the name of this instantiation of WellNewtonSolver
   * @param parent the parent group of this instantiation of WellNewtonSolver
   */
  explicit WellNewtonSolver( string const & name,
                             Group * const parent );

  /**
   * @brief Move constructor for WellNewtonSolver
   */
  WellNewtonSolver( WellNewtonSolver && ) = default;

  /**
   * @brief Destructor for WellNewtonSolver
   */
  virtual ~WellNewtonSolver() override;

  /**
   * @brief Deleted constructor
   */
  WellNewtonSolver() = delete;

  /**
   * @brief Deleted copy constructor
   */
  WellNewtonSolver( WellNewtonSolver const & ) = delete;

  /**
   * @brief Deleted copy assignment operator
   */
  WellNewtonSolver & operator=( WellNewtonSolver const & ) = delete;

  /**
   * @brief Deleted move assignment operator
   */
  WellNewtonSolver & operator=( WellNewtonSolver && ) = delete;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName() { return "WellNewtonSolver"; }



  /**
   * @brief Generate mesh targets from target regions
   * @param meshBodies the group of mesh bodies
   */
  void generateMeshTargetsFromTargetRegions( Group const & meshBodies );


  /**
   *
   * @brief Getter for system matrix
   * @return a reference to linear system matrix of this solver
   */
  ParallelMatrix & getSystemMatrix() { return m_matrix; }

  /**
   * @brief Getter for system rhs vector
   * @return a reference to linear system right-hand side of this solver
   */
  ParallelMatrix const & getSystemMatrix() const { return m_matrix; }

  /**
   * @brief Getter for system rhs vector
   * @return a reference to linear system right-hand side of this solver
   */
  ParallelVector & getSystemRhs() { return m_rhs; }

  /**
   * @brief Getter for system rhs vector
   * @return a reference to linear system right-hand side of this solver
   */
  ParallelVector const & getSystemRhs() const { return m_rhs; }

  /**
   * @brief Getter for system solution vector
   * @return a reference to solution vector of this solver
   */
  ParallelVector & getSystemSolution() { return m_solution; }

  /**
   * @brief Getter for system solution vector
   * @return a reference to solution vector of this solver
   */
  ParallelVector const & getSystemSolution() const { return m_solution; }

  /**
   * @brief Getter for degree-of-freedom manager
   * @return a reference to degree-of-freedom manager of this solver
   */
  DofManager & getDofManager() { return m_dofManager; }

  /**
   * @brief Getter for degree-of-freedom manager
   * @return a reference to degree-of-freedom manager of this solver
   */
  DofManager const & getDofManager() const { return m_dofManager; }

  /**
   * @brief Getter for local matrix
   * @return a reference to linear system matrix of this solver
   */
  CRSMatrix< real64, globalIndex > & getLocalMatrix() { return m_localMatrix; }

  /**
   * @brief Getter for local matrix
   * @return a reference to linear system matrix of this solver
   */
  CRSMatrixView< real64 const, globalIndex const > getLocalMatrix() const { return m_localMatrix.toViewConst(); }


  template< typename T >
  void  setupSystem( T & well, DomainPartition & domain,
                     std::string const & meshBodyName,
                     MeshLevel const & meshLevel,
                     WellElementRegion & wellElementRegion,
                     bool const setSparsity =true );

  template< typename T >
  bool
  solveNonlinearSystem( T & well, real64 const & time_n,
                        real64 const & stepDt,
                        integer const cycleNumber,
                        DomainPartition & domain,
                        MeshLevel & mesh,
                        ElementRegionManager & elemManager,
                        WellElementSubRegion & subRegion );


  /**
   * @brief Populate degree-of-freedom manager with fields relevant to this solver
   * @param domain the domain containing the mesh and fields
   * @param dofManager degree-of-freedom manager associated with the linear system
   */
  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const;

  /**
   * @brief Set up the linear system (DOF indices and sparsity patterns)
   * @param domain the domain containing the mesh and fields
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   * @param setSparsity flag to indicate if the sparsity pattern should be set
   *
   * @note While the function is virtual, the base class implementation should be
   *       sufficient for most single-physics solvers.
   */


  /**
   * @brief Set the sparsity pattern of the linear system matrix
   * @param domain the domain containing the mesh and fields
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param pattern the sparsity pattern to be filled
   */
  virtual void
  setSparsityPattern( DomainPartition & domain,
                      DofManager & dofManager,
                      CRSMatrix< real64, globalIndex > & localMatrix,
                      SparsityPattern< globalIndex > & pattern );

  /**
   * @brief Create a preconditioner for this solver's linear system.
   * @param domain the domain containing the mesh and fields
   * @return the newly created preconditioner object
   */
  virtual std::unique_ptr< PreconditionerBase< LAInterface > >
  createPreconditioner( DomainPartition & domain ) const;


  /**
   * @brief Output the assembled linear system for debug purposes.
   * @param time beginning-of-step time
   * @param cycleNumber event cycle number
   * @param nonlinearIteration current nonlinear iteration number
   * @param matrix system matrix
   * @param rhs system right-hand side vector
   */
  void
  debugOutputSystem( real64 const & time,
                     integer const cycleNumber,
                     integer const nonlinearIteration,
                     ParallelMatrix const & matrix,
                     ParallelVector const & rhs ) const;

  /**
   * @brief Output the linear system solution for debug purposes.
   * @param time beginning-of-step time
   * @param cycleNumber event cycle number
   * @param nonlinearIteration current nonlinear iteration number
   * @param solution system solution vector
   */
  void
  debugOutputSolution( real64 const & time,
                       integer const cycleNumber,
                       integer const nonlinearIteration,
                       ParallelVector const & solution ) const;

  /**
   * @brief Update the convergence information and write then into a CSV file
   * @param time_n the time at the beginning of the step
   * @param dt the desired timestep
   * @param cycleNumber event cycle number
   * @param iteration current iteration
   */
  virtual void
  updateAndWriteConvergenceStep( real64 const & time_n,
                                 real64 const & dt,
                                 integer const cycleNumber,
                                 integer const iteration );



  /**
   * @brief function to apply a linear system solver to the assembled system.
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   *
   * This function calls the linear solver package to perform a single linear solve on the block
   * system. The derived physics solver is required to specify the call, as no default is provided.
   *
   * @note This function must be overridden in the derived physics solver in order to use an implict
   * solution method such as LinearImplicitStep() or NonlinearImplicitStep().
   */
  virtual void
  solveLinearSystem( DofManager const & dofManager,
                     ParallelMatrix & matrix,
                     ParallelVector & rhs,
                     ParallelVector & solution );



  /**
   * @brief creates a child group of of this WellNewtonSolver instantiation
   * @param childKey the key of the child type
   * @param childName the name of the child
   * @return a pointer to the child group
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /**
   * @brief Type alias for catalog interface used by this class. See CatalogInterface.
   */
  using CatalogInterface = dataRepository::CatalogInterface< WellNewtonSolver, string const &, Group * const >;

  /**
   * @brief Get the singleton catalog for WellNewtonSolver.
   * @return reference to the catalog object
   */
  static CatalogInterface::CatalogType & getCatalog();

  /**
   * @brief Structure to hold scoped key names
   */
  struct viewKeyStruct
  {
    /// @return string for the cflFactor wrapper
    static constexpr char const * cflFactorString() { return "cflFactor"; }

    /// @return string for the initialDt wrapper
    static constexpr char const * initialDtString() { return "initialDt"; }

    /// @return string for the minDtIncreaseInterval wrapper
    static constexpr char const * minDtIncreaseIntervalString() { return "minDtIncreaseInterval"; }

    /// @return string for the discretization wrapper
    static constexpr char const * discretizationString() { return "discretization"; }

    /// @return string for the nextDt targetRegions wrapper
    static constexpr char const * targetRegionsString() { return "targetRegions"; }

    /// @return string for the writeLinearSystem wrapper
    static constexpr char const * writeLinearSystemString() { return "writeLinearSystem"; }

    /// @return string for the usePhysicsScaling wrapper
    static constexpr char const * usePhysicsScalingString() { return "usePhysicsScaling"; }

    /// @return string for the allowNonConvergedLinearSolverSolution wrapper
    static constexpr char const * allowNonConvergedLinearSolverSolutionString() { return "allowNonConvergedLinearSolverSolution"; }

    /// @return string for the writeStatistics wrapper
    static constexpr char const * writeStatisticsCSVString() { return "writeStatistics"; }

    /// @return string for the numTimestepsSinceLastDtCut wrapper
    static constexpr char const * numTimestepsSinceLastDtCutString() { return "numTimestepsSinceLastDtCut"; }

    /// string key for the esitmate well solution flag
    static constexpr char const * activeCoupledIterationsString() { return "activeCoupledIterations"; }
    static constexpr char const * estimateWellSolutionString() { return "estimateWellSolution"; }
    /// string key for the enable iso thermal estimator flag
    static constexpr char const * enableIsoThermalEstimatorString() { return "enableIsoThermalEstimator"; }
  };

  /**
   * @brief Structure to hold scoped key names
   */
  struct groupKeyStruct
  {
    /// @return string for the linearSolverParameters wrapper
    static constexpr char const * linearSolverParametersString() { return "LinearSolverParameters"; }

    /// @return string for the nonlinearSolverParameters wrapper
    static constexpr char const * nonlinearSolverParametersString() { return "NonlinearSolverParameters"; }

    /// @return string for the solverStatistics wrapper
    static constexpr char const * solverStatisticsString() { return "SolverStatistics"; }
  };

  /**
   * @brief getter for the timestamp of the system setup
   * @return the timestamp of the last time systemSetup was called
   */
  Timestamp getSystemSetupTimestamp() const { return m_systemSetupTimestamp; }


  /**
   * @brief set the timestamp of the system setup
   * @param[in] timestamp the new timestamp of system setup
   */
  void setSystemSetupTimestamp( Timestamp timestamp );


  /**
   * @brief accessor for the linear solver parameters.
   * @return the linear solver parameter list
   */
  LinearSolverParameters & getLinearSolverParameters()
  {
    return m_linearSolverParameters.get();
  }

  /**
   * @brief const accessor for the linear solver parameters.
   * @return the linear solver parameter list
   */
  LinearSolverParameters const & getLinearSolverParameters() const
  {
    return m_linearSolverParameters.get();
  }

  /**
   * @brief accessor for the nonlinear solver parameters.
   * @return the nonlinear solver parameter list
   */
  NonlinearSolverParameters & getNonlinearSolverParameters()
  {
    return m_nonlinearSolverParameters;
  }

  /**
   * @brief const accessor for the nonlinear solver parameters.
   * @return the nonlinear solver parameter list
   */
  NonlinearSolverParameters const & getNonlinearSolverParameters() const
  {
    return m_nonlinearSolverParameters;
  }

  /**
   * @brief synchronize the nonlinear solver parameters.
   */
  virtual void
  synchronizeNonlinearSolverParameters()
  { /* empty here, overriden in CoupledSolver */ }

  /**
   * @brief Get position of a given region within solver's target region list
   * @param regionName the region name to find
   * @return index within target regions list
   */
  localIndex targetRegionIndex( string const & regionName ) const;

  /**
   * @brief return the list of target regions
   * @return the array of region names
   */
  string_array const & getTargetRegionNames() const {return m_targetRegionNames;}



  /**
   * @brief function to set the value of m_assemblyCallback
   * @param func the function to set m_assemblyCallback to
   * @param funcType the type of the function
   * @return true if the function was successfully set, false otherwise
   *
   * This is used to provide a callback function for to be called in the assembly step.
   */
  virtual bool registerCallback( void * func, const std::type_info & funcType ) final override;

  /**
   * @return An IterationsStatistics for the "root" solver.
   * Otherwise return an empty IterationsStatistics
   */
  IterationsStatistics & getIterationStats()
  {
    return m_solverStatistics.m_iterationsStats;
  }
  /**
   * @return An IterationsStatistics for the "root" solver.
   * Otherwise return an empty IterationsStatistics
   * (const version)
   */
  IterationsStatistics const & getIterationStats() const
  {
    return m_solverStatistics.m_iterationsStats;
  }
  /**
   * @return A ConvergenceStatistics for all sub-solvers
   */
  ConvergenceStatistics & getConvergenceStats()
  {
    return m_solverStatistics.m_convergenceStats;
  }
  /**
   * @return A ConvergenceStatistics for all sub-solvers (const version)
   */
  ConvergenceStatistics const & getConvergenceStats() const
  {
    return m_solverStatistics.m_convergenceStats;
  }

  /**
   * @brief accessor for the solver statistics.
   * @return reference to m_solverStatistics
   */
  SolverStatistics & getSolverStatistics() { return m_solverStatistics; }

  /**
   * @brief const accessor for the solver statistics.
   * @return reference to m_solverStatistics
   */
  SolverStatistics const & getSolverStatistics() const { return m_solverStatistics; }



  /**
   * @brief Detect oscillations in the solution
   * @return true if oscillations are detected, false otherwise
   */
  bool detectOscillations() const;


  /**
   * @brief Set thermal effects enable
   * @param[in] true/false
   */
  void enableThermalEffects ( bool enable ) { m_thermalEffectsEnabled = enable; };

  /**
   * @brief Are thermal effects enabled
   * @return true if thermal effects are enabled, false otherwise
   */
  bool thermalEffectsEnabled() const { return m_thermalEffectsEnabled; }

  /**
   * @brief Is isoThermalEstimator  enabled
   * @return true if isoThermalEstimator is enabled, false otherwise
   */
  bool isoThermalEstimatorEnabled() const { return m_enableIsoThermalEstimator; }

  bool getNumActiveCoupledIterations() const { return m_activeCoupledIterations; }

protected:

  virtual void postInputInitialization() override;

  /// behavior in case of linear solver failure
  integer m_allowNonConvergedLinearSolverSolution;



  /// Data structure to handle degrees of freedom
  DofManager m_dofManager;

  /// System matrix
  ParallelMatrix m_matrix;

  /// System right-hand side vector
  ParallelVector m_rhs;

  /// System solution vector
  ParallelVector m_solution;

  /// Diagonal scaling vector D (Ahat = D * A * D, bhat = D * b, x = D * xhat)
  ParallelVector m_scaling;

  /// Flag to decide whether to apply physics-based scaling to the linear system
  integer m_usePhysicsScaling;

  /// Local system matrix and rhs
  CRSMatrix< real64, globalIndex > m_localMatrix;

  /// Custom linear solver for the "native" solver type
  std::unique_ptr< LinearSolverBase< LAInterface > > m_linearSolver;

  /// Custom preconditioner for the "native" iterative solver
  std::unique_ptr< PreconditionerBase< LAInterface > > m_precond;

  /// flag for debug output of matrix, rhs, and solution
  integer m_writeLinearSystem;

  /// Parameter for outputing statistics information
  StatsOutputType m_writeStatisticsCSV;

  /// Linear solver parameters
  LinearSolverParametersInput m_linearSolverParameters;

  /// Result of the last linear solver
  LinearSolverResult m_linearSolverResult;

  /// Nonlinear solver parameters
  NonlinearSolverParameters m_nonlinearSolverParameters;

  /// Solver statistics
  SolverStatistics m_solverStatistics;

  /// Timestamp of the last call to setup system
  Timestamp m_systemSetupTimestamp;

  /// Callback function for assembly step
  std::function< void( CRSMatrix< real64, globalIndex >, array1d< real64 > ) > m_assemblyCallback;

  /// Timers for the aggregate profiling of the solver
  stdMap< std::string, std::chrono::system_clock::duration > m_timers;

  /// History of the solution vector, used for oscillation detection
  ArrayOfArrays< real64 > m_solutionHistory;

private:
  /// List of names of regions the solver will be applied to
  string_array m_targetRegionNames;

  /// Map containing the array of target regions (value) for each MeshBody (key).
  map< std::pair< string, string >, string_array > m_meshTargets;

  /// Number of coupled iterations to activate estimator solve
  integer m_activeCoupledIterations;

  /// Flag to enable thermal effects in wellbore calculations
  bool m_thermalEffectsEnabled;
  integer m_enableIsoThermalEstimator;


  /**
   * @brief output information about the cycle to the log
   * @param cycleNumber the current cycle number
   * @param numOfSubSteps the number of substeps taken
   * @param subStepDts the time step size for each substep
   */
  void logEndOfCycleInformation( integer const cycleNumber,
                                 integer const numOfSubSteps,
                                 stdVector< real64 > const & subStepDts ) const;
};

template< typename T >
void WellNewtonSolver::setupSystem( T & well, DomainPartition & domain,
                                    std::string const & meshBodyName,
                                    MeshLevel const & meshLevel,
                                    WellElementRegion & wellElementRegion,
                                    bool const setSparsity )
{
  GEOS_MARK_FUNCTION;

  map< std::pair< string, string >, string_array > meshTargets;
  string_array regions;

  meshTargets.clear();
  regions.clear();
  regions.emplace_back( wellElementRegion.getName() );
  auto const key = std::make_pair( meshBodyName, meshLevel.getName() );
  meshTargets[key] = std::move( regions );

  m_dofManager.setDomain( domain );
  m_dofManager.addField( well.wellElementDofName(),
                         FieldLocation::Elem,
                         well.numDofPerWellElement(),
                         meshTargets );

  m_dofManager.addCoupling( well.wellElementDofName(),
                            well.wellElementDofName(),
                            DofManager::Connector::Node );

  m_dofManager.reorderByRank();
  if( setSparsity )
  {
    SparsityPattern< globalIndex > pattern;
    setSparsityPattern( domain, m_dofManager, m_localMatrix, pattern );
    m_localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
  }
  m_localMatrix.setName( this->getName() + "/matrix" );

  m_rhs.setName( this->getName() + "/rhs" );
  m_rhs.create( m_dofManager.numLocalDofs(), MPI_COMM_GEOS );

  m_solution.setName( this->getName() + "/solution" );
  m_solution.create( m_dofManager.numLocalDofs(), MPI_COMM_GEOS );
}


template< typename T >
bool WellNewtonSolver::solveNonlinearSystem( T & well, real64 const & time_n,
                                             real64 const & stepDt,
                                             integer const cycleNumber,
                                             DomainPartition & domain,
                                             MeshLevel & mesh,
                                             ElementRegionManager & elemManager,
                                             WellElementSubRegion & subRegion )
{
  integer const maxNewtonIter = m_nonlinearSolverParameters.m_maxIterNewton;
  integer dtAttempt = m_nonlinearSolverParameters.m_numTimeStepAttempts;
  integer configurationLoopIter = m_nonlinearSolverParameters.m_numConfigurationAttempts;
  integer const minNewtonIter = m_nonlinearSolverParameters.m_minIterNewton;
  real64 const newtonTol = m_nonlinearSolverParameters.m_newtonTol;

// keep residual from previous iteration in case we need to do a line search

  integer newtonIter = 0;
  real64 scaleFactor = 1.0;

  bool isNewtonConverged = false;

  for( newtonIter = 0; newtonIter < maxNewtonIter; ++newtonIter )
  {
    if( m_nonlinearSolverParameters.getLogLevel() > 4 )
      GEOS_LOG_LEVEL_RANK_0( logInfo::NonlinearSolver,
                             GEOS_FMT( " Well: {}   Est Attempt: NewtonIter: {:2}", subRegion.getName(), stepDt, newtonIter ));

    {
      Timer timer( m_timers.get_inserted( "assemble" ) );

// We sync the nonlinear convergence history. The coupled solver parameters are the one being
// used. We want to propagate the info to subsolvers. It can be important for solvers that
// have special treatment for specific iterations.
      synchronizeNonlinearSolverParameters();

// zero out matrix/rhs before assembly
      m_localMatrix.zero();
      m_rhs.zero();

      arrayView1d< real64 > const localRhs = m_rhs.open();

// call assemble to fill the matrix and the rhs
      well.assembleSystem( time_n,
                           stepDt,
                           cycleNumber,
                           elemManager,
                           subRegion,
                           m_dofManager,
                           m_localMatrix.toViewConstSizes(),
                           localRhs );

// apply boundary conditions to system
      well.applyWellBoundaryConditions( time_n,
                                        stepDt,
                                        elemManager,
                                        subRegion,
                                        m_dofManager,
                                        localRhs,
                                        m_localMatrix.toViewConstSizes() );

      m_rhs.close();

      if( m_assemblyCallback )
      {
// Make a copy of LA objects and ship off to the callback
        array1d< real64 > localRhsCopy( m_rhs.localSize() );
        localRhsCopy.setValues< parallelDevicePolicy<> >( m_rhs.values() );
        m_assemblyCallback( m_localMatrix, std::move( localRhsCopy ) );
      }
    }

    // well.outputSingleWellDebug( time_n, stepDt, 0, newtonIter, 0,
    //                             mesh, subRegion, dofManager, m_localMatrix.toViewConstSizes(), m_rhs.values()  );
    real64 residualNorm = 0;
    {
      Timer timer( m_timers.get_inserted( "convergence check" ) );

// get residual norm
      residualNorm = well.calculateWellResidualNorm( time_n, stepDt, m_nonlinearSolverParameters, subRegion, m_dofManager, m_rhs.values() );
      if( m_nonlinearSolverParameters.getLogLevel() > 4 )
        GEOS_LOG_LEVEL_RANK_0( logInfo::Convergence,
                               GEOS_FMT( "        ( R ) = ( {:4.2e} )", residualNorm ) );
    }
    //auto iterInfo = currentIter( time_n, dt );
    //outputSingleWellDebug( time_n, stepDt, 0, newtonIter, 0,
    //                       mesh, subRegion, dofManager, m_localMatrix.toViewConstSizes(), m_rhs.values()  );
    // if the residual norm is less than the Newton tolerance we denote that we have
    // converged and break from the Newton loop immediately.
    if( residualNorm < newtonTol && newtonIter >= minNewtonIter )
    {
      isNewtonConverged = true;
      break;
    }

// if the residual norm is above the max allowed residual norm, we break from
// the Newton loop to avoid crashes due to Newton divergence
    if( residualNorm > m_nonlinearSolverParameters.m_maxAllowedResidualNorm )
    {
      string const maxAllowedResidualNormString = NonlinearSolverParameters::viewKeysStruct::maxAllowedResidualNormString();
      if( m_nonlinearSolverParameters.getLogLevel() > 4 )
        GEOS_LOG_LEVEL_RANK_0( logInfo::Convergence,
                               GEOS_FMT( "    The residual norm is above the {} of {}. Newton loop terminated.",
                                         maxAllowedResidualNormString,
                                         m_nonlinearSolverParameters.m_maxAllowedResidualNorm )  );
      isNewtonConverged = false;
      break;
    }

    {
      Timer timer( m_timers.get_inserted( "linear solver total" ) );

// TODO: Trilinos currently requires this, re-evaluate after moving to Tpetra-based solvers
      if( m_precond )
      {
        m_precond->clear();
      }

      {
        Timer timer_setup( m_timers.get_inserted( "linear solver create" ) );

// Compose parallel LA matrix/rhs out of local LA matrix/rhs
//
        m_matrix.create( m_localMatrix.toViewConst(), m_dofManager.numLocalDofs(), MPI_COMM_GEOS );
      }

// Output the linear system matrix/rhs for debugging purposes
      //string tag = "_"+std::to_string( my_ctime );  tjb
      //debugOutputSystem( time_n, cycleNumber, newtonIter, m_matrix, m_rhs, tag );

      debugOutputSystem( time_n, cycleNumber, newtonIter, m_matrix, m_rhs );
// Solve the linear system
      solveLinearSystem( m_dofManager, m_matrix, m_rhs, m_solution );

// Increment the solver statistics for reporting purposes
      getIterationStats().updateNonlinearIteration( m_linearSolverResult.numIterations );

// Output the linear system solution for debugging purposes
      debugOutputSolution( time_n, cycleNumber, newtonIter, m_solution );
      //debugOutputSolution( time_n, cycleNumber, newtonIter, m_solution, tag );
    }

    {
      Timer timer( m_timers.get_inserted( "apply solution" ) );

// Compute the scaling factor for the Newton update
      scaleFactor = well.scalingForWellSystemSolution( subRegion, m_dofManager, m_solution.values() );
      if( m_nonlinearSolverParameters.getLogLevel() > 4 )
        GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                               GEOS_FMT( "        {}: Global solution scaling factor = {}", getName(), scaleFactor ) );

      real64 minPressure = 0.0, minDensity = 0.0, minTotalDensity = 0.0;
      bool const solutionLogActive = isLogLevelActive< logInfo::Solution >( getLogLevel() );
      bool const solutionDetailsLogActive = isLogLevelActive< logInfo::SolutionDetails >( getLogLevel() );
      ElementsReporterBuffer rankNegPressureIds{ solutionLogActive, solutionDetailsLogActive ? 16 : 0 };
      ElementsReporterBuffer rankNegDensityIds{ solutionLogActive, solutionDetailsLogActive ? 16 : 0 };
      // output only total density sum, not cell details
      ElementsReporterBuffer rankTotalNegDensityIds{ solutionLogActive, 0 };
      if( !well.checkWellSystemSolution( subRegion, m_dofManager, m_solution.values(), scaleFactor, minPressure, minDensity, minTotalDensity, rankNegPressureIds, rankNegDensityIds,
                                         rankTotalNegDensityIds ) )
      {
// TODO try chopping (similar to line search)
        if( m_nonlinearSolverParameters.getLogLevel() > 4 )
          GEOS_LOG_RANK_0( GEOS_FMT( "    {}: Solution check failed. Newton loop terminated.", getName()) );
        break;
      }

// apply the system solution to the fields/variables
      well.applyWellSystemSolution( m_dofManager, m_solution.values(), scaleFactor, stepDt, domain, mesh, subRegion );
    }

    {
      Timer timer( m_timers.get_inserted( "update state" ) );

      // update derived variables (constitutive models)
      well.updateWellState( domain.getMeshBody( mesh.getParent().getParent().getName() ), elemManager, subRegion );
    }

  }
  return isNewtonConverged;
}


/**
 * @brief String for the stats output type
 */
ENUM_STRINGS( WellNewtonSolver::StatsOutputType,
              "none",
              "iteration",
              "convergence",
              "all" );

} // namespace geos


#endif /* GEOS_PHYSICSSOLVERS_WELLNEWTONSOLVER_HPP_ */
