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
 * @file WellSolverBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLSOLVERBASE_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLSOLVERBASE_HPP_

#include "physicsSolvers/PhysicsSolverBase.hpp"
#include "physicsSolvers/fluidFlow/wells/WellPropWriter.hpp"

namespace geos
{

class DomainPartition;
class WellControls;
class WellElementSubRegion;

/**
 * @class WellSolverBase
 *
 * Base class for well solvers.
 * Provides some common features
 */
class WellSolverBase : public PhysicsSolverBase
{
public:

  /// String used to form the solverName used to register single-physics solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "well"; }

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  WellSolverBase( const string & name,
                  Group * const parent );

  /// default destructor
  virtual ~WellSolverBase() override;

  /// deleted default constructor
  WellSolverBase() = delete;

  /// deleted copy constructor
  WellSolverBase( WellSolverBase const & ) = delete;

  /// default move constructor
  WellSolverBase( WellSolverBase && ) = default;

  /// deleted assignment operator
  WellSolverBase & operator=( WellSolverBase const & ) = delete;

  /// deleted move operator
  WellSolverBase & operator=( WellSolverBase && ) = delete;

  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /// Expand catalog for schema generation
  virtual void expandObjectCatalogs() override;


  /**
   * @brief setter for the name of the flow solver (needed to use the flow kernels like UpdateFluid)
   * @param name the name of the flow solver
   */
  void setFlowSolverName( string const & name ) { m_flowSolverName = name; }

  /**
   * @brief getter for the name of the flow solver (used in UpdateState)
   * @return a string containing the name of the flow solver
   */
  string const & getFlowSolverName() const { return m_flowSolverName; }

  /**
   * @brief getter for the number of degrees of freedom per well element
   * @return the number of dofs
   */
  localIndex numDofPerWellElement() const { return m_numDofPerWellElement; }

  /**
   * @brief getter for the number of degrees of freedom per mesh element
   * @return the number of dofs
   */
  localIndex numDofPerResElement() const { return m_numDofPerResElement; }

  /**
   * @brief getter for iso/thermal switch
   * @return True if thermal
   */
  integer isThermal() const { return m_isThermal; }

  /**
   * @brief get the name of DOF defined on well elements
   * @return name of the DOF field used by derived solver type
   */
  virtual string wellElementDofName() const = 0;

  /**
   * @brief get the name of DOF defined on well elements
   * @return name of the DOF field used by derived solver type
   */
  virtual string resElementDofName() const = 0;

  /**
   * @brief const getter for the number of fluid components
   * @return the number of fluid components
   */
  virtual localIndex numFluidComponents() const = 0;

  /**
   * @brief Get the number of fluid phases
   * @return the number of phases
   */
  virtual localIndex numFluidPhases() const = 0;

  /**
   * @brief getter for the well controls associated to this well subRegion
   * @param subRegion the well subRegion whose controls are requested
   * @return a reference to the controls
   */
  WellControls & getWellControls( WellElementSubRegion const & subRegion );

  /**
   * @brief const getter for the well controls associated to this well subRegion
   * @param subRegion the well subRegion whose controls are requested
   * @return a reference to the const controls
   */
  WellControls const & getWellControls( WellElementSubRegion const & subRegion ) const;


  /**
   * @brief Open and close perfs based on user defined perf status table
   * @param time_n evaluation time
   * @param domain  the domain
   */
  void setPerforationStatus( real64 const & time_n, DomainPartition & domain );

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void registerDataOnMesh( Group & meshBodies ) override;

  void selectWellConstraint( real64 const & time_n,
                             real64 const & dt,
                             integer const coupledIterationNumber,
                             DomainPartition & domain );


  void setupWellDofs( DomainPartition & domain );

  void setupWellSystem ( DomainPartition & domain,
                         DofManager & dofManager,
                         CRSMatrix< real64, globalIndex > & localMatrix,
                         ParallelVector & rhs,
                         ParallelVector & solution,
                         bool const setSparsity = true );

  virtual void setupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override;

  virtual void implicitStepSetup( real64 const & time_n,
                                  real64 const & dt,
                                  DomainPartition & domain ) override;

  virtual void implicitStepComplete( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                     real64 const & GEOS_UNUSED_PARAM( dt ),
                                     DomainPartition & GEOS_UNUSED_PARAM( domain ) ) override {}

  virtual void applyBoundaryConditions( real64 const GEOS_UNUSED_PARAM( time_n ),
                                        real64 const GEOS_UNUSED_PARAM( dt ),
                                        DomainPartition & GEOS_UNUSED_PARAM( domain ),
                                        DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                        CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ),
                                        arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) ) override {}

  virtual void applyWellBoundaryConditions( real64 const GEOS_UNUSED_PARAM( time_n ),
                                            real64 const GEOS_UNUSED_PARAM( dt ),
                                            ElementRegionManager & GEOS_UNUSED_PARAM( elemManager ),
                                            WellElementSubRegion & GEOS_UNUSED_PARAM( subRegion ),
                                            DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                            arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ),
                                            CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ) )  {}

  /**@}*/

  /**
   * @brief function to assemble the linear system matrix and rhs
   * @param time the time at the beginning of the step
   * @param dt the desired timestep
   * @param domain the domain partition
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */

  void assembleWellSystem( real64 const time,
                           real64 const dt,
                           ElementRegionManager const & elementRegionManager,
                           WellElementSubRegion & subRegion,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs );

  virtual void assembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

  virtual void assembleWellFluxTerms( real64 const & GEOS_UNUSED_PARAM( time ),
                                      real64 const & GEOS_UNUSED_PARAM( dt ),
                                      WellElementSubRegion const & GEOS_UNUSED_PARAM( subRegion ),
                                      DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                      CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ),
                                      arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) ) = 0;


  /**
   * @brief assembles the flux terms for all connections between well elements
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void assembleFluxTerms( real64 const & time_n,
                                  real64 const & dt,
                                  DomainPartition & domain,
                                  DofManager const & dofManager,
                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                  arrayView1d< real64 > const & localRhs ) = 0;



  virtual void assembleWellAccumulationTerms( real64 const & time,
                                              real64 const & dt,
                                              WellElementSubRegion & subRegion,
                                              DofManager const & dofManager,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs ) = 0;


  /**
   * @brief assembles the accumulation term for all the well elements
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void assembleAccumulationTerms( real64 const & time_n,
                                          real64 const & dt,
                                          DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs ) = 0;

  virtual void assembleWellConstraintTerms( real64 const & GEOS_UNUSED_PARAM( time ),
                                            real64 const & GEOS_UNUSED_PARAM( dt ),
                                            WellElementSubRegion const & GEOS_UNUSED_PARAM( subRegion ),
                                            DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                            CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ),
                                            arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) ) = 0;

  virtual void assembleWellPressureRelations ( real64 const & GEOS_UNUSED_PARAM( time ),
                                               real64 const & GEOS_UNUSED_PARAM( dt ),
                                               WellElementSubRegion const & GEOS_UNUSED_PARAM( subRegion ),
                                               DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                               CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ),
                                               arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) ) = 0;
  /**
   * @brief assembles the pressure relations at all connections between well elements except at
   * the well head
   * @param time_n time at the beginning of the time step
   * @param dt the time step size
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void assemblePressureRelations( real64 const & time_n,
                                          real64 const & dt,
                                          DomainPartition const & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs ) = 0;

  virtual void outputSingleWellDebug( real64 const GEOS_UNUSED_PARAM( time ),
                                      real64 const GEOS_UNUSED_PARAM( dt ),
                                      integer GEOS_UNUSED_PARAM( num_timesteps ),
                                      integer GEOS_UNUSED_PARAM( current_newton_iteration ),
                                      integer GEOS_UNUSED_PARAM( num_timestep_cuts ),
                                      MeshLevel & GEOS_UNUSED_PARAM( mesh ),
                                      WellElementSubRegion & GEOS_UNUSED_PARAM( subRegion ),
                                      DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                      CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ),
                                      arrayView1d< const real64 > const & GEOS_UNUSED_PARAM( localRhs ) ) = 0;

  virtual void outputWellDebug( real64 const time,
                                real64 const dt,
                                integer num_timesteps,
                                integer current_newton_iteration,
                                integer num_timestep_cuts,
                                DomainPartition & domain,
                                DofManager const & dofManager,
                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                arrayView1d< real64 > const & localRhs ) = 0;
  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive
   * models)
   * @param elemManager the element region manager
   * @param subRegion the well subRegion containing the well elements and their associated fields
   */
  virtual real64 updateWellState( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion ) = 0;
  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param domain the domain containing the mesh and fields
   */
  virtual void updateState( DomainPartition & domain ) override;
  virtual void saveState( WellElementSubRegion & subRegion )   = 0;
  /**
   * @brief Initialize all the primary and secondary variables in all the wells
   * @param domain the domain containing the well manager to access individual wells
   */
  virtual void initializeWells( DomainPartition & domain, real64 const & time_n ) = 0;
  virtual void initializeWell( DomainPartition & domain, MeshLevel & mesh, WellElementSubRegion & subRegion, real64 const & time_n ) = 0;

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive
   * models)
   * @param elemManager the element region manager
   * fields
   */
  virtual real64 updateSubRegionState( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion ) = 0;


  virtual void computeWellPerforationRates( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                            real64 const & GEOS_UNUSED_PARAM( dt ),
                                            ElementRegionManager const & GEOS_UNUSED_PARAM( elemManager ),
                                            WellElementSubRegion & GEOS_UNUSED_PARAM( subRegion ) ){}

  /**
   * @brief Recompute the perforation rates for all the wells
   * @param domain the domain containing the mesh and fields
   */
  virtual void computePerforationRates( real64 const & time_n,
                                        real64 const & dt,
                                        DomainPartition & domain ) = 0;

  virtual real64
  calculateWellResidualNorm( real64 const & GEOS_UNUSED_PARAM( time_n ),
                             real64 const & GEOS_UNUSED_PARAM( dt ),
                             WellElementSubRegion const & GEOS_UNUSED_PARAM( subRegion ),
                             DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                             arrayView1d< real64 const > const & GEOS_UNUSED_PARAM( localRhs ) ) = 0;

  virtual real64
    scalingForWellSystemSolution( ElementSubRegionBase & GEOS_UNUSED_PARAM( subRegion ),
                                  DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                  arrayView1d< real64 const > const & GEOS_UNUSED_PARAM( localSolution ) ) = 0;

  virtual bool
    checkWellSystemSolution( ElementSubRegionBase & GEOS_UNUSED_PARAM( subRegion ),
                             DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                             arrayView1d< real64 const > const & GEOS_UNUSED_PARAM( localSolution ),
                             real64 const GEOS_UNUSED_PARAM( scalingFactor ) ) = 0;
  virtual void
    applyWellSystemSolution( DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                             arrayView1d< real64 const > const & GEOS_UNUSED_PARAM( localSolution ),
                             real64 const GEOS_UNUSED_PARAM( scalingFactor ),
                             real64 const GEOS_UNUSED_PARAM( dt ),
                             DomainPartition & GEOS_UNUSED_PARAM( domain ),
                             MeshLevel & GEOS_UNUSED_PARAM( mesh ),
                             WellElementSubRegion & GEOS_UNUSED_PARAM( subRegion ) ) = 0;

  bool solveNonlinearSystem( real64 const & time_n,
                             real64 const & stepDt,
                             integer const cycleNumber,
                             DomainPartition & domain,
                             MeshLevel & mesh,
                             ElementRegionManager & elementRegionManager,
                             WellElementSubRegion & subregion,
                             DofManager const & dofManager );

  /**
   * @brief Function to perform line search
   * @param time_n time at the beginning of the step
   * @param dt the perscribed timestep
   * @param cycleNumber the current cycle number
   * @param domain the domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   * @param scaleFactor the scaling factor to apply to the solution
   * @param lastResidual (in) target value below which to reduce residual norm, (out) achieved
   * residual norm
   * @return return true if line search succeeded, false otherwise
   *
   * This function implements a nonlinear newton method for implicit problems. It requires that
   * the
   * other functions in the solver interface are implemented in the derived physics solver. The
   * nonlinear loop includes a simple line search algorithm, and will cut the timestep if
   * convergence is not achieved according to the parameters in linearSolverParameters member.
   */
  bool
  lineSearch1( real64 const & time_n,
               real64 const & dt,
               integer const cycleNumber,
               DomainPartition & domain,
               ElementRegionManager & elemManager,
               WellElementSubRegion & subRegion,
               MeshLevel & mesh,
               DofManager const & dofManager,
               CRSMatrixView< real64, globalIndex const > const & localMatrix,
               ParallelVector & rhs,
               ParallelVector & solution,
               real64 const scaleFactor,
               real64 & lastResidual );
  /**
   * @brief function to set the next time step size
   * @param[in] currentTime the current time
   * @param[in] currentDt the current time step size
   * @param[in] domain the domain object
   * @return the prescribed time step size
   */
  virtual real64 setNextDt( real64 const & currentTime,
                            real64 const & currentDt,
                            DomainPartition & domain ) override;

  /**
   * @brief Utility function to keep the well variables during a time step (used in
   * poromechanics simulations)
   * @param[in] keepVariablesConstantDuringInitStep flag to tell the solver to freeze its
   * primary variables during a time step
   * @detail This function is meant to be called by a specific task before/after the
   * initialization step
   */
  void setKeepVariablesConstantDuringInitStep( bool const keepVariablesConstantDuringInitStep )
  { m_keepVariablesConstantDuringInitStep = keepVariablesConstantDuringInitStep; }

  struct viewKeyStruct : PhysicsSolverBase::viewKeyStruct
  {
    static constexpr char const * isThermalString() { return "isThermal"; }
    static constexpr char const * writeCSVFlagString() { return "writeCSV"; }
    static constexpr char const * timeStepFromTablesFlagString() { return "timeStepFromTables"; }
    static constexpr char const * writeSegDebugFlagString() { return "writeSegDebug"; }


    static constexpr char const * fluidNamesString() { return "fluidNames"; }
  };


  std::tuple< integer, integer, integer > currentIter( real64 const time, real64 const dt );
private:

  /**
   * @brief This function generates various discretization information for later use.
   * @param domain the domain parition
   */
  void precomputeData( DomainPartition & domain );

protected:

  virtual void postInputInitialization() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  virtual void initializePostSubGroups() override;


  /**
   * @brief Make sure that the well constraints are compatible
   * @param time_n the time at the beginning of the time step
   * @param dt the time step dt
   * @param subRegion the well subRegion
   */
  virtual void validateWellConstraints( real64 const & time_n,
                                        real64 const & dt,
                                        WellElementSubRegion const & subRegion ) = 0;

  virtual void printRates( real64 const & time_n,
                           real64 const & dt,
                           DomainPartition & domain ) = 0;

  virtual bool evaluateConstraints( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                    real64 const & GEOS_UNUSED_PARAM( stepDt ),
                                    integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                    integer const GEOS_UNUSED_PARAM( coupledIterationNumber ),
                                    DomainPartition & GEOS_UNUSED_PARAM( domain ),
                                    MeshLevel & GEOS_UNUSED_PARAM( mesh ),
                                    ElementRegionManager & GEOS_UNUSED_PARAM( elemManager ),
                                    WellElementSubRegion & GEOS_UNUSED_PARAM( subRegion ),
                                    DofManager const & GEOS_UNUSED_PARAM( dofManager ) ) { return false;};



  /// name of the flow solver
  string m_flowSolverName;

  /// the max number of fluid phases
  integer m_numPhases;

  /// the number of fluid components
  integer m_numComponents;

  /// the number of Degrees of Freedom per well element
  integer m_numDofPerWellElement;

  /// the number of Degrees of Freedom per reservoir element
  integer m_numDofPerResElement;

  /// flag indicating whether thermal formulation is used
  integer m_isThermal;

  /// rates output
  integer m_writeCSV;
  string const m_ratesOutputDir;

  // flag to enable time step selection base on rates/bhp tables coordinates
  integer m_timeStepFromTables;

  /// flag to freeze the initial state during initialization in coupled problems
  bool m_keepVariablesConstantDuringInitStep;
  /// flag to write detailed segment properties
  integer m_writeSegDebug;

  integer m_globalNumTimeSteps;
  real64 m_currentTime;
  real64 m_currentDt;
  real64 m_prevTime;
  real64 m_prevDt;
  integer m_numTimeStepCuts;
  integer m_currentNewtonIteration;

  std::map< std::string, WellPropWriter > m_wellPropWriter;
  std::map< std::string, WellPropWriter > m_wellPropWriter_eot;



  /// name of the fluid constitutive model used as a reference for component/phase description
  string m_referenceFluidModelName;

  /// flag to use the estimator
  integer m_estimateSolution;


  /// @brief  DofManagers for each wells estimator
  /// @details This DofManager is used to store the DOF numbers for the estimator
  /// @note This DofManager is used in the assembly of the estimators linear system
  std::map< std::string, DofManager >   m_estimatorDoFManager;
  integer my_ctime; //tjb

  real64 m_nextDt;
};

}

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLSOLVERBASE_HPP_
