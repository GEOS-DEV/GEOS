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

#ifndef GEOSX_PHYSICSSOLVERS_SOLVERBASE_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLVERBASE_HPP_

#include "codingUtilities/traits.hpp"
#include "common/DataTypes.hpp"
#include "dataRepository/ExecutableGroup.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/utilities/LinearSolverResult.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "physicsSolvers/NonlinearSolverParameters.hpp"
#include "physicsSolvers/LinearSolverParameters.hpp"


#include <limits>

namespace geosx
{

class DomainPartition;

class SolverBase : public ExecutableGroup
{
public:

  explicit SolverBase( string const & name,
                       Group * const parent );

  SolverBase( SolverBase && ) = default;

  virtual ~SolverBase() override;

  SolverBase() = delete;
  SolverBase( SolverBase const & ) = delete;
  SolverBase & operator=( SolverBase const & ) = delete;
  SolverBase & operator=( SolverBase && ) = delete;

  static string catalogName() { return "SolverBase"; }


  virtual void registerDataOnMesh( Group & MeshBodies ) override;


  virtual void initialize_postMeshGeneration() override;

  void generateMeshTargetsFromTargetRegions( Group const & meshBodies );

  /**
   * This method is called when its host event is triggered
   */
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**
   * @brief Getter for system matrix
   * @return a reference to linear system matrix of this solver
   */
  ParallelMatrix & getSystemMatrix() { return m_matrix; }
  ParallelMatrix const & getSystemMatrix() const { return m_matrix; }

  /**
   * @brief Getter for system rhs vector
   * @return a reference to linear system right-hand side of this solver
   */
  ParallelVector & getSystemRhs() { return m_rhs; }
  ParallelVector const & getSystemRhs() const { return m_rhs; }

  /**
   * @brief Getter for system solution vector
   * @return a reference to solution vector of this solver
   */
  ParallelVector & getSystemSolution() { return m_solution; }
  ParallelVector const & getSystemSolution() const { return m_solution; }

  /**
   * @brief Getter for degree-of-freedom manager
   * @return a reference to degree-of-freedom manager of this solver
   */
  DofManager & getDofManager() { return m_dofManager; }
  DofManager const & getDofManager() const { return m_dofManager; }

  /**
   * @brief Getter for local matrix
   * @return a reference to linear system matrix of this solver
   */
  CRSMatrix< real64, globalIndex > & getLocalMatrix() { return m_localMatrix; }
  CRSMatrixView< real64 const, globalIndex const > getLocalMatrix() const { return m_localMatrix.toViewConst(); }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  /**
   * @brief entry function to perform a solver step
   * @param time_n time at the beginning of the step
   * @param dt the perscribed timestep
   * @param cycleNumber the current cycle number
   * @param domain the domain object
   * @return return the timestep that was achieved during the step.
   *
   * This function is the entry point to perform a solver step. The choice of time integration
   * method is determined in this function, and the appropriate step function is called.
   */
  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain );



  /**
   * @brief entry function to perform a solver step
   * @param [in]  time_n time at the beginning of the step
   * @param [in]  dt the perscribed timestep
   * @param [out] return the timestep that was achieved during the step.
   *
   * T
   */
  virtual void setNextDt( real64 const & currentDt,
                          real64 & nextDt );

  /**
   * @brief entry function to perform a solver step
   * @param [in]  time_n time at the beginning of the step
   * @param [in]  dt the perscribed timestep
   * @param [out] return the timestep that was achieved during the step.
   *
   * T
   */
  void setNextDtBasedOnNewtonIter( real64 const & currentDt,
                                   real64 & nextDt );


  /**
   * @brief Entry function for an explicit time integration step
   * @param time_n time at the beginning of the step
   * @param dt the perscribed timestep
   * @param cycleNumber the current cycle number
   * @param domain the domain object
   * @return return the timestep that was achieved during the step.
   */
  virtual real64 explicitStep( real64 const & time_n,
                               real64 const & dt,
                               integer const cycleNumber,
                               DomainPartition & domain );

  /**
   * @brief Function for a nonlinear implicit integration step
   * @param time_n time at the beginning of the step
   * @param dt the perscribed timestep
   * @param cycleNumber the current cycle number
   * @param domain the domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   * @return return the timestep that was achieved during the step.
   *
   * This function implements a nonlinear newton method for implicit problems. It requires that the
   * other functions in the solver interface are implemented in the derived physics solver. The
   * nonlinear loop includes a simple line search algorithm, and will cut the timestep if
   * convergence is not achieved according to the parameters in linearSolverParameters member.
   */
  virtual real64 nonlinearImplicitStep( real64 const & time_n,
                                        real64 const & dt,
                                        integer const cycleNumber,
                                        DomainPartition & domain );

  /**
   * @brief Function to perform line search
   * @param time_n time at the beginning of the step
   * @param dt the perscribed timestep
   * @param cycleNumber the current cycle number
   * @param domain the domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   * @param lastResidual (in) target value below which to reduce residual norm, (out) achieved residual norm
   * @return return true if line search succeeded, false otherwise
   *
   * This function implements a nonlinear newton method for implicit problems. It requires that the
   * other functions in the solver interface are implemented in the derived physics solver. The
   * nonlinear loop includes a simple line search algorithm, and will cut the timestep if
   * convergence is not achieved according to the parameters in linearSolverParameters member.
   */
  virtual bool
  lineSearch( real64 const & time_n,
              real64 const & dt,
              integer const cycleNumber,
              DomainPartition & domain,
              DofManager const & dofManager,
              CRSMatrixView< real64, globalIndex const > const & localMatrix,
              ParallelVector & rhs,
              ParallelVector & solution,
              real64 const scaleFactor,
              real64 & lastResidual );

  /**
   * @brief Function for a linear implicit integration step
   * @param time_n time at the beginning of the step
   * @param dt the perscribed timestep
   * @param cycleNumber the current cycle number
   * @param domain the domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   * @return return the timestep that was achieved during the step.
   *
   * This function implements a single linear step. Similar to the nonlinear step, however it
   * assumes that the solution is achieved in a iteration. The use of this function requires that
   * the other functions in the solver interface are implemented in the derived physics solver. The
   * nonlinear loop includes a simple line search algorithm, and will cut the timestep if
   * convergence is not achieved according to the parameters in linearSolverParameters member.
   */
  virtual real64 linearImplicitStep( real64 const & time_n,
                                     real64 const & dt,
                                     integer const cycleNumber,
                                     DomainPartition & domain );

  /**
   * @brief function to perform setup for implicit timestep
   * @param time_n the time at the beginning of the step
   * @param dt the desired timestep
   * @param domain the domain partition
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   *
   * This function should contain any step level initialization required to perform an implicit
   * step.
   *
   * @note This function must be overridden in the derived physics solver in order to use an implict
   * solution method such as LinearImplicitStep() or NonlinearImplicitStep().
   */
  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain );

  /**
   * @brief Populate degree-of-freedom manager with fields relevant to this solver
   * @param dofManager degree-of-freedom manager associated with the linear system
   */
  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const;

  /**
   * @brief Set up the linear system (DOF indices and sparsity patterns)
   * @param domain the domain containing the mesh and fields
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   *
   * @note While the function is virtual, the base class implementation should be
   *       sufficient for most single-physics solvers.
   */
  virtual void
  setupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               ParallelVector & rhs,
               ParallelVector & solution,
               bool const setSparsity = true );

  /**
   * @brief function to assemble the linear system matrix and rhs
   * @param time the time at the beginning of the step
   * @param dt the desired timestep
   * @param domain the domain partition
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   * @return the residual for convergence evaluation
   *
   * This function assembles the residual and the jacobian of the residual wrt the primary
   * variables. In a stand alone physics solver, this function will fill a single block in the
   * block system. However the capability to query the block system structure for any coupled blocks
   * may be implemented to fill in off diagonal blocks of the system to enable coupling between
   * solvers.
   *
   * @note This function must be overridden in the derived physics solver in order to use an implict
   * solution method such as LinearImplicitStep() or NonlinearImplicitStep().
   */
  virtual void
  assembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs );

  /**
   * @brief apply boundary condition to system
   * @param domain the domain partition
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param time the time at the beginning of the step
   * @param dt the desired timestep
   *
   * This function applies all boundary conditions to the linear system. This is essentially a
   * completion of the system assembly, but is separated for use in coupled solvers.
   */
  virtual void
  applyBoundaryConditions( real64 const time,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs );

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
   * @brief calculate the norm of the global system residual
   * @param rhs the system right-hand side vector
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain partition
   * @return norm of the residual
   *
   * This function returns the norm of global residual vector, which is suitable for comparison with
   * a tolerance.
   */
  virtual real64
  calculateResidualNorm( DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs );

  /**
   * @brief function to apply a linear system solver to the assembled system.
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param solution the solver parameters to set the parameters of the linear system
   *
   * This function calls the linear solver package to perform a single linear solve on the block
   * system. The derived physics solver is required to specify the call, as no default is provided.
   *
   * @note This function must be overridden in the derived physics solver in order to use an implict
   * solution method such as LinearImplicitStep() or NonlinearImplicitStep().
   */
  virtual void
  solveSystem( DofManager const & dofManager,
               ParallelMatrix & matrix,
               ParallelVector & rhs,
               ParallelVector & solution );

  /**
   * @brief Function to check system solution for physical consistency and constraint violation
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param scalingFactor factor to scale the solution prior to application
   * @param objectManager the object manager that holds the fields we wish to apply the solution to
   * @return true if solution can be safely applied without violating physical constraints, false otherwise
   *
   * @note This function must be overridden in the derived physics solver in order to use an implict
   * solution method such as LinearImplicitStep() or NonlinearImplicitStep().
   *
   */
  virtual bool
  checkSystemSolution( DomainPartition const & domain,
                       DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor );

  /**
   * @brief Function to determine if the solution vector should be scaled back in order to maintain a known constraint.
   * @param[in] domain The domain partition.
   * @param[in] dofManager degree-of-freedom manager associated with the linear system
   * @param[in] solution the solution vector
   * @return The factor that should be used to scale the solution vector values when they are being applied.
   */
  virtual real64
  scalingForSystemSolution( DomainPartition const & domain,
                            DofManager const & dofManager,
                            arrayView1d< real64 const > const & localSolution );

  /**
   * @brief Function to apply the solution vector to the state
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param scalingFactor factor to scale the solution prior to application
   * @param objectManager the object manager that holds the fields we wish to apply the solution to
   *
   * This function performs 2 operations:
   * 1) extract the solution vector for the "blockSystem" parameter, and applies the
   *    contents of the solution vector to the primary variable field data,
   * 2) perform a synchronization of the primary field variable such that all ghosts are updated,
   *
   * The "scalingFactor" parameter allows for the scaled application of the solution vector. For
   * instance, a line search may apply a negative scaling factor to remove part of the previously
   * applied solution.
   *
   * @note This function must be overridden in the derived physics solver in order to use an implict
   * solution method such as LinearImplicitStep() or NonlinearImplicitStep().
   *
   */
  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain );

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param domain the domain containing the mesh and fields
   */
  virtual void updateState( DomainPartition & domain );

  /**
   * @brief reset state of physics back to the beginning of the step.
   * @param domain
   *
   * This function essentially abandons the step, and resets all variables back to thier values at
   * the beginning of the step. This is useful for cases where convergence was not achieved, and
   * a cut in timestep was required.
   *
   * @note This function must be overridden in the derived physics solver in order to use an implict
   * solution method such as LinearImplicitStep() or NonlinearImplicitStep().
   */
  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain );

  /**
   * @brief perform cleanup for implicit timestep
   * @param time_n the time at the beginning of the step
   * @param dt the desired timestep
   * @param domain the domain partition
   *
   * This function performs whatever tasks are required to complete an implicit timestep. For
   * example, the acceptance of the solution will occur during this step, and deallocation of
   * temporaries will be be performed in this function.
   *
   * @note This function must be overridden in the derived physics solver in order to use an implict
   * solution method such as LinearImplicitStep() or NonlinearImplicitStep().
   */
  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain );


  /*
   * Returns the requirement for the next time-step to the event executing the solver.
   */
  virtual real64 getTimestepRequest( real64 const GEOSX_UNUSED_PARAM( time ) ) override
  {return m_nextDt;};
  /**@}*/

  real64 GetTimestepRequest()
  {return m_nextDt;};

  virtual Group * createChild( string const & childKey, string const & childName ) override;

  using CatalogInterface = dataRepository::CatalogInterface< SolverBase, string const &, Group * const >;
  static CatalogInterface::CatalogType & getCatalog();

  struct viewKeyStruct
  {
    static constexpr char const * cflFactorString() { return "cflFactor"; }
    static constexpr char const * initialDtString() { return "initialDt"; }
    static constexpr char const * maxStableDtString() { return "maxStableDt"; }
    static constexpr char const * discretizationString() { return "discretization"; }
    static constexpr char const * targetRegionsString() { return "targetRegions"; }
    static constexpr char const * meshTargetsString() { return "meshTargets"; }

  };

  struct groupKeyStruct
  {
    static constexpr char const * linearSolverParametersString() { return "LinearSolverParameters"; }
    static constexpr char const * nonlinearSolverParametersString() { return "NonlinearSolverParameters"; }
  };


  /**
   * @brief return the value of the gravity vector specified in PhysicsSolverManager
   * @return the value of the gravity vector
   *
   * @note if the solver is instantiated outside of a simulation (for instance for a unit test)
   *       and therefore does not have a parent of type PhysicsSolverManager, this function returns
   *       {0.0,0.0,-9.81}
   */
  R1Tensor const gravityVector() const;

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
   * @brief Get position of a given region within solver's target region list
   * @param regionName the region name to find
   * @return index within target regions list
   */
  localIndex targetRegionIndex( string const & regionName ) const;


  template< typename LAMBDA >
  void forMeshTargets( Group const & meshBodies, LAMBDA && lambda ) const
  {
    for( auto const & target: m_meshTargets )
    {
      string const meshBodyName = target.first;
      arrayView1d< string const > const & regionNames = target.second.toViewConst();
      MeshBody const & meshBody = meshBodies.getGroup< MeshBody >( meshBodyName );
      meshBody.forMeshLevels( [&]( MeshLevel const & meshLevel )
      {
        lambda( meshBodyName, meshLevel, regionNames );
      } );
    }
  }

  template< typename LAMBDA >
  void forMeshTargets( Group & meshBodies, LAMBDA && lambda ) const
  {
    for( auto const & target: m_meshTargets )
    {
      string const meshBodyName = target.first;
      arrayView1d< string const > const & regionNames = target.second.toViewConst();
      MeshBody & meshBody = meshBodies.getGroup< MeshBody >( meshBodyName );
      meshBody.forMeshLevels( [&]( MeshLevel & meshLevel )
      {
        lambda( meshBodyName, meshLevel, regionNames );
      } );
    }
  }

  template< typename LAMBDA >
  void forDiscretizationOnMeshTargets( Group const & meshBodies, LAMBDA && lambda ) const
  {
    for( auto const & target: m_meshTargets )
    {
      string const meshBodyName = target.first;
      arrayView1d< string const > const & regionNames = target.second.toViewConst();
      MeshBody const & meshBody = meshBodies.getGroup< MeshBody >( meshBodyName );
      MeshLevel const & meshLevel = meshBody.getMeshLevel(this->m_discretizationName);
      lambda( meshBodyName, meshLevel, regionNames );
    }
  }

  template< typename LAMBDA >
  void forDiscretizationOnMeshTargets( Group & meshBodies, LAMBDA && lambda ) const
  {
    for( auto const & target: m_meshTargets )
    {
      string const meshBodyName = target.first;
      arrayView1d< string const > const & regionNames = target.second.toViewConst();
      MeshBody & meshBody = meshBodies.getGroup< MeshBody >( meshBodyName );
      MeshLevel & meshLevel = meshBody.getMeshLevel(this->m_discretizationName);
      lambda( meshBodyName, meshLevel, regionNames );
    }
  }

  string getDiscretizationName() const {return m_discretizationName;}

  virtual bool registerCallback( void * func, const std::type_info & funcType ) final override;

  /**
   * @brief Return PySolver type.
   * @return Return PySolver type.
   */
#if defined(GEOSX_USE_PYGEOSX)
  virtual PyTypeObject * getPythonType() const override;
#endif

protected:


  static real64 eisenstatWalker( real64 const newNewtonNorm,
                                 real64 const oldNewtonNorm,
                                 real64 const weakestTol );

  /**
   * @brief Get the Constitutive Name object
   *
   * @tparam CONSTITUTIVE_BASE_TYPE the base type of the constitutive model.
   * @param subregion the element subregion on which the constitutive model is registered
   * @return the name name of the constitutive model of type @p CONSTITUTIVE_BASE_TYPE registered on the @p subregion.
   */
  template< typename CONSTITUTIVE_BASE_TYPE >
  static string getConstitutiveName( ElementSubRegionBase const & subRegion );

  /**
   * @brief This function sets constitutive name fields on an
   *  ElementSubRegionBase, and calls the base function it overrides.
   * @param subRegion The ElementSubRegionBase that will have constitutive
   *  names set.
   */
  virtual void setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const { GEOSX_UNUSED_VAR( subRegion ); }

  template< typename BASETYPE = constitutive::ConstitutiveBase, typename LOOKUP_TYPE >
  static BASETYPE const & getConstitutiveModel( dataRepository::Group const & dataGroup, LOOKUP_TYPE const & key );

  template< typename BASETYPE = constitutive::ConstitutiveBase, typename LOOKUP_TYPE >
  static BASETYPE & getConstitutiveModel( dataRepository::Group & dataGroup, LOOKUP_TYPE const & key );

  real64 m_cflFactor;
  real64 m_maxStableDt;
  real64 m_nextDt;

  /// name of the FV discretization object in the data repository
  string m_discretizationName;

  /// Data structure to handle degrees of freedom
  DofManager m_dofManager;

  /// System matrix, rhs and solution
  ParallelMatrix m_matrix;
  ParallelVector m_rhs;
  ParallelVector m_solution;

  /// Local system matrix and rhs
  CRSMatrix< real64, globalIndex > m_localMatrix;

  /// Custom preconditioner for the "native" iterative solver
  std::unique_ptr< PreconditionerBase< LAInterface > > m_precond;

  /// Linear solver parameters
  LinearSolverParametersInput m_linearSolverParameters;

  /// Result of the last linear solve
  LinearSolverResult m_linearSolverResult;

  /// Nonlinear solver parameters
  NonlinearSolverParameters m_nonlinearSolverParameters;

  std::function< void( CRSMatrix< real64, globalIndex >, array1d< real64 > ) > m_assemblyCallback;

  /// Map containing the array of target regions (value) for each MeshBody (key).
  map< string, array1d< string > > m_meshTargets;
  /// List of names of regions the solver will be applied to
  array1d< string > m_targetRegionNames;

private:

  /**
   * @brief This function sets constitutive name fields on an
   *  ElementSubRegionBase, and DOES NOT call the base function it overrides.
   * @param subRegion The ElementSubRegionBase that will have constitutive
   *  names set.
   */
  virtual void setConstitutiveNames( ElementSubRegionBase & subRegion ) const { GEOSX_UNUSED_VAR( subRegion ); }


};

template< typename CONSTITUTIVE_BASE_TYPE >
string SolverBase::getConstitutiveName( ElementSubRegionBase const & subRegion )
{
  string validName;
  dataRepository::Group const & constitutiveModels = subRegion.getConstitutiveModels();

  constitutiveModels.forSubGroups< CONSTITUTIVE_BASE_TYPE >( [&]( dataRepository::Group const & model )
  {
    GEOSX_ERROR_IF( !validName.empty(), "A valid constitutive model was already found." );
    validName = model.getName();
  } );
  return validName;
}

template< typename BASETYPE, typename LOOKUP_TYPE >
BASETYPE const & SolverBase::getConstitutiveModel( dataRepository::Group const & dataGroup, LOOKUP_TYPE const & key )
{
  dataRepository::Group const & constitutiveModels =
    dataGroup.getGroup( constitutive::ConstitutiveManager::groupKeyStruct::constitutiveModelsString() );

  return constitutiveModels.getGroup< BASETYPE >( key );
}

template< typename BASETYPE, typename LOOKUP_TYPE >
BASETYPE & SolverBase::getConstitutiveModel( dataRepository::Group & dataGroup, LOOKUP_TYPE const & key )
{
  Group & constitutiveModels =
    dataGroup.getGroup( constitutive::ConstitutiveManager::groupKeyStruct::constitutiveModelsString() );

  return constitutiveModels.getGroup< BASETYPE >( key );
}

} // namespace geosx


#endif /* GEOSX_PHYSICSSOLVERS_SOLVERBASE_HPP_ */
