/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef SOLVERBASE_HPP_
#define SOLVERBASE_HPP_



#include <string>
#include <limits>

#include "dataRepository/Group.hpp"
#include "codingUtilities/GeosxTraits.hpp"
#include "common/DataTypes.hpp"
#include "dataRepository/ExecutableGroup.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "physicsSolvers/SystemSolverParameters.hpp"

#include "linearAlgebraInterface/src/InterfaceTypes.hpp" // LAI
#include "linearAlgebraInterface/src/DofManager.hpp"

namespace geosx
{

class DomainPartition;
class SystemSolverParameters;

namespace dataRepository
{
namespace keys
{
string const courant = "courant";
string const maxDt   = "maxDt";
}
}

class SolverBase : public ExecutableGroup
{
public:

  explicit SolverBase( std::string const & name,
                       Group * const parent );

  SolverBase( SolverBase && ) = default;

  virtual ~SolverBase() override;

  SolverBase() = delete;
  SolverBase( SolverBase const & ) = delete;
  SolverBase& operator=( SolverBase const & ) = delete;
  SolverBase& operator=( SolverBase&& ) = delete;

  static string CatalogName() { return "SolverBase"; }

//  virtual void Registration( dataRepository::WrapperCollection& domain );



  /**
   * This method is called when it's host event is triggered
   */
  virtual void Execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        dataRepository::Group * const domain ) override;

  /**
   * @brief Getter for system matrix
   * @return a reference to linear system matrix of this solver
   */
  ParallelMatrix       & getSystemMatrix()       { return m_matrix; }
  ParallelMatrix const & getSystemMatrix() const { return m_matrix; }

  /**
   * @brief Getter for system rhs vector
   * @return a reference to linear system right-hand side of this solver
   */
  ParallelVector       & getSystemRhs()       { return m_rhs; }
  ParallelVector const & getSystemRhs() const { return m_rhs; }

  /**
   * @brief Getter for system solution vector
   * @return a reference to solution vector of this solver
   */
  ParallelVector       & getSystemSolution()       { return m_solution; }
  ParallelVector const & getSystemSolution() const { return m_solution; }

  /**
   * @brief Getter for degree-of-freedom manager
   * @return a reference to degree-of-freedom manager of this solver
   */
  DofManager       & getDofManager()       { return m_dofManager; }
  DofManager const & getDofManager() const { return m_dofManager; }

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
  virtual real64 SolverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition * const domain );

  /**
   * @brief Entry function for an explicit time integration step
   * @param time_n time at the beginning of the step
   * @param dt the perscribed timestep
   * @param cycleNumber the current cycle number
   * @param domain the domain object
   * @return return the timestep that was achieved during the step.
   */
  virtual real64 ExplicitStep( real64 const & time_n,
                               real64 const & dt,
                               integer const cycleNumber,
                               DomainPartition * const domain );

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
   * convergence is not achieved according to the parameters in systemSolverParameters member.
   */
  virtual real64 NonlinearImplicitStep( real64 const & time_n,
                                        real64 const & dt,
                                        integer const cycleNumber,
                                        DomainPartition * const domain,
                                        DofManager const & dofManager,
                                        ParallelMatrix & matrix,
                                        ParallelVector & rhs,
                                        ParallelVector & solution );

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
   * convergence is not achieved according to the parameters in systemSolverParameters member.
   */
  virtual bool
  LineSearch( real64 const & time_n,
              real64 const & dt,
              integer const cycleNumber,
              DomainPartition * const domain,
              DofManager const & dofManager,
              ParallelMatrix & matrix,
              ParallelVector & rhs,
              ParallelVector & solution,
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
   * convergence is not achieved according to the parameters in systemSolverParameters member.
   */
  virtual real64 LinearImplicitStep( real64 const & time_n,
                                     real64 const & dt,
                                     integer const cycleNumber,
                                     DomainPartition * const domain,
                                     DofManager & dofManager,
                                     ParallelMatrix & matrix,
                                     ParallelVector & rhs,
                                     ParallelVector & solution );

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
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition * const domain,
                     DofManager & dofManager,
                     ParallelMatrix & matrix,
                     ParallelVector & rhs,
                     ParallelVector & solution );

  /**
   * @brief Populate degree-of-freedom manager with fields relevant to this solver
   * @param dofManager degree-of-freedom manager associated with the linear system
   */
  virtual void
  SetupDofs( DomainPartition const * const domain,
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
  SetupSystem( DomainPartition * const domain,
               DofManager & dofManager,
               ParallelMatrix & matrix,
               ParallelVector & rhs,
               ParallelVector & solution );

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
  AssembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition * const domain,
                  DofManager const & dofManager,
                  ParallelMatrix & matrix,
                  ParallelVector & rhs );

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
  ApplyBoundaryConditions( real64 const time,
                           real64 const dt,
                           DomainPartition * const domain,
                           DofManager const & dofManager,
                           ParallelMatrix & matrix,
                           ParallelVector & rhs );

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
  CalculateResidualNorm( DomainPartition const * const domain,
                         DofManager const & dofManager,
                         ParallelVector const & rhs );

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
  SolveSystem( DofManager const & dofManager,
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
  CheckSystemSolution( DomainPartition const * const domain,
                       DofManager const & dofManager,
                       ParallelVector const & solution,
                       real64 const scalingFactor );

  /**
   * @brief Function to apply the solution vector to the state
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param scalingFactor factor to scale the solution prior to application
   * @param objectManager the object manager that holds the fields we wish to apply the solution to
   *
   * This function performs 3 operations:
   * 1) extract the solution vector for the "blockSystem" parameter, and applies the
   *    contents of the solution vector to the primary variable field data,
   * 2) perform a synchronization of the primary field variable such that all ghosts are updated,
   * 3) update secondary variables/state to ensure consistency.
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
  ApplySystemSolution( DofManager const & dofManager,
                       ParallelVector const & solution,
                       real64 const scalingFactor,
                       DomainPartition * const domain );

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
  ResetStateToBeginningOfStep( DomainPartition * const domain );

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
  ImplicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition * const domain );

  /**@}*/


  virtual Group * CreateChild( string const & childKey, string const & childName ) override;
  virtual void ExpandObjectCatalogs() override;

  using CatalogInterface = cxx_utilities::CatalogInterface< SolverBase, std::string const &, Group * const >;
  static CatalogInterface::CatalogType& GetCatalog();

  struct viewKeyStruct
  {
    constexpr static auto verboseLevelString = "verboseLevel";
    constexpr static auto gravityVectorString = "gravityVector";
    constexpr static auto cflFactorString = "cflFactor";
    constexpr static auto maxStableDtString = "maxStableDt";
    static constexpr auto discretizationString = "discretization";
    constexpr static auto targetRegionsString = "targetRegions";

  } viewKeys;

  struct groupKeyStruct
  {
    constexpr static auto systemSolverParametersString = "SystemSolverParameters";
  } groupKeys;



  R1Tensor const & getGravityVector() const { return m_gravityVector; }
  R1Tensor       & getGravityVector()       { return m_gravityVector; }
  R1Tensor const * globalGravityVector() const;

  integer verboseLevel() const { return m_verboseLevel; }

  /**
   * accessor for the system solver parameters.
   * @return
   */

  SystemSolverParameters * getSystemSolverParameters()
  {
    return &m_systemSolverParameters;
  }

  SystemSolverParameters const * getSystemSolverParameters() const
  {
    return &m_systemSolverParameters;
  }

  string getDiscretization() const {return m_discretizationName;}

  string_array const & getTargetRegions() const {return m_targetRegions;}


  template<bool CONST>
  using SubregionFunc = std::function<void ( add_const_if_t<ElementSubRegionBase, CONST> * )>;

  template<typename MESH, typename LAMBDA>
  typename std::enable_if<std::is_same<typename std::remove_cv<MESH>::type, MeshLevel>::value &&
                          std::is_convertible<LAMBDA, SubregionFunc<std::is_const<MESH>::value>>::value,
                          void>::type
  applyToSubRegions( MESH * const mesh, LAMBDA && lambda ) const
  {
    mesh->getElemManager()->forElementSubRegions( m_targetRegions, std::forward<LAMBDA>(lambda) );
  }

  template<bool CONST>
  using SubregionFuncComplete = std::function<void ( localIndex, localIndex,
                                                     add_const_if_t<ElementRegionBase, CONST> *,
                                                     add_const_if_t<ElementSubRegionBase, CONST> * )>;

  template<typename MESH, typename LAMBDA>
  typename std::enable_if<std::is_same<typename std::remove_cv<MESH>::type, MeshLevel>::value &&
                          std::is_convertible<LAMBDA, SubregionFuncComplete<std::is_const<MESH>::value>>::value,
                          void>::type
  applyToSubRegions( MESH * const mesh, LAMBDA && lambda ) const
  {
    mesh->getElemManager()->forElementSubRegionsComplete( m_targetRegions, std::forward<LAMBDA>(lambda) );
  }

protected:

  void PostProcessInput() override;

  void SetLinearSolverParameters();

  string getDiscretizationName() const {return m_discretizationName;}

  template<typename BASETYPE>
  static BASETYPE const * GetConstitutiveModel( dataRepository::Group const * dataGroup, string const & name );

  template<typename BASETYPE>
  static BASETYPE * GetConstitutiveModel( dataRepository::Group * dataGroup, string const & name );

  integer m_verboseLevel = 0;
  R1Tensor m_gravityVector;
  SystemSolverParameters m_systemSolverParameters;

  real64 m_cflFactor;
  real64 m_maxStableDt;

  /// name of the FV discretization object in the data repository
  string m_discretizationName;

  string_array m_targetRegions;


  /// Data structure to handle degrees of freedom
  DofManager m_dofManager;

  /// System matrix, rhs and solution
  ParallelMatrix m_matrix;
  ParallelVector m_rhs;
  ParallelVector m_solution;

  /// Linear solver parameters
  LinearSolverParameters m_linearSolverParameters;

};

template<typename BASETYPE>
BASETYPE const * SolverBase::GetConstitutiveModel( dataRepository::Group const * dataGroup, string const & name )
{
  Group const * const constitutiveModels =
    dataGroup->GetGroup( constitutive::ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOS_ERROR_IF( constitutiveModels == nullptr, "Target group does not contain constitutive models" );

  BASETYPE const * const model = constitutiveModels->GetGroup<BASETYPE>( name );
  GEOS_ERROR_IF( model == nullptr, "Target group does not contain model " << name );

  return model;
}

template<typename BASETYPE>
BASETYPE * SolverBase::GetConstitutiveModel( dataRepository::Group * dataGroup, string const & name )
{
  Group * const constitutiveModels =
    dataGroup->GetGroup( constitutive::ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOS_ERROR_IF( constitutiveModels == nullptr, "Target group does not contain constitutive models" );

  BASETYPE * const model = constitutiveModels->GetGroup<BASETYPE>( name );
  GEOS_ERROR_IF( model == nullptr, "Target group does not contain model " << name );

  return model;
}

} /* namespace ANST */


#endif /* SOLVERBASE_HPP_ */
