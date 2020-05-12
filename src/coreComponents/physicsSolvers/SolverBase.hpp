/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLVERBASE_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLVERBASE_HPP_



#include <string>
#include <limits>

#include "dataRepository/Group.hpp"
#include "codingUtilities/traits.hpp"
#include "common/DataTypes.hpp"
#include "dataRepository/ExecutableGroup.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "physicsSolvers/SystemSolverParameters.hpp"
#include "physicsSolvers/NonlinearSolverParameters.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/DofManager.hpp"

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
  SolverBase & operator=( SolverBase const & ) = delete;
  SolverBase & operator=( SolverBase && ) = delete;

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
  ParallelMatrix & getSystemMatrix()       { return m_matrix; }
  ParallelMatrix const & getSystemMatrix() const { return m_matrix; }

  /**
   * @brief Getter for system rhs vector
   * @return a reference to linear system right-hand side of this solver
   */
  ParallelVector & getSystemRhs()       { return m_rhs; }
  ParallelVector const & getSystemRhs() const { return m_rhs; }

  /**
   * @brief Getter for system solution vector
   * @return a reference to solution vector of this solver
   */
  ParallelVector & getSystemSolution()       { return m_solution; }
  ParallelVector const & getSystemSolution() const { return m_solution; }

  /**
   * @brief Getter for degree-of-freedom manager
   * @return a reference to degree-of-freedom manager of this solver
   */
  DofManager & getDofManager()       { return m_dofManager; }
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
   * @brief entry function to perform a solver step
   * @param [in]  time_n time at the beginning of the step
   * @param [in]  dt the perscribed timestep
   * @param [out] return the timestep that was achieved during the step.
   *
   * T
   */
  virtual void SetNextDt( real64 const & currentDt,
                          real64 & nextDt );

  /**
   * @brief entry function to perform a solver step
   * @param [in]  time_n time at the beginning of the step
   * @param [in]  dt the perscribed timestep
   * @param [out] return the timestep that was achieved during the step.
   *
   * T
   */
  void SetNextDtBasedOnNewtonIter( real64 const & currentDt,
                                   real64 & nextDt );


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
              ParallelVector const & solution,
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
   * @brief Function to determine if the solution vector should be scaled back in order to maintain a known constraint.
   * @param[in] domain The domain partition.
   * @param[in] dofManager degree-of-freedom manager associated with the linear system
   * @param[in] solution the solution vector
   * @return The factor that should be used to scale the solution vector values when they are being applied.
   */
  virtual real64
  ScalingForSystemSolution( DomainPartition const * const domain,
                            DofManager const & dofManager,
                            ParallelVector const & solution );

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


  /*
   * Returns the requirement for the next time-step to the event executing the solver.
   */
  virtual real64 GetTimestepRequest( real64 const GEOSX_UNUSED_PARAM( time ) ) override
  {return m_nextDt;};
  /**@}*/

  real64 GetTimestepRequest()
  {return m_nextDt;};

  virtual Group * CreateChild( string const & childKey, string const & childName ) override;
  virtual void ExpandObjectCatalogs() override;

  using CatalogInterface = dataRepository::CatalogInterface< SolverBase, std::string const &, Group * const >;
  static CatalogInterface::CatalogType & GetCatalog();

  struct viewKeyStruct
  {
    constexpr static auto cflFactorString = "cflFactor";
    constexpr static auto initialDtString = "initialDt";
    constexpr static auto maxStableDtString = "maxStableDt";
    static constexpr auto discretizationString = "discretization";
    constexpr static auto targetRegionsString = "targetRegions";

  } viewKeys;

  struct groupKeyStruct
  {
    constexpr static auto systemSolverParametersString = "SystemSolverParameters";
    constexpr static auto nonlinearSolverParametersString = "NonlinearSolverParameters";
  } groupKeys;


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


  NonlinearSolverParameters & getNonlinearSolverParameters()
  {
    return m_nonlinearSolverParameters;
  }

  NonlinearSolverParameters const & getNonlinearSolverParameters() const
  {
    return m_nonlinearSolverParameters;
  }

  string getDiscretization() const { return m_discretizationName; }

  arrayView1d< string const > const & targetRegionNames() const { return m_targetRegionNames; }

  /**
   * @brief Get position of a given region within solver's target region list
   * @param regionName the region name to find
   * @return index within target regions list
   */
  localIndex targetRegionIndex( string const & regionName ) const;

  template< typename REGIONTYPE = ElementRegionBase, typename ... REGIONTYPES, typename LAMBDA >
  void forTargetRegions( MeshLevel const & mesh, LAMBDA && lambda ) const
  {
    mesh.getElemManager()->
      template forElementRegions< REGIONTYPE, REGIONTYPES... >( targetRegionNames(), std::forward< LAMBDA >( lambda ) );
  }

  template< typename REGIONTYPE = ElementRegionBase, typename ... REGIONTYPES, typename LAMBDA >
  void forTargetRegions( MeshLevel & mesh, LAMBDA && lambda ) const
  {
    mesh.getElemManager()->
      template forElementRegions< REGIONTYPE, REGIONTYPES... >( targetRegionNames(), std::forward< LAMBDA >( lambda ) );
  }

  template< typename REGIONTYPE = ElementRegionBase, typename ... REGIONTYPES, typename LAMBDA >
  void forTargetRegionsComplete( MeshLevel const & mesh, LAMBDA && lambda ) const
  {
    mesh.getElemManager()->
      template forElementRegionsComplete< REGIONTYPE, REGIONTYPES... >( targetRegionNames(), std::forward< LAMBDA >( lambda ) );
  }

  template< typename REGIONTYPE = ElementRegionBase, typename ... REGIONTYPES, typename LAMBDA >
  void forTargetRegionsComplete( MeshLevel & mesh, LAMBDA && lambda ) const
  {
    mesh.getElemManager()->
      template forElementRegionsComplete< REGIONTYPE, REGIONTYPES... >( targetRegionNames(), std::forward< LAMBDA >( lambda ) );
  }

  template< typename SUBREGIONTYPE = ElementSubRegionBase, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forTargetSubRegions( MeshLevel const & mesh, LAMBDA && lambda ) const
  {
    mesh.getElemManager()->
      template forElementSubRegions< SUBREGIONTYPE, SUBREGIONTYPES... >( targetRegionNames(), std::forward< LAMBDA >( lambda ) );
  }

  template< typename SUBREGIONTYPE = ElementSubRegionBase, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forTargetSubRegions( MeshLevel & mesh, LAMBDA && lambda ) const
  {
    mesh.getElemManager()->
      template forElementSubRegions< SUBREGIONTYPE, SUBREGIONTYPES... >( targetRegionNames(), std::forward< LAMBDA >( lambda ) );
  }

  template< typename SUBREGIONTYPE = ElementSubRegionBase, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forTargetSubRegionsComplete( MeshLevel const & mesh, LAMBDA && lambda ) const
  {
    mesh.getElemManager()->
      template forElementSubRegionsComplete< SUBREGIONTYPE, SUBREGIONTYPES... >( targetRegionNames(), std::forward< LAMBDA >( lambda ) );
  }

  template< typename SUBREGIONTYPE = ElementSubRegionBase, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forTargetSubRegionsComplete( MeshLevel & mesh, LAMBDA && lambda ) const
  {
    mesh.getElemManager()->
      template forElementSubRegionsComplete< SUBREGIONTYPE, SUBREGIONTYPES... >( targetRegionNames(), std::forward< LAMBDA >( lambda ) );
  }

protected:

  virtual void PostProcessInput() override;

  void SetLinearSolverParameters();

  string getDiscretizationName() const {return m_discretizationName;}

  template< typename BASETYPE = constitutive::ConstitutiveBase, typename LOOKUP_TYPE >
  static BASETYPE const & GetConstitutiveModel( dataRepository::Group const & dataGroup, LOOKUP_TYPE const & key );

  template< typename BASETYPE = constitutive::ConstitutiveBase, typename LOOKUP_TYPE >
  static BASETYPE & GetConstitutiveModel( dataRepository::Group & dataGroup, LOOKUP_TYPE const & key );

  /**
   * @brief Partially validates constitutive model names input.
   * @param[in,out] modelNames reference to input array of model names
   * @param[in] allowEmpty if @p true, empty array is not considered an error
   * @return flag indicating whether at least one model has been provided
   *
   * Checks that number of model names is equal to the number of solver's target regions.
   * Additionally, currently admits a single-element list, which is interpreted as one model
   * used for all target regions (the list is resized and populated accordingly).
   * If @p allowEmpty is true and the input is empty, returns false, which the solver can
   * interpret as a signal this type of model is disabled for the run (for optional models).
   */
  bool CheckModelNames( array1d< string > & modelNames,
                        string const & attribute,
                        bool const allowEmpty = false ) const;

  /**
   * @brief Populate array of constitutive model indices from list of model names.
   * @tparam MODEL_TYPE Base class of constitutive models to check against
   * @param elemRegionManager reference to element manager
   * @param modelNames list of model names
   *
   * This function is typically called from solver's InitializePreSubGroups() method,
   * after constitutive models have been set up but before they are used.
   * Looks up each model by name and type in each subregion of target regions.
   */
  template< typename MODEL_TYPE = constitutive::ConstitutiveBase >
  void ValidateModelMapping( ElementRegionManager const & elemRegionManager,
                             arrayView1d< string const > const & modelNames ) const;

  SystemSolverParameters m_systemSolverParameters;

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

  /// Linear solver parameters
  LinearSolverParameters m_linearSolverParameters;

  /// Nonlinear solver parameters
  NonlinearSolverParameters m_nonlinearSolverParameters;

private:

  /// List of names of regions the solver will be applied to
  array1d< string > m_targetRegionNames;

};

template< typename BASETYPE, typename LOOKUP_TYPE >
BASETYPE const & SolverBase::GetConstitutiveModel( dataRepository::Group const & dataGroup, LOOKUP_TYPE const & key )
{
  Group const * const constitutiveModels =
    dataGroup.GetGroup( constitutive::ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOSX_ERROR_IF( constitutiveModels == nullptr, "Target group does not contain constitutive models" );

  BASETYPE const * const model = constitutiveModels->GetGroup< BASETYPE >( key );
  GEOSX_ERROR_IF( model == nullptr, "Target group does not contain model " << key );

  return *model;
}

template< typename BASETYPE, typename LOOKUP_TYPE >
BASETYPE & SolverBase::GetConstitutiveModel( dataRepository::Group & dataGroup, LOOKUP_TYPE const & key )
{
  Group * const constitutiveModels =
    dataGroup.GetGroup( constitutive::ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOSX_ERROR_IF( constitutiveModels == nullptr, "Target group does not contain constitutive models" );

  BASETYPE * const model = constitutiveModels->GetGroup< BASETYPE >( key );
  GEOSX_ERROR_IF( model == nullptr, "Target group does not contain model " << key );

  return *model;
}

template< typename MODEL_TYPE >
void SolverBase::ValidateModelMapping( ElementRegionManager const & elemRegionManager,
                                       arrayView1d< string const > const & modelNames ) const
{
  GEOSX_ERROR_IF_NE( modelNames.size(), m_targetRegionNames.size() );
  for( localIndex k = 0; k < modelNames.size(); ++k )
  {
    ElementRegionBase const & region = *elemRegionManager.GetRegion( m_targetRegionNames[k] );
    for( localIndex esr = 0; esr < region.numSubRegions(); ++esr )
    {
      ElementSubRegionBase const & subRegion = *region.GetSubRegion( esr );
      MODEL_TYPE const * const model = subRegion.GetConstitutiveModels()->GetGroup< MODEL_TYPE >( modelNames[k] );
      GEOSX_ERROR_IF( model == nullptr,
                      getName() << ": constitutive model " << modelNames[k] << " not found in " << region.getName() << '/' << subRegion.getName() );
    }
  }
}

} // namespace geosx


#endif /* GEOSX_PHYSICSSOLVERS_SOLVERBASE_HPP_ */
