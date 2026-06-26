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
 * @file WellManager.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELL_MANAGER_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELL_MANAGER_HPP_

#include "physicsSolvers/PhysicsSolverBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWell.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWell.hpp"
namespace geos
{

class DomainPartition;
class WellControls;
class WellElementSubRegion;

namespace dataRepository
{
namespace keys
{
static constexpr auto compositionalMultiphaseWell = "CompositionalMultiphaseWell";
static constexpr auto singlePhaseWell = "SinglePhaseWell";
}
}
/**
 * @class WellManager
 *
 * Base class for well solvers.
 * Provides some common features
 */
class WellManager : public PhysicsSolverBase
{
public:

  /// String used to form the solverName used to register single-physics solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "well"; }

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  WellManager( const string & name,
               Group * const parent );

  /// default destructor
  virtual ~WellManager() override = default;

  /// deleted default constructor
  WellManager() = delete;

  /// deleted copy constructor
  WellManager( WellManager const & ) = delete;

  /// default move constructor
  WellManager( WellManager && ) = default;

  /// deleted assignment operator
  WellManager & operator=( WellManager const & ) = delete;

  /// deleted move operator
  WellManager & operator=( WellManager && ) = delete;

  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /// Expand catalog for schema generation
  virtual void expandObjectCatalogs() override;

  /**
   * @brief setter for the name of the flow solver (needed to use the flow kernels like UpdateFluid)
   * @param name the name of the flow solver
   */
  void setFlowSolverName( string const & name ) { m_flowSolverName = name; }

  /**
   * @brief setter for compositional flag
   * @param compositional the compositional flag
   */
  void setCompositional( bool const & isCompositional ) { m_isCompositional = isCompositional; }

  /**
   * @brief getter for compositional flag
   * @return the compositional flag
   */
  bool isCompositional() const { return m_isCompositional; }


  /**
   * @brief getter for the name of the flow solver (used in UpdateState)
   * @return a string containing the name of the flow solver
   */
  string const & getFlowSolverName() const { return m_flowSolverName; }


  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName() { return "WellManager"; }
  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual void registerDataOnMesh( Group & meshBodies ) override;
  /**
   * @brief Get a well solver for a given well element sub-region
   * @param subRegion the well subRegion whose well solver is requested
   * @return a reference to the well solver
   */
  WellControls & getWell( WellElementSubRegion const & subRegion );

  /**
   * @brief Get a well solver for a given well element sub-region
   * @param wellControlsName name of well
   * @return a reference to the well solver
   */
  WellControls & getWell( std::string const & wellControlsName );

  /**
   * @brief Get a well solver for a given well element sub-region
   * @param wellControlsName name of well
   * @return a reference to the well solver
   */
  WellControls const & getWell( std::string const & wellControlsName ) const;
/**
 * @brief get the name of DOF defined on well elements
 * @return name of the DOF field used by derived solver type
 */
  string wellElementDofName() const;

  struct viewKeyStruct : PhysicsSolverBase::viewKeyStruct
  {
    static constexpr char const * dofFieldString() { return "wellVars"; }
    static constexpr char const * isThermalString() { return "isThermal"; }
    static constexpr char const * useMassFlagString() {return "useMass"; }
    /// @return string for the nextDt targetRegions wrapper
    static constexpr char const * targetRegionsString() { return "targetRegions"; }
    static constexpr char const * timeStepFromTablesFlagString() { return "timeStepFromTables"; }
    static constexpr char const * useTotalMassEquationString() { return "useTotalMassEquation"; }
    static constexpr char const * allowLocalCompDensChoppingString() { return CompositionalMultiphaseBase::viewKeyStruct::allowLocalCompDensChoppingString(); }


  };

  /**
   * @brief getter for the number of degrees of freedom per well element
   * @return the number of dofs
   */
  localIndex numDofPerWellElement() const;

  /**
   * @brief getter for the number of degrees of freedom per mesh element
   * @return the number of dofs
   */
  localIndex numDofPerResElement() const;

  /**
   * @brief getter for iso/thermal switch
   * @return True if thermal
   */
  integer isThermal() const;


  /**
   * @brief get the name of DOF defined on well elements
   * @return name of the DOF field used by derived solver type
   */
  virtual string resElementDofName() const;

  /**
   * @brief const getter for the number of fluid components
   * @return the number of fluid components
   */
  virtual localIndex numFluidComponents() const;

  /**
   * @brief const getter for the number of fluid phases
   * @return the number of fluid phases
   */
  virtual localIndex numFluidPhases() const;

  /**
   * @brief const getter for well total mass equation usage
   * @return true if total mass equation is used
   */
  integer useTotalMassEquation() const { return m_useTotalMassEquation; }

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
   * @brief getter for the compositional multiphase well associated to this well subRegion
   * @param subRegion the well subRegion whose controls are requested
   * @return a reference to the well
   */
  CompositionalMultiphaseWell & getCompositionalMultiphaseWell( WellElementSubRegion const & subRegion );

  /**
   * @brief const getter for the compositional multiphase well associated to this well subRegion
   * @param subRegion the well subRegion whose controls are requested
   * @return a reference to the const well
   */
  CompositionalMultiphaseWell const & getCompositionalMultiphaseWell( WellElementSubRegion const & subRegion ) const;

  /**
   * @brief Selects the active well constraint  based on current conditions
   * @param[in] currentTime the current time
   * @param[in] currentDt the current time step size
   * @param[in] coupledIterationNumber the current coupled iteration number
   * @param[in] domain the domain object
   * @return the prescribed time step size
   */
  void selectWellConstraint( real64 const & time_n,
                             real64 const & dt,
                             integer const coupledIterationNumber,
                             DomainPartition & domain );
  /* PhysicsSolverBase interfaces */

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

  virtual void setupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override;

  /**
   * @brief function to perform setup for implicit timestep
   * @param time_n the time at the beginning of the step
   * @param dt the desired timestep
   * @param domain the domain partition
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
                     DomainPartition & domain ) override;


  /**
   * @brief function to assemble the linear system matrix and rhs
   * @param time the time at the beginning of the step
   * @param dt the desired timestep
   * @param domain the domain partition
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
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
                  arrayView1d< real64 > const & localRhs ) override;


  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  virtual void applyBoundaryConditions( real64 const GEOS_UNUSED_PARAM( time_n ),
                                        real64 const GEOS_UNUSED_PARAM( dt ),
                                        DomainPartition & GEOS_UNUSED_PARAM( domain ),
                                        DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                        CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ),
                                        arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) ) override {}

  /**
   * @brief calculate the norm of the global system residual
   * @param time the time at the beginning of the step
   * @param dt the desired timestep
   * @param domain the domain partition
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localRhs the system right-hand side vector
   * @return norm of the residual
   *
   * This function returns the norm of global residual vector, which is suitable for comparison with
   * a tolerance.
   */
  virtual real64
  calculateResidualNorm( real64 const & time,
                         real64 const & dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;
  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param domain the domain containing the mesh and fields
   */
  virtual void updateState( DomainPartition & domain ) override;

  /**
   * @brief Function to determine if the solution vector should be scaled back in order to maintain a known constraint.
   * @param[in] domain The domain partition.
   * @param[in] dofManager degree-of-freedom manager associated with the linear system
   * @param[in] localSolution the solution vector
   * @return The factor that should be used to scale the solution vector values when they are being applied.
   */
  virtual real64
  scalingForSystemSolution( DomainPartition & domain,
                            DofManager const & dofManager,
                            arrayView1d< real64 const > const & localSolution )override;

  /**
   * @brief Function to check system solution for physical consistency and constraint violation
   * @param domain the domain partition
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localSolution the solution vector
   * @param scalingFactor factor to scale the solution prior to application
   * @return true if solution can be safely applied without violating physical constraints, false otherwise
   *
   * @note This function must be overridden in the derived physics solver in order to use an implict
   * solution method such as LinearImplicitStep() or NonlinearImplicitStep().
   *
   */
  virtual bool
  checkSystemSolution( DomainPartition & domain,
                       DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor ) override;

  /**
   * @brief Function to apply the solution vector to the state
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localSolution the solution vector
   * @param scalingFactor factor to scale the solution prior to application
   * @param dt the timestep
   * @param domain the domain partition
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
                       real64 const dt,
                       DomainPartition & domain ) override;

  /**
   * @brief Sets all the negative component densities (if any) to zero.
   * @param domain the physical domain object
   */
  void chopNegativeDensities( DomainPartition & domain );
#if 0
  /**
   * @brief calculate the norm of the global system residual
   * @param time the time at the beginning of the step
   * @param dt the desired timestep
   * @param domain the domain partition
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localRhs the system right-hand side vector
   * @return norm of the residual
   *
   * This function returns the norm of global residual vector, which is suitable for comparison with
   * a tolerance.
   */
  virtual real64
  calculateResidualNorm( real64 const & time,
                         real64 const & dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs );
  /**
   * @brief Function to check system solution for physical consistency and constraint violation
   * @param domain the domain partition
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localSolution the solution vector
   * @param scalingFactor factor to scale the solution prior to application
   * @return true if solution can be safely applied without violating physical constraints, false otherwise
   *
   * @note This function must be overridden in the derived physics solver in order to use an implict
   * solution method such as LinearImplicitStep() or NonlinearImplicitStep().
   *
   */
  virtual bool
  checkSystemSolution( DomainPartition & domain,
                       DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor );

  /**
   * @brief Function to determine if the solution vector should be scaled back in order to maintain a known constraint.
   * @param[in] domain The domain partition.
   * @param[in] dofManager degree-of-freedom manager associated with the linear system
   * @param[in] localSolution the solution vector
   * @return The factor that should be used to scale the solution vector values when they are being applied.
   */
  virtual real64
  scalingForSystemSolution( DomainPartition & domain,
                            DofManager const & dofManager,
                            arrayView1d< real64 const > const & localSolution );

  /**
   * @brief Function to apply the solution vector to the state
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localSolution the solution vector
   * @param scalingFactor factor to scale the solution prior to application
   * @param dt the timestep
   * @param domain the domain partition
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
                       real64 const dt,
                       DomainPartition & domain );


  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param domain the domain containing the mesh and fields
   */
  virtual void updateState( DomainPartition & domain );

  /**
   * @brief perform cleanup for implicit timestep
   * @param time the time at the beginning of the step
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
                        DomainPartition & domain ) override;

#endif

  /**
   * @brief Utility function to keep the well variables during a time step (used in
   * poromechanics simulations)
   * @param[in] keepVariablesConstantDuringInitStep flag to tell the solver to freeze its
   * primary variables during a time step
   * @detail This function is meant to be called by a specific task before/after the
   * initialization step
   */
  void setKeepVariablesConstantDuringInitStep( bool const keepVariablesConstantDuringInitStep );


protected:
  //virtual void postInputInitialization() override;

  virtual void initializePostSubGroups() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  virtual void postRestartInitialization() override final;


private:

  /// name of the flow solver
  string m_flowSolverName;

  /// flag indicating whether mass or molar formulation should be used
  integer m_useMass;

  /// flag indicating whether total mass equation should be used
  integer m_useTotalMassEquation;

  /// flag indicating whether thermal formulation is used
  integer m_isThermal;

  /// flag indicating whether compositional formulation is used
  bool m_isCompositional;


  /// number of phases
  integer m_numFluidPhases;

  /// number of components
  integer m_numFluidComponents;

  /// number of degrees of freedom per well element
  integer m_numDofPerWellElement;

  /// number of degrees of freedom per reservoir element
  integer m_numDofPerResElement;

  /// minimum value of the scaling factor obtained by enforcing maxCompFracChange
  real64 m_minScalingFactor;

  /// flag indicating whether local (cell-wise) chopping of negative compositions is allowed
  integer m_allowCompDensChopping;

  // flag to enable time step selection base on rates/bhp tables coordinates
  integer m_timeStepFromTables;
};

}

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELL_MANAGER_HPP_
