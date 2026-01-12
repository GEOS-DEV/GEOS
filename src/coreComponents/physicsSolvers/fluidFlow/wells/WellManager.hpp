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
#include "physicsSolvers/fluidFlow/wells/WellSolverBase.hpp"

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
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName() { return "WellManager"; }
  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  /**
   * @brief Get a well solver for a given well element sub-region
   * @param subRegion the well subRegion whose well solver is requested
   * @return a reference to the well solver
   */
  WellSolverBase & getWell( WellElementSubRegion const & subRegion );

  /* PhysicsSolverBase interfaces */
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

#if 0
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
                  arrayView1d< real64 > const & localRhs );

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
};

}

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELL_MANAGER_HPP_
