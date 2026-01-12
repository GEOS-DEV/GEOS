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
 * @file DLFlowSolverBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FINITEVOLUME_DLFLOWSOLVERBASE_HPP_
#define GEOS_PHYSICSSOLVERS_FINITEVOLUME_DLFLOWSOLVERBASE_HPP_

#include "physicsSolvers/DLPhysicsSolverBase.hpp"
#include "common/Units.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"

// Shared memory
#include <sys/mman.h>  // Provides functions for memory management (e.g., mmap, munmap).
#include <fcntl.h>     // Provides file control options (e.g., open, O_CREAT).
#include <unistd.h>    // Provides access to the POSIX operating system API (e.g., read, write, close, fork, exec).
#include <cstring>     // Provides functions for string manipulation (e.g., memcpy, memmove, memset, strlen).
#include <semaphore.h> // Provides POSIX semaphore functionality for inter-process synchronization (e.g., sem_open, sem_wait, sem_post).

namespace geos
{

/**
 * @class DLFlowSolverBase
 *
 * Base class for finite volume fluid flow solvers.
 * Provides some common features
 */
class DLFlowSolverBase : public DLPhysicsSolverBase
{
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

public:

  /// String used to form the solverName used to register single-physics solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "flow"; }
  // Shared Memory
  struct SharedMetadataStruct
  {
    float m_t_time_n;
    float m_t_stepDt;
    float m_t_step_no;
    float m_shared_solution[2200]; // Example array to hold solution response data
    float m_shared_dofXCoords[2200];
    float m_shared_dofYCoords[2200];
    float m_shared_dofZCoords[2200];
    float m_shared_strainTrace[2200];
    float m_shared_prevSolution[2200];
  };
  // TODO: Shared memory variables and sizes to be user defined or auto calculated according to the python server needs
  // TODO: The whole logic of shared memory management to be encapsulated, with help of DLSharedMemoryManager class
  // DLSharedMemoryManager is already a member of DLPhysicsSolverBase (m_sharedMemoryManager)
  const std::string m_mem_name_meta_data = "/mem_name_meta_data";
  const std::string m_st_semaphore_python_ready = "/sem_name_python_ready";
  const std::string m_st_semaphore_cpp_ready = "/sem_name_cpp_ready";
  SharedMetadataStruct *m_shared_metadata;
  sem_t *m_sem_python_output_ready;
  sem_t *m_sem_cpp_input_ready;

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  DLFlowSolverBase( const string & name,
                    Group * const parent );


  /// deleted default constructor
  DLFlowSolverBase() = delete;

  /// deleted copy constructor
  DLFlowSolverBase( DLFlowSolverBase const & ) = delete;

  /// default move constructor
  DLFlowSolverBase( DLFlowSolverBase && ) = default;

  /// deleted assignment operator
  DLFlowSolverBase & operator=( DLFlowSolverBase const & ) = delete;

  /// deleted move operator
  DLFlowSolverBase & operator=( DLFlowSolverBase && ) = delete;

  virtual void registerDataOnMesh( Group & MeshBodies ) override;

  struct viewKeyStruct : DLPhysicsSolverBase::viewKeyStruct
  {
    // misc inputs
    static constexpr char const * isThermalString() { return "isThermal"; }
    static constexpr char const * inputTemperatureString() { return "temperature"; }
    static constexpr char const * allowNegativePressureString() { return "allowNegativePressure"; }
    static constexpr char const * maxAbsolutePresChangeString() { return "maxAbsolutePressureChange"; }
    static constexpr char const * maxSequentialPresChangeString() { return "maxSequentialPressureChange"; }
    static constexpr char const * maxSequentialTempChangeString() { return "maxSequentialTemperatureChange"; }

    static constexpr char const * fluidNamesString() { return "fluidNames"; }
    static constexpr char const * solidNamesString() { return "solidNames"; }
    static constexpr char const * permeabilityNamesString() { return "permeabilityNames"; }
    static constexpr char const * solidInternalEnergyNamesString() { return "solidInternalEnergyNames"; }
    static constexpr char const * thermalConductivityNamesString() { return "thermalConductivityNames"; }
  };

  /**
   * @brief Prepare the stencil weights by removing the contribution of the hydraulic aperture before
   * the aperture is updated
   * @param[in] domain the domain partition
   */
  void prepareStencilWeights( DomainPartition & domain ) const;

  /**
   * @brief Update the stencil weights by adding the contribution of the hydraulic aperture after
   * the aperture is updated
   * @param[in] domain the domain partition
   */
  void updateStencilWeights( DomainPartition & domain ) const;

  void enableFixedStressPoromechanicsUpdate() { m_isFixedStressPoromechanicsUpdate = true; }

  void enableJumpStabilization() { m_isJumpStabilized = true; }

  void updatePorosityAndPermeability( CellElementSubRegion & subRegion ) const;

  virtual void updatePorosityAndPermeability( SurfaceElementSubRegion & subRegion ) const;

  /**
   * @brief Utility function to save the iteration state (useful for sequential simulations)
   * @param[in] domain the domain partition
   */
  virtual void saveSequentialIterationState( DomainPartition & domain ) override;

  integer & isThermal() { return m_isThermal; }

  /**
   * @return The unit in which we evaluate the amount of fluid per element (Mass or Mole).
   */
  virtual units::Unit getMassUnit() const { return units::Unit::Mass; }

  /**
   * @brief Function to activate the flag allowing negative pressure
   */
  void allowNegativePressure() { m_allowNegativePressure = 1; }

  /**
   * @brief Utility function to keep the flow variables during a time step (used in poromechanics simulations)
   * @param[in] keepVariablesConstantDuringInitStep flag to tell the solver to freeze its primary variables during a time step
   * @detail This function is meant to be called by a specific task before/after the initialization step
   */
  void setKeepVariablesConstantDuringInitStep( bool const keepVariablesConstantDuringInitStep )
  { m_keepVariablesConstantDuringInitStep = keepVariablesConstantDuringInitStep; }

  virtual bool checkSequentialSolutionIncrements( DomainPartition & domain ) const override;

  void enableLaggingFractureStencilWeightsUpdate(){ m_isLaggingFractureStencilWeightsUpdate = 1; };

  real64 sumAquiferFluxes( BoundaryStencil const & stencil,
                           AquiferBoundaryCondition::KernelWrapper const & aquiferBCWrapper,
                           ElementViewConst< arrayView1d< real64 const > > const & pres,
                           ElementViewConst< arrayView1d< real64 const > > const & presOld,
                           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                           real64 const & timeAtBeginningOfStep,
                           real64 const & dt );

  /**
   * @brief assembles the flux terms for all cells for the hydrofracture case
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   * @param dR_dAper
   */
  virtual void assembleHydrofracFluxTerms( real64 const time_n,
                                           real64 const dt,
                                           DomainPartition const & domain,
                                           DofManager const & dofManager,
                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                           arrayView1d< real64 > const & localRhs,
                                           CRSMatrixView< real64, localIndex const > const & dR_dAper )
  {
    GEOS_UNUSED_VAR ( time_n, dt, domain, dofManager, localMatrix, localRhs, dR_dAper );
    GEOS_ERROR( "Poroelastic fluxes with conforming fractures not yet implemented." );
  }

  void initializeState( DomainPartition & domain );

  virtual void initializeFluidState( MeshLevel & mesh, string_array const & regionNames ) { GEOS_UNUSED_VAR( mesh, regionNames ); }

  virtual void initializeThermalState( MeshLevel & mesh, string_array const & regionNames ) { GEOS_UNUSED_VAR( mesh, regionNames ); }

  /**
   * @brief For each equilibrium initial condition, loop over all the target cells and compute the min/max elevation
   * @param[in] domain the domain partition
   * @param[in] equilNameToEquilId the map from the name of the initial condition to the initial condition index (used in min/maxElevation)
   * @param[out] maxElevation the max elevation for each initial condition
   * @param[out] minElevation the min elevation for each initial condition
   */
  void findMinMaxElevationInEquilibriumTarget( DomainPartition & domain, // cannot be const...
                                               stdMap< string, localIndex > const & equilNameToEquilId,
                                               arrayView1d< real64 > const & maxElevation,
                                               arrayView1d< real64 > const & minElevation ) const;

  /**
   * @brief For each source flux boundary condition, loop over all the target cells and sum the owned cells
   * @param[in] time the time at the beginning of the time step
   * @param[in] dt the time step size
   * @param[in] domain the domain partition
   * @param[in] bcNameToBcId the map from the name of the boundary condition to the boundary condition index
   * @param[out] bcAllSetsSize the total number of owned cells for each source flux boundary condition
   */
  void computeSourceFluxSizeScalingFactor( real64 const & time,
                                           real64 const & dt,
                                           DomainPartition & domain, // cannot be const...
                                           stdMap< string, localIndex > const & bcNameToBcId,
                                           arrayView1d< globalIndex > const & bcAllSetsSize ) const;

  integer numberOfDofsPerCell() const { return m_numDofPerCell; }

  virtual void initializeSharedMemories() override;

  virtual void shareDLModelInputs( real64 const & time_n,
                                   real64 const & dt,
                                   integer const cycleNumber,
                                   DomainPartition & domain ) override;

  virtual void readDLModelOutputs( real64 const & time_n,
                                   real64 const & dt,
                                   integer const cycleNumber,
                                   DomainPartition & domain ) override;

protected:

  /**
   * @brief Increment the cumulative flux from each aquifer
   * @param[in] time the time at the beginning of the time step
   * @param[in] dt the time step size
   * @param[in] domain the domain partition
   *
   * For now this function is here because it can be used for both single-phase flow and multiphase flow
   * This may have to be revisited when aquifer BC is implemented for hybrid FVM
   */
  virtual void saveAquiferConvergedState( real64 const & time,
                                          real64 const & dt,
                                          DomainPartition & domain );

  /**
   * @brief Utility function to save the converged state
   * @param[in] subRegion the element subRegion
   */
  virtual void saveConvergedState( ElementSubRegionBase & subRegion ) const;

  /**
   * @brief Helper function to compute/report the elements with small pore volumes
   * @param[in] domain the domain partition
   */
  virtual void validatePoreVolumes( DomainPartition const & domain ) const;

  virtual void precomputeData( MeshLevel & mesh,
                               string_array const & regionNames );

  virtual void initializePreSubGroups() override;

  /**
   * @brief Checks the validity of the discretization name for the FiniteVolume method (errors if issues are detected)
   */
  void checkDiscretizationName() const;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  virtual void computeHydrostaticEquilibrium( DomainPartition & domain ) { GEOS_UNUSED_VAR( domain ); }

  void initializePorosityAndPermeability( MeshLevel & mesh, string_array const & regionNames );

  void initializeHydraulicAperture( MeshLevel & mesh, string_array const & regionNames );

  void saveInitialPressureAndTemperature( MeshLevel & mesh, string_array const & regionNames );

  virtual void setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const override;

  /// the number of Degrees of Freedom per cell
  integer m_numDofPerCell;

  /// flag to determine whether or not this is a thermal simulation
  integer m_isThermal;

  /// the input temperature
  real64 m_inputTemperature;

  /// flag to freeze the initial state during initialization in coupled problems
  bool m_keepVariablesConstantDuringInitStep;

  /// enable the fixed stress poromechanics update of porosity
  bool m_isFixedStressPoromechanicsUpdate;

  /// enable pressure jump stabilzation for fixed-stress poromechanics
  bool m_isJumpStabilized;

  /// flag if negative pressure is allowed
  integer m_allowNegativePressure;

  /// maximum (absolute) pressure change in a Newton iteration
  real64 m_maxAbsolutePresChange;

  /// maximum (absolute) pressure change in a sequential iteration
  real64 m_sequentialPresChange;
  real64 m_maxSequentialPresChange;

  /// maximum (absolute) temperature change in a sequential iteration
  real64 m_sequentialTempChange;
  real64 m_maxSequentialTempChange;

  /**
   * @brief Class used for displaying boundary warning message
   */
  class BCMessage
  {
public:
    static string pressureConflict( string_view regionName, string_view subRegionName,
                                    string_view setName, string_view fieldName );

    static string temperatureConflict( string_view regionName, string_view subRegionName,
                                       string_view setName, string_view fieldName );

    static string missingPressure( string_view regionName, string_view subRegionName,
                                   string_view setName, string_view fieldName );

    static string missingTemperature( string_view regionName, string_view subRegionName,
                                      string_view setName, string_view fieldName );

    static string conflictingComposition( int comp, string_view componentName,
                                          string_view regionName, string_view subRegionName,
                                          string_view setName, string_view fieldName );

    static string invalidComponentIndex( int comp,
                                         string_view fsName, string_view fieldName );

    static string notAppliedOnRegion( int componentIndex, string_view componentName,
                                      string_view regionName, string_view subRegionName,
                                      string_view setName, string_view fieldName );
private:
    static string generateMessage( string_view baseMessage,
                                   string_view fieldName, string_view setName );
  };

private:
  virtual void setConstitutiveNames( ElementSubRegionBase & subRegion ) const override;

  // flag to determine whether or not to apply lagging update for the fracture stencil weights
  integer m_isLaggingFractureStencilWeightsUpdate;

};


}

#endif //GEOS_PHYSICSSOLVERS_FINITEVOLUME_DLFLOWSOLVERBASE_HPP_
