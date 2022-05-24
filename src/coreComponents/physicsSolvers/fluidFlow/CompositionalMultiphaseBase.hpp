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

/**
 * @file CompositionalMultiphaseBase.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASE_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASE_HPP_

#include "common/DataLayouts.hpp"
#include "constitutive/fluid/layouts.hpp"
#include "constitutive/relativePermeability/layouts.hpp"
#include "constitutive/capillaryPressure/layouts.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"

namespace geosx
{

//START_SPHINX_INCLUDE_00
/**
 * @class CompositionalMultiphaseBase
 *
 * A compositional multiphase solver
 */
class CompositionalMultiphaseBase : public FlowSolverBase
{
public:

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  CompositionalMultiphaseBase( const string & name,
                               Group * const parent );

  /// deleted default constructor
  CompositionalMultiphaseBase() = delete;

  /// deleted copy constructor
  CompositionalMultiphaseBase( CompositionalMultiphaseBase const & ) = delete;

  /// default move constructor
  CompositionalMultiphaseBase( CompositionalMultiphaseBase && ) = default;

  /// deleted assignment operator
  CompositionalMultiphaseBase & operator=( CompositionalMultiphaseBase const & ) = delete;

  /// deleted move operator
  CompositionalMultiphaseBase & operator=( CompositionalMultiphaseBase && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~CompositionalMultiphaseBase() override = default;

//START_SPHINX_INCLUDE_01

  virtual void registerDataOnMesh( Group & meshBodies ) override;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual real64
  solverStep( real64 const & time_n,
              real64 const & dt,
              integer const cycleNumber,
              DomainPartition & domain ) override;

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  assembleSystem( real64 const time_n,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual void
  applyBoundaryConditions( real64 const time_n,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  virtual void
  solveLinearSystem( DofManager const & dofManager,
                     ParallelMatrix & matrix,
                     ParallelVector & rhs,
                     ParallelVector & solution ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  /**
   * @brief Recompute component fractions from primary variables (component densities)
   * @param dataGroup the group storing the required fields
   */
  void updateComponentFraction( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Recompute phase volume fractions (saturations) from constitutive and primary variables
   * @param dataGroup the group storing the required fields
   */
  void updatePhaseVolumeFraction( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param dataGroup the group storing the required fields
   */
  void updateFluidModel( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Update all relevant relperm models using current values of phase volume fraction
   * @param dataGroup the group storing the required fields
   */
  void updateRelPermModel( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Update all relevant capillary pressure models using current values of phase volume fraction
   * @param dataGroup the group storing the required fields
   */
  void updateCapPressureModel( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Update all relevant solid internal energy models using current values of temperature
   * @param dataGroup the group storing the required fields
   */
  void updateSolidInternalEnergyModel( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Recompute phase mobility from constitutive and primary variables
   * @param dataGroup the group storing the required field
   */
  virtual void updatePhaseMobility( ObjectManagerBase & dataGroup ) const = 0;

  void updateFluidState( ObjectManagerBase & dataGroup ) const;

  virtual void updateState( DomainPartition & domain ) override final;

  /**
   * @brief Getter for the number of fluid components (species)
   * @return the number of components
   */
  integer numFluidComponents() const { return m_numComponents; }

  /**
   * @brief Getter for the number of fluid phases
   * @return the number of phases
   */
  integer numFluidPhases() const { return m_numPhases; }

  /**
   * @brief Getter for the name of the reference fluid model name
   * @return the name of the reference fluid
   */
  string referenceFluidModelName() const { return m_referenceFluidModelName; }

  /**
   * @brief assembles the accumulation and volume balance terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   */
  void assembleAccumulationAndVolumeBalanceTerms( DomainPartition & domain,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief assembles the flux terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void
  assembleFluxTerms( real64 const dt,
                     DomainPartition const & domain,
                     DofManager const & dofManager,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs ) const = 0;


  /**@}*/

  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {
    static constexpr char const * elemDofFieldString() { return "compositionalVariables"; }

    // inputs

    static constexpr char const * inputTemperatureString() { return "temperature"; }
    static constexpr char const * useMassFlagString() { return "useMass"; }
    static constexpr char const * isThermalString()  { return "isThermal"; }
    static constexpr char const * relPermNamesString() { return "relPermNames"; }
    static constexpr char const * capPressureNamesString() { return "capPressureNames"; }
    static constexpr char const * thermalConductivityNamesString() { return "thermalConductivityNames"; }
    static constexpr char const * solidInternalEnergyNamesString() { return "solidInternalEnergyNames"; }
    static constexpr char const * maxCompFracChangeString() { return "maxCompFractionChange"; }
    static constexpr char const * allowLocalCompDensChoppingString() { return "allowLocalCompDensityChopping"; }

    // reservoir region statistics

    static constexpr char const * phasePoreVolumeString() { return "phasePoreVolume"; }

  };

  /**
   * @brief Initialize all variables from initial conditions
   * @param domain the domain containing the mesh and fields
   *
   * Initialize all variables from initial conditions. This calculating primary variable values
   * from prescribed intermediate values (i.e. global densities from global fractions)
   * and any applicable hydrostatic equilibration of the domain
   */
  void initializeFluidState( MeshLevel & mesh, arrayView1d< string const > const & regionNames );

  /**
   * @brief Compute the hydrostatic equilibrium using the compositions and temperature input tables
   */
  void computeHydrostaticEquilibrium();

  /**
   * @brief Backup current values of all constitutive fields that participate in the accumulation term
   * @param domain the domain containing the mesh and fields
   */
  void backupFields( MeshLevel & mesh, arrayView1d< string const > const & regionNames ) const;

  /**
   * @brief Function to perform the Application of Dirichlet type BC's
   * @param time current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param localMatrix local system matrix
   * @param localRhs local system right-hand side vector
   */
  void applyDirichletBC( real64 const time,
                         real64 const dt,
                         DofManager const & dofManager,
                         DomainPartition & domain,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief Apply source flux boundary conditions to the system
   * @param time current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param localMatrix local system matrix
   * @param localRhs local system right-hand side vector
   */
  void applySourceFluxBC( real64 const time,
                          real64 const dt,
                          DofManager const & dofManager,
                          DomainPartition & domain,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief Apply aquifer boundary conditions to the system
   * @param time current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param localMatrix local system matrix
   * @param localRhs local system right-hand side vector
   */
  virtual void applyAquiferBC( real64 const time,
                               real64 const dt,
                               DofManager const & dofManager,
                               DomainPartition & domain,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) const = 0;


  /**
   * @brief Sets all the negative component densities (if any) to zero.
   * @param domain the physical domain object
   */
  void chopNegativeDensities( DomainPartition & domain );

  /**
   * @brief Utility function to freeze the flow variables during a time step
   * @param[in] freezeFlowVariablesDuringStep flag to tell the solver to freeze its primary variables during a time step
   * @detail This function is meant to be called by a specific task before/after the initialization step
   */
  void freezeFlowVariablesDuringStep( bool freezeFlowVariablesDuringStep )
  { m_freezeFlowVariablesDuringStep = freezeFlowVariablesDuringStep; }

  virtual void initializePostInitialConditionsPreSubGroups() override;

protected:

  virtual void postProcessInput() override;

  virtual void initializePreSubGroups() override;

  virtual void computeRegionStatistics( MeshLevel & mesh,
                                        arrayView1d< string const > const & regionNames ) const override;

  /**
   * @brief Function to fix the initial state during the initialization step in coupled problems
   * @param time current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param localMatrix local system matrix
   * @param localRhs local system right-hand side vector
   *
   * @detail This function is meant to be called when the flag m_freezeFlowVariablesDuringStep is on
   *         The main use case is the initialization step in coupled problems during which we solve an elastic problem for a fixed pressure
   */
  void freezeFlowVariablesDuringStep( real64 const time,
                                      real64 const dt,
                                      DofManager const & dofManager,
                                      DomainPartition & domain,
                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                      arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief Utility function that checks the consistency of the constitutive models
   * @param[in] domain the domain partition
   * This function will produce an error if one of the constitutive models
   * (fluid, relperm) is incompatible with the reference fluid model.
   */
  void validateConstitutiveModels( DomainPartition const & domain ) const;

  /**
   * @brief Initialize the aquifer boundary condition (gravity vector, water phase index)
   * @param[in] cm reference to the global constitutive model manager
   */
  void initializeAquiferBC( constitutive::ConstitutiveManager const & cm ) const;


  /// flag to specify whether the sparsity pattern needs to be rebuilt
  bool m_systemSetupDone;

  /// the max number of fluid phases
  integer m_numPhases;

  /// the number of fluid components
  integer m_numComponents;

  /// the input temperature
  real64 m_inputTemperature;

  /// flag indicating whether mass or molar formulation should be used
  integer m_useMass;

  /// flag to determine whether or not to apply capillary pressure
  integer m_hasCapPressure;

  /// flag to determine whether or not this is a thermal simulation
  integer m_isThermal;

  /// flag to fix the initial state during initialization in coupled problems
  integer m_freezeFlowVariablesDuringStep;

  /// maximum (absolute) change in a component fraction between two Newton iterations
  real64 m_maxCompFracChange;

  /// minimum value of the scaling factor obtained by enforcing maxCompFracChange
  real64 m_minScalingFactor;

  /// flag indicating whether local (cell-wise) chopping of negative compositions is allowed
  integer m_allowCompDensChopping;

  /// name of the fluid constitutive model used as a reference for component/phase description
  string m_referenceFluidModelName;

private:

  virtual void setConstitutiveNames( ElementSubRegionBase & subRegion ) const override;

};

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASE_HPP_
