/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

/**
 * @file CompositionalMultiphaseFlow.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_COMPOSITIONALMULTIPHASEFLOW_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_COMPOSITIONALMULTIPHASEFLOW_HPP_

#include "physicsSolvers/FiniteVolume/FlowSolverBase.hpp"

class Epetra_FECrsGraph;

namespace geosx
{

namespace dataRepository
{
class ManagedGroup;
}
class BoundaryConditionBase;
class FiniteElementBase;
class DomainPartition;

/**
 * @class CompositionalMultiphaseFlow
 *
 * A compositional multiphase solver
 */
class CompositionalMultiphaseFlow : public FlowSolverBase
{
public:

  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  CompositionalMultiphaseFlow( const string& name,
                               ManagedGroup * const parent );

  /// deleted default constructor
  CompositionalMultiphaseFlow() = delete;

  /// deleted copy constructor
  CompositionalMultiphaseFlow( CompositionalMultiphaseFlow const & ) = delete;

  /// default move constructor
  CompositionalMultiphaseFlow( CompositionalMultiphaseFlow && ) = default;

  /// deleted assignment operator
  CompositionalMultiphaseFlow & operator=( CompositionalMultiphaseFlow const & ) = delete;

  /// deleted move operator
  CompositionalMultiphaseFlow & operator=( CompositionalMultiphaseFlow && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~CompositionalMultiphaseFlow() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName() { return "CompositionalMultiphaseFlow"; }

  virtual void FillDocumentationNode() override;

  virtual void FillOtherDocumentationNodes( dataRepository::ManagedGroup * const rootGroup ) override;

  virtual void InitializePreSubGroups( ManagedGroup * const rootGroup ) override;

  virtual void IntermediateInitializationPreSubGroups( ManagedGroup * const rootGroup ) override;

  virtual void FinalInitializationPreSubGroups( dataRepository::ManagedGroup * const rootGroup ) override;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual real64 SolverStep( real64 const& time_n,
                             real64 const& dt,
                             integer const cycleNumber,
                             DomainPartition * domain ) override;

  virtual void ImplicitStepSetup( real64 const& time_n,
                                  real64 const& dt,
                                  DomainPartition * const domain,
                                  systemSolverInterface::EpetraBlockSystem * const blockSystem ) override;


  virtual void AssembleSystem( DomainPartition * const domain,
                               systemSolverInterface::EpetraBlockSystem * const blockSystem,
                               real64 const time_n,
                               real64 const dt ) override;

  virtual void ApplyBoundaryConditions( DomainPartition * const domain,
                                        systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                        real64 const time_n,
                                        real64 const dt ) override;

  virtual real64
  CalculateResidualNorm(systemSolverInterface::EpetraBlockSystem const *const blockSystem,
                        DomainPartition *const domain) override;

  virtual void SolveSystem( systemSolverInterface::EpetraBlockSystem * const blockSystem,
                            SystemSolverParameters const * const params ) override;

  virtual void
  ApplySystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                       real64 const scalingFactor,
                       DomainPartition * const domain ) override;

  virtual void ResetStateToBeginningOfStep( DomainPartition * const domain ) override;

  virtual void ImplicitStepComplete( real64 const & time,
                                     real64 const & dt,
                                     DomainPartition * const domain ) override;

  /**@}*/

  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {
    using ViewKey = dataRepository::ViewKey;

    // inputs
    ViewKey temperature = { "temperature" };

    ViewKey blockLocalDofNumber    = { "blockLocalDofNumber_CompositionalMultiphaseFlow" };

    // primary solution field
    ViewKey pressure      = { "pressure" };
    ViewKey deltaPressure = { "deltaPressure" };
    ViewKey facePressure  = { "facePressure" };

    ViewKey globalCompDensity      = { "globalCompDensity" };
    ViewKey deltaGlobalCompDensity = { "deltaGlobalCompDensity" };

    // intermediate values for constitutive model input
    ViewKey globalCompMassFraction                     = { "globalCompMassFraction" };
    ViewKey dGlobalCompMassFraction_dGlobalCompDensity = { "dGlobalCompMassFraction_dGlobalCompDensity" };

    // these are used to store last converged time step values
    ViewKey phaseVolumeFraction        = { "phaseVolumeFraction" };
    ViewKey phaseDensity               = { "phaseDensity" };
    ViewKey phaseComponentMassFraction = { "phaseComponentMassFraction" };
    ViewKey phaseViscosity             = { "phaseViscosity" };
    ViewKey phaseRelativePermeability  = { "phaseRelativePermeability" };
    ViewKey porosity                   = { "porosity" };

  } viewKeysCompMultiphaseFlow;

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysCompMultiphaseFlow;


private:

  /**
   * @brief Resize the allocated multidimensional fields
   * @param domain the domain containing the mesh and fields
   *
   * Resize fields along dimensions 1 and 2 (0 is the size of containing object, i.e. element subregion)
   * once the number of phases/components is known (e.g. component fractions)
   */
  void ResizeFields( DomainPartition * domain );

  /**
   * @brief Recompute component mass fractions from primary variables (component densities)
   * @param domain the domain containing the mesh and fields
   */
  void UpdateComponentFraction( DomainPartition * domain );

  /**
   * @brief Update all relevant constitutive models using current values of pressure and composition
   * @param domain the domain containing the mesh and fields
   */
  void UpdateConstitutiveModels( DomainPartition * domain );

  /**
   * @brief Initialize all variables from initial conditions
   * @param domain the domain containing the mesh and fields
   *
   * Initialize all variables from initial conditions. This calculating primary variable values
   * from prescribed intermediate values (i.e. global densities from global mass fractions)
   * and any applicable hydrostatic equilibration of the domain
   */
  void InitializeFluidState( DomainPartition * domain );

  /**
   * @brief Backup current values of all constitutive fields that participate in the accumulation term
   * @param domain the domain containing the mesh and fields
   */

  void BackupFields( DomainPartition * domain );

  /**
   * @brief Set up the linear system (DOF indices and sparsity patterns)
   * @param domain the domain containing the mesh and fields
   * @param blockSystem the linear system object
   */
  void SetupSystem ( DomainPartition * const domain,
                     systemSolverInterface::EpetraBlockSystem * const blockSystem );

  /**
   * @brief set the sparsity pattern for the linear system
   * @param domain the domain partition
   * @param sparsity the sparsity pattern matrix
   */
  void SetSparsityPattern( DomainPartition const * const domain,
                           Epetra_FECrsGraph * const sparsity );

  /**
   * @brief sets the dof indices for this solver
   * @param meshLevel the mesh object (single level only)
   * @param numLocalRows the number of local rows on this partition
   * @param numGlobalRows the number of global rows in the problem
   * @param offset the DOF offset for this solver in the case of a non-block system
   *
   * This function sets the number of global rows, and sets the dof numbers for
   * this solver. dof numbers are referred to trilinosIndices currently.
   */
  void SetNumRowsAndTrilinosIndices( MeshLevel * const meshLevel,
                                     localIndex & numLocalRows,
                                     globalIndex & numGlobalRows,
                                     localIndex offset );

  /**
   * @brief assembles the accumulation terms for all cells
   * @param domain the physical domain object
   * @param blockSystem the entire block system
   * @param time_n previous time value
   * @param dt time step
   */
  void AssembleAccumulationTerms( DomainPartition * const domain,
                                  systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                  real64 const time_n,
                                  real64 const dt );

  /**
   * @brief assembles the flux terms for all cells
   * @param domain the physical domain object
   * @param blockSystem the entire block system
   * @param time_n previous time value
   * @param dt time step
   */
  void AssembleFluxTerms( DomainPartition * const domain,
                          systemSolverInterface::EpetraBlockSystem * const blockSystem,
                          real64 const time_n,
                          real64 const dt );

  /**
   * @brief assembles the volume balance terms for all cells
   * @param domain the physical domain object
   * @param blockSystem the entire block system
   * @param time_n previous time value
   * @param dt time step
   */
  void AssembleVolumeBalanceTerms( DomainPartition * const domain,
                                   systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                   real64 const time_n,
                                   real64 const dt );

  /**
   * @brief Function to perform the Application of Dirichlet type BC's
   * @param object the target ObjectManager for the application of the BC.
   * @param time current time
   * @param blockSystem the entire block system
   */
  void ApplyDirichletBC_implicit( DomainPartition * object,
                                  real64 const time, real64 const dt,
                                  systemSolverInterface::EpetraBlockSystem * const blockSystem );

  /**
   * @brief Function to perform the application of Dirichlet BCs on faces
   * @param domain the domain
   * @param time current time
   * @param blockSystem the entire block system
   */
  void ApplyFaceDirichletBC_implicit( DomainPartition * domain,
                                      real64 const time, real64 const dt,
                                      systemSolverInterface::EpetraBlockSystem * const blockSystem );

  /// the max number of fluid phases
  localIndex m_numPhases;

  /// the number of fluid components
  localIndex m_numComponents;

  /// the number of Degrees of Freedom per cell
  localIndex m_numDofPerCell;

  /// the (uniform) temperature
  real64 m_temperature;
};

} // namespace geosx


#endif //SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_COMPOSITIONALMULTIPHASEFLOW_HPP_
