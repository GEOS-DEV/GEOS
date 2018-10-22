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

  virtual bool
  CheckSystemSolution(systemSolverInterface::EpetraBlockSystem const * const blockSystem, real64 const scalingFactor,
                      DomainPartition * const domain) override;

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
    // inputs
    static constexpr auto temperatureString = "temperature";
    static constexpr auto useMassFlagString = "useMass";

    static constexpr auto blockLocalDofNumberString    = "blockLocalDofNumber_CompositionalMultiphaseFlow";

    // primary solution field
    static constexpr auto pressureString      = "pressure";
    static constexpr auto deltaPressureString = "deltaPressure";
    static constexpr auto facePressureString  = "facePressure";
    static constexpr auto bcPressureString    = "bcPressure";

    static constexpr auto globalCompDensityString      = "globalCompDensity";
    static constexpr auto deltaGlobalCompDensityString = "deltaGlobalCompDensity";

    // intermediate values for constitutive model input
    static constexpr auto globalCompFractionString                     = "globalCompFraction";
    static constexpr auto dGlobalCompFraction_dGlobalCompDensityString = "dGlobalCompFraction_dGlobalCompDensity";

    // intermediate values for saturations
    static constexpr auto phaseVolumeFractionString                     = "phaseVolumeFraction";
    static constexpr auto dPhaseVolumeFraction_dPressureString          = "dPhaseVolumeFraction_dPressure";
    static constexpr auto dPhaseVolumeFraction_dGlobalCompDensityString = "dPhaseVolumeFraction_dGlobalCompDensity";

    // these are used to store last converged time step values
    static constexpr auto phaseVolumeFractionOldString     = "phaseVolumeFractionOld";
    static constexpr auto phaseDensityOldString            = "phaseDensityOld";
    static constexpr auto phaseComponentFractionOldString  = "phaseComponentFractionOld";
    static constexpr auto porosityOldString                = "porosityOld";

    // these are allocated on faces for BC application until we can get constitutive models on faces
    static constexpr auto phaseViscosityString             = "phaseViscosity";
    static constexpr auto phaseRelativePermeabilityString  = "phaseRelativePermeability";

    using ViewKey = dataRepository::ViewKey;

    // inputs
    ViewKey temperature = { temperatureString };
    ViewKey useMassFlag = { useMassFlagString };

    ViewKey blockLocalDofNumber    = { blockLocalDofNumberString };

    // primary solution field
    ViewKey pressure      = { pressureString };
    ViewKey deltaPressure = { deltaPressureString };
    ViewKey facePressure  = { facePressureString };
    ViewKey bcPressure    = { bcPressureString };

    ViewKey globalCompDensity      = { globalCompDensityString };
    ViewKey deltaGlobalCompDensity = { deltaGlobalCompDensityString };

    // intermediate values for constitutive model input
    ViewKey globalCompFraction                     = { globalCompFractionString };
    ViewKey dGlobalCompFraction_dGlobalCompDensity = { dGlobalCompFraction_dGlobalCompDensityString };

    // intermediate values for saturations
    ViewKey phaseVolumeFraction                     = { phaseVolumeFractionString };
    ViewKey dPhaseVolumeFraction_dPressure          = { dPhaseVolumeFraction_dPressureString };
    ViewKey dPhaseVolumeFraction_dGlobalCompDensity = { dPhaseVolumeFraction_dGlobalCompDensityString };

    // these are used to store last converged time step values
    ViewKey phaseVolumeFractionOld     = { phaseVolumeFractionOldString };
    ViewKey phaseDensityOld            = { phaseDensityOldString };
    ViewKey phaseComponentFractionOld  = { phaseComponentFractionOldString };
    ViewKey porosityOld                = { porosityOldString };

    // these are allocated on faces for BC application until we can get constitutive models on faces
    ViewKey phaseViscosity             = { phaseViscosityString };
    ViewKey phaseRelativePermeability  = { phaseRelativePermeabilityString };

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
   * @brief Recompute component fractions from primary variables (component densities)
   * @param domain the domain containing the mesh and fields
   */
  void UpdateComponentFraction( DomainPartition * domain );

  /**
   * @brief Recompute phase volume fractions (saturations) from constitutive and primary variables
   * @param domain the domain containing the mesh and fields
   */
  void UpdatePhaseVolumeFraction( DomainPartition * domain );

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param dataGroup the group storing the required fields
   * @param targetSet the set to call model updates on
   */
  void UpdateFluidModel( ManagedGroup * dataGroup,
                         set<localIndex> const & targetSet,
                         string const & pressureFieldName = viewKeyStruct::pressureString,
                         string const & deltaPressureFieldName = viewKeyStruct::deltaPressureString,
                         string const & compFracFieldName = viewKeyStruct::globalCompFractionString );

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param dataGroup the group storing the required fields
   * @param setName name of the set to call model updates on
   */
  void UpdateFluidModel( ManagedGroup * dataGroup, string const & setName );

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param domain the domain containing the mesh and fields
   */
  void UpdateFluidModels( DomainPartition * domain );

  /**
 * @brief Update all relevant solid models using current values of pressure
 * @param dataGroup the group storing the required fields
 * @param targetSet the set to call model updates on
 */
  void UpdateSolidModel( ManagedGroup * dataGroup, set<localIndex> const & targetSet );

  /**
   * @brief Update all relevant solid models using current values of pressure
   * @param dataGroup the group storing the required fields
   * @param setName name of the set to call model updates on
   */
  void UpdateSolidModel( ManagedGroup * dataGroup, string const & setName );

  /**
   * @brief Update all relevant solid models using current values of pressure
   * @param domain the domain containing the mesh and fields
   */
  void UpdateSolidModels( DomainPartition * domain );

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
   * from prescribed intermediate values (i.e. global densities from global fractions)
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

  /// flag indicating whether mass or molar formulation should be used
  integer m_useMass;
};

} // namespace geosx


#endif //SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_COMPOSITIONALMULTIPHASEFLOW_HPP_
