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

#include <constitutive/RelPerm/RelativePermeabilityBase.hpp>
#include "physicsSolvers/FiniteVolume/FlowSolverBase.hpp"

class Epetra_FECrsGraph;

namespace geosx
{

namespace dataRepository
{
class ManagedGroup;

namespace keys
{
string const compositionalMultiphaseFlow = "CompositionalMultiphaseFlow";
}
}
class BoundaryConditionBase;
class FiniteElementBase;
class DomainPartition;

namespace constitutive
{
class MultiFluidBase;
}

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
  static string CatalogName() { return dataRepository::keys::compositionalMultiphaseFlow; }

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

  /**
   * @brief Recompute component fractions from primary variables (component densities)
   * @param domain the domain containing the mesh and fields
   */
  void UpdateComponentFraction( ManagedGroup * dataGroup );

  /**
   * @brief Recompute component fractions from primary variables (component densities)
   * @param domain the domain containing the mesh and fields
   */
  void UpdateComponentFractionAll( DomainPartition * domain );

  /**
   * @brief Recompute phase volume fractions (saturations) from constitutive and primary variables
   * @param domain the domain containing the mesh and fields
   */
  void UpdatePhaseVolumeFraction( ManagedGroup * dataGroup );

  /**
   * @brief Recompute phase volume fractions (saturations) from constitutive and primary variables
   * @param domain the domain containing the mesh and fields
   */
  void UpdatePhaseVolumeFractionAll( DomainPartition * domain );

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param dataGroup the group storing the required fields
   */
  void UpdateFluidModel( ManagedGroup * dataGroup );

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param domain the domain containing the mesh and fields
   */
  void UpdateFluidModelAll( DomainPartition * domain );

  /**
   * @brief Update all relevant solid models using current values of pressure
   * @param dataGroup the group storing the required fields
   */
  void UpdateSolidModel( ManagedGroup * dataGroup );

  /**
   * @brief Update all relevant solid models using current values of pressure
   * @param domain the domain containing the mesh and fields
   */
  void UpdateSolidModelAll( DomainPartition * domain );

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param dataGroup the group storing the required fields
   */
  void UpdateRelPermModel( ManagedGroup * dataGroup );

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param domain the domain containing the mesh and fields
   */
  void UpdateRelPermModelAll( DomainPartition * domain );

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param domain the domain containing the mesh and fields
   */
  void UpdateState( ManagedGroup * dataGroup );

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param domain the domain containing the mesh and fields
   */
  void UpdateStateAll( DomainPartition * domain );

  /**
   * @brief Get the number of fluid components (species)
   * @return the number of components
   */
  localIndex numFluidComponents() const;

  /**
   * @brief Get the number of fluid phases
   * @return the number of phases
   */
  localIndex numFluidPhases() const;

  /**
   * @brief assembles the accumulation terms for all cells
   * @param domain the physical domain object
   * @param blockSystem the entire block system
   * @param time_n previous time value
   * @param dt time step
   */
  void AssembleAccumulationTerms( DomainPartition const * const domain,
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
  void AssembleFluxTerms( DomainPartition const * const domain,
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
  void AssembleVolumeBalanceTerms( DomainPartition const * const domain,
                                   systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                   real64 const time_n,
                                   real64 const dt );

  /**@}*/

  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {
    // inputs
    static constexpr auto temperatureString = "temperature";
    static constexpr auto useMassFlagString = "useMass";

    static constexpr auto relPermNameString  = "relPermName";
    static constexpr auto relPermIndexString = "relPermIndex";

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

    ViewKey relPermName  = { relPermNameString };
    ViewKey relPermIndex = { relPermIndexString };

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
   * @brief Extract the fluid model used by this solver from a group
   * @param dataGroup target group (e.g. subregion, face/edge/node manager, etc.)
   * @return
   */
  constitutive::MultiFluidBase * GetFluidModel( ManagedGroup * const dataGroup ) const;

  /**
   * @brief Extract the fluid model used by this solver from a group (const version)
   * @param dataGroup target group (e.g. subregion, face/edge/node manager, etc.)
   * @return
   */
  constitutive::MultiFluidBase const * GetFluidModel( ManagedGroup const * const dataGroup ) const;

  /**
   * @brief Extract the solid model used by this solver from a group
   * @param dataGroup target group (e.g. subregion, face/edge/node manager, etc.)
   * @return
   */
  constitutive::ConstitutiveBase * GetSolidModel( ManagedGroup * const dataGroup ) const;

  /**
   * @brief Extract the solid model used by this solver from a group (const version)
   * @param dataGroup target group (e.g. subregion, face/edge/node manager, etc.)
   * @return
   */
  constitutive::ConstitutiveBase const * GetSolidModel( ManagedGroup const * const dataGroup ) const;

  /**
   * @brief Extract the relative permeability model used by this solver from a group
   * @param dataGroup target group (e.g. subregion, face/edge/node manager, etc.)
   * @return
   */
  constitutive::RelativePermeabilityBase * GetRelPermModel( ManagedGroup * const dataGroup ) const;

  /**
   * @brief Extract the relative permeability model used by this solver from a group (const version)
   * @param dataGroup target group (e.g. subregion, face/edge/node manager, etc.)
   * @return
   */
  constitutive::RelativePermeabilityBase const * GetRelPermModel( ManagedGroup const * const dataGroup ) const;

  /**
   * @brief Resize the allocated multidimensional fields
   * @param domain the domain containing the mesh and fields
   *
   * Resize fields along dimensions 1 and 2 (0 is the size of containing object, i.e. element subregion)
   * once the number of phases/components is known (e.g. component fractions)
   */
  void ResizeFields( DomainPartition * const domain );

  /**
   * @brief Initialize all variables from initial conditions
   * @param domain the domain containing the mesh and fields
   *
   * Initialize all variables from initial conditions. This calculating primary variable values
   * from prescribed intermediate values (i.e. global densities from global fractions)
   * and any applicable hydrostatic equilibration of the domain
   */
  void InitializeFluidState( DomainPartition * const domain );

  /**
   * @brief Backup current values of all constitutive fields that participate in the accumulation term
   * @param domain the domain containing the mesh and fields
   */
  void BackupFields( DomainPartition * const domain );

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
   * @brief Function to perform the Application of Dirichlet type BC's
   * @param object the target ObjectManager for the application of the BC.
   * @param time current time
   * @param blockSystem the entire block system
   */
  void ApplyDirichletBC_implicit( DomainPartition * const object,
                                  real64 const time, real64 const dt,
                                  systemSolverInterface::EpetraBlockSystem * const blockSystem );

  /**
   * @brief Function to perform the application of Dirichlet BCs on faces
   * @param domain the domain
   * @param time current time
   * @param blockSystem the entire block system
   */
  void ApplyFaceDirichletBC_implicit( DomainPartition * const domain,
                                      real64 const time, real64 const dt,
                                      systemSolverInterface::EpetraBlockSystem * const blockSystem );

  /// the max number of fluid phases
  localIndex m_numPhases;

  /// the number of fluid components
  localIndex m_numComponents;

  /// the (uniform) temperature
  real64 m_temperature;

  /// flag indicating whether mass or molar formulation should be used
  integer m_useMass;

  /// name of the rel perm constitutive model
  string m_relPermName;

  /// index of the rel perm constitutive model
  localIndex m_relPermIndex;
};

} // namespace geosx


#endif //SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_COMPOSITIONALMULTIPHASEFLOW_HPP_
