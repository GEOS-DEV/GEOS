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

/**
 * @file CompositionalMultiphaseWell.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_WELLS_COMPOSITIONALMULTIPHASEWELLSOLVER_HPP_
#define SRC_COMPONENTS_CORE_SRC_WELLS_COMPOSITIONALMULTIPHASEWELLSOLVER_HPP_

#include <constitutive/RelPerm/RelativePermeabilityBase.hpp>
#include "WellSolverBase.hpp"

class Epetra_FECrsGraph;

namespace geosx
{

namespace dataRepository
{
class ManagedGroup;

namespace keys
{
string const compositionalMultiphaseWell = "CompositionalMultiphaseWell";
}
}

namespace constitutive
{
class MultiFluidBase;
}

/**
 * @class CompositionalMultiphaseWell
 *
 * A compositional multiphase well solver
 */
class CompositionalMultiphaseWell : public WellSolverBase
{
public:

  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  CompositionalMultiphaseWell( const string& name,
                                     ManagedGroup * const parent );

  /// deleted default constructor
  CompositionalMultiphaseWell() = delete;

  /// deleted copy constructor
  CompositionalMultiphaseWell( CompositionalMultiphaseWell const & ) = delete;

  /// default move constructor
  CompositionalMultiphaseWell( CompositionalMultiphaseWell && ) = default;

  /// deleted assignment operator
  CompositionalMultiphaseWell & operator=( CompositionalMultiphaseWell const & ) = delete;

  /// deleted move operator
  CompositionalMultiphaseWell & operator=( CompositionalMultiphaseWell && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~CompositionalMultiphaseWell() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName() { return dataRepository::keys::compositionalMultiphaseWell; }

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
   * @brief assembles the accumulation terms for all segments
   * @param domain the physical domain object
   * @param blockSystem the entire block system
   * @param time_n previous time value
   * @param dt time step
   */
  void AssembleAccumulationTerms( DomainPartition * const domain,
                                  Epetra_FECrsMatrix * const jacobian,
                                  Epetra_FEVector * const residual,
                                  real64 const time_n,
                                  real64 const dt );

  /**
   * @brief assembles the flux terms for all connections
   * @param domain the physical domain object
   * @param blockSystem the entire block system
   * @param time_n previous time value
   * @param dt time step
   */
  void AssembleFluxTerms( DomainPartition * const domain,
                          Epetra_FECrsMatrix * const jacobian,
                          Epetra_FEVector * const residual,
                          real64 const time_n,
                          real64 const dt );

  /**
   * @brief assembles the volume balance terms for all perforations
   * @param domain the physical domain object
   * @param blockSystem the entire block system
   * @param time_n previous time value
   * @param dt time step
   */
  void AssembleVolumeBalanceTerms( DomainPartition * const domain,
                                   Epetra_FECrsMatrix * const jacobian,
                                   Epetra_FEVector * const residual,
                                   real64 const time_n,
                                   real64 const dt );


  /**
   * @brief assembles the perforation rate terms 
   * @param domain the physical domain object
   * @param blockSystem the entire block system
   * @param time_n previous time value
   * @param dt time step
   */

  void AssembleSourceTerms( DomainPartition * const domain,
                            systemSolverInterface::EpetraBlockSystem * const blockSystem,
                            real64 const time_n,
                            real64 const dt );

  void CheckWellControlSwitch( DomainPartition * const domain );
  
  /**@}*/

  struct viewKeyStruct : WellSolverBase::viewKeyStruct
  {
    // inputs
    static constexpr auto temperatureString = "temperature";
    static constexpr auto useMassFlagString = "useMass";

    static constexpr auto relPermNameString  = "relPermName";
    static constexpr auto relPermIndexString = "relPermIndex";
    
    // primary solution field
    static constexpr auto pressureString               = "pressure";
    static constexpr auto deltaPressureString          = "deltaPressure";
    static constexpr auto globalCompDensityString      = "globalCompDensity";
    static constexpr auto deltaGlobalCompDensityString = "deltaGlobalCompDensity";
    static constexpr auto mixtureVelocityString        = "mixtureVelocity";
    static constexpr auto deltaMixtureVelocityString   = "deltaMixtureVelocity";

    static constexpr auto phaseFlowRateString = "phaseFlowRate";
    static constexpr auto bhpString           = "bhp";

    // intermediate values for constitutive model input
    static constexpr auto globalComponentFracString        = "globalComponentFraction";
    static constexpr auto dGlobalComponentFrac_dPresString = "dGlobalComponentFraction_dPres";
    static constexpr auto dGlobalComponentFrac_dCompString = "dGlobalComponentFraction_dComp";

    using ViewKey = dataRepository::ViewKey;

    // inputs
    ViewKey temperature = { temperatureString };
    ViewKey useMassFlag = { useMassFlagString };

    ViewKey relPermName  = { relPermNameString };
    ViewKey relPermIndex = { relPermIndexString };
    
    // primary solution field
    ViewKey pressure               = { pressureString };
    ViewKey deltaPressure          = { deltaPressureString };
    ViewKey globalCompDensity      = { globalCompDensityString };
    ViewKey deltaGlobalCompDensity = { deltaGlobalCompDensityString };
    ViewKey mixtureVelocity        = { mixtureVelocityString };
    ViewKey deltaMixtureVelovity   = { deltaMixtureVelocityString };
    
    // well controls
    ViewKey phaseFlowRate = { phaseFlowRateString };
    ViewKey bhp           = { bhpString };

    // global composition to input injection stream
    ViewKey globalComponentFrac        = { globalComponentFracString };
    ViewKey dGlobalComponentFrac_dPres = { dGlobalComponentFrac_dPresString };
    ViewKey dGlobalComponentFrac_dComp = { dGlobalComponentFrac_dCompString };    

  } viewKeysCompMultiphaseWell;

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysCompMultiphaseWell;

protected:
  virtual void InitializePreSubGroups( ManagedGroup * const rootGroup ) override;

  virtual void InitializePostInitialConditions_PreSubGroups( dataRepository::ManagedGroup * const rootGroup ) override;


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
   * @brief Initialize all well variables from initial conditions and injection stream
   * @param domain the domain containing the mesh and fields
   */

  void InitializeWellState( DomainPartition * const domain );
  
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


#endif //SRC_COMPONENTS_CORE_SRC_WELLS_COMPOSITIONALMULTIPHASEWELLS_HPP_
