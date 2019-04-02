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
class CompositionalMultiphaseFlow;
class Well;
  
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

  struct ColOffset
  {
    static constexpr integer DPRES = 0;
    static constexpr integer DRATE = 1;
  };

  struct RowOffset
  {
    static constexpr integer CONTROL = 0;
    static constexpr integer MASSBAL = 1;
  };
  
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
  
  virtual void RegisterDataOnMesh(ManagedGroup * const meshBodies) override;
  

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void ImplicitStepSetup( real64 const& time_n,
                                  real64 const& dt,
                                  DomainPartition * const domain,
                                  systemSolverInterface::EpetraBlockSystem * const blockSystem ) override;


  virtual void AssembleSystem( DomainPartition * const domain,
                               systemSolverInterface::EpetraBlockSystem * const blockSystem,
                               real64 const time_n,
                               real64 const dt ) override;

  virtual real64
  CalculateResidualNorm(systemSolverInterface::EpetraBlockSystem const *const blockSystem,
                        DomainPartition *const domain) override;
    
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
   * @param well the well containing all the primary and dependent fields
   */
  void UpdateComponentFraction( Well * well );

  /**
   * @brief Recompute phase volume fractions (saturations) from constitutive and primary variables
   * @param well the well containing all the primary and dependent fields
   */
  void UpdatePhaseVolumeFraction( Well * well );

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param well the well containing all the primary and dependent fields
   */
  void UpdateFluidModel( Well * well );
  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param well the well containing all the primary and dependent fields
   */
  void UpdateRelPermModel( Well * well );

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param well the well containing all the primary and dependent fields
   */
  void UpdateState( Well * well );

  
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
   * @brief assembles the accumulation terms for all well elements
   * @param domain the physical domain object
   * @param jacobian the entire jacobian matrix of the system
   * @param residual the entire residual of the system
   * @param time_n previous time value
   * @param dt time step
   */
  void AssembleAccumulationTerms( DomainPartition * const domain,
                                  Epetra_FECrsMatrix * const jacobian,
                                  Epetra_FEVector * const residual,
                                  real64 const time_n,
                                  real64 const dt );

  /**
   * @brief assembles the flux terms for all connections between well elements
   * @param domain the physical domain object
   * @param jacobian the entire jacobian matrix of the system
   * @param residual the entire residual of the system
   * @param time_n previous time value
   * @param dt time step
   */
  void AssembleFluxTerms( DomainPartition * const domain,
                          Epetra_FECrsMatrix * const jacobian,
                          Epetra_FEVector * const residual,
                          real64 const time_n,
                          real64 const dt );

  /**
   * @brief assembles the perforation rate terms 
   * @param domain the physical domain object
   * @param jacobian the entire jacobian matrix of the system
   * @param residual the entire residual of the system
   * @param time_n previous time value
   * @param dt time step
   */

  void AssembleSourceTerms( DomainPartition * const domain,
                            Epetra_FECrsMatrix * const jacobian,
                            Epetra_FEVector * const residual,
                            real64 const time_n,
                            real64 const dt );

  /**
   * @brief assembles the volume balance terms for all well elements
   * @param domain the physical domain object
   * @param jacobian the entire jacobian matrix of the system
   * @param residual the entire residual of the system
   * @param time_n previous time value
   * @param dt time step
   */
  void AssembleVolumeBalanceTerms( DomainPartition * const domain,
                                   Epetra_FECrsMatrix * const jacobian,
                                   Epetra_FEVector * const residual,
                                   real64 const time_n,
                                   real64 const dt );

  /**
   * @brief assembles the momentum at all connections except the first global connection
   * @param domain the physical domain object
   * @param jacobian the entire jacobian matrix of the system
   * @param residual the entire residual of the system
   */
  void FormPressureRelations( DomainPartition * const domain,
                              Epetra_FECrsMatrix * const jacobian,
                              Epetra_FEVector * const residual );

  /**
   * @brief assembles the control equation for the first global connection
   * @param domain the physical domain object
   * @param jacobian the entire jacobian matrix of the system
   * @param residual the entire residual of the system
   */
  void FormControlEquation( DomainPartition * const domain,
                            Epetra_FECrsMatrix * const jacobian,
                            Epetra_FEVector * const residual );
  
  /**
   * @brief set the sparsity pattern for the linear system
   * @param domain the domain partition
   * @param sparsity the sparsity pattern matrix
   */
  void SetSparsityPattern( DomainPartition const * const domain,
                           Epetra_FECrsGraph * const sparsity,
                           globalIndex firstWellElemDofNumber,
                           localIndex numDofPerResElement ) override;

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
  void SetNumRowsAndTrilinosIndices( DomainPartition const * const domain,
                                     localIndex & numLocalRows,
                                     globalIndex & numGlobalRows,
                                     localIndex offset ) override;
  
  /**@}*/

  struct viewKeyStruct : WellSolverBase::viewKeyStruct
  {
    // degrees of freedom numbers on the well elements
    static constexpr auto dofNumberString = "segmentLocalDofNumber";
    
    // inputs
    static constexpr auto temperatureString = "segmentTemperature";
    static constexpr auto useMassFlagString = "useMass";

    static constexpr auto relPermNameString  = "segmentRelPermName";
    static constexpr auto resRelPermIndexString  = "elementRelPermIndex";
    
    // primary solution field
    static constexpr auto pressureString               = "segmentPressure";
    static constexpr auto deltaPressureString          = "deltaSegmentPressure";
    static constexpr auto globalCompDensityString      = "segmentGlobalCompDensity";
    static constexpr auto deltaGlobalCompDensityString = "deltaSegmentGlobalCompDensity";
    static constexpr auto mixtureConnRateString        = "segmentMixtureConnectionRate";
    static constexpr auto deltaMixtureConnRateString   = "deltaSegmentMixtureConnectionRate";

    // saturations
    static constexpr auto phaseVolumeFractionString = "segmentPhaseVolumeFraction";
    static constexpr auto dPhaseVolumeFraction_dPressureString = "dSegmentPhaseVolumeFraction_dPres";
    static constexpr auto dPhaseVolumeFraction_dGlobalCompDensityString = "dSegmentPhaseVolumeFraction_dComp";

    // intermediate values for constitutive model input
    static constexpr auto globalCompFractionString                     = "segmentGlobalComponentFraction";
    static constexpr auto dGlobalCompFraction_dGlobalCompDensityString = "dSegmentGlobalComponentFraction_dGlobalCompDensity";
    
    // perforation rates and derivatives
    static constexpr auto compPerforationRateString = "compPerforationRate";
    static constexpr auto dCompPerforationRate_dPresString = "dCompPerforationRate_dPres";
    static constexpr auto dCompPerforationRate_dCompString = "dCompPerforationRate_dComp";
    static constexpr auto sumCompPerforationRateString = "sumCompPerforationRate";
    
    using ViewKey = dataRepository::ViewKey;

    // degrees of freedom numbers of the well elements
    ViewKey dofNumber = { dofNumberString };
    
    // inputs
    ViewKey temperature = { temperatureString };
    ViewKey useMassFlag = { useMassFlagString };

    ViewKey relPermName  = { relPermNameString };
    ViewKey resRelPermIndex = { resRelPermIndexString };
    
    // primary solution field
    ViewKey pressure               = { pressureString };
    ViewKey deltaPressure          = { deltaPressureString };
    ViewKey globalCompDensity      = { globalCompDensityString };
    ViewKey deltaGlobalCompDensity = { deltaGlobalCompDensityString };
    ViewKey mixtureRate            = { mixtureConnRateString };
    ViewKey deltaMixtureRate       = { deltaMixtureConnRateString };
    
    // saturation
    ViewKey phaseVolFrac        = { phaseVolumeFractionString };
    ViewKey dPhaseVolFrac_dPres = { dPhaseVolumeFraction_dPressureString };
    ViewKey dPhaseVolFrac_dComp = { dPhaseVolumeFraction_dGlobalCompDensityString };
    
    // global composition to input injection stream
    ViewKey globalComponentFrac        = { globalCompFractionString };
    ViewKey dGlobalComponentFrac_dComp = { dGlobalCompFraction_dGlobalCompDensityString };

    // perforation rates
    ViewKey compPerforationRate        = { compPerforationRateString };
    ViewKey dCompPerforationRate_dPres = { dCompPerforationRate_dPresString };
    ViewKey dCompPerforationRate_dComp = { dCompPerforationRate_dCompString };
    ViewKey sumCompPerforationRate     = { sumCompPerforationRateString };
    
  } viewKeysCompMultiphaseWell;

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysCompMultiphaseWell;

protected:
  virtual void InitializePreSubGroups( ManagedGroup * const rootGroup ) override;

  virtual void InitializePostInitialConditions_PreSubGroups( dataRepository::ManagedGroup * const rootGroup ) override;


private:
  
  /**
   * @brief Backup current values of all constitutive fields that participate in the accumulation term
   * @param domain the domain containing the mesh and fields
   */
  void BackupFields( DomainPartition * const domain );

  /**
   * @brief Setup stored reservoir views into domain data for the current step
   */
  void ResetViews( DomainPartition * const domain ) override;

  /**
   * @brief Compute the perforation rates for this well
   * @param well the well with its perforations
   */
  void ComputeAllPerforationRates( Well * well );

  /**
   * @brief Save all the rates and pressures in the well for reporting purposes
   * @param well the well with its perforations
   */
  void RecordWellData( Well * well );

  /**
   * @brief Initialize all the primary and secondary variables in all the wells
   * @param domain the domain containing the well manager to access individual wells
   */
  void InitializeWells( DomainPartition * const domain );
  
  /**
   * @brief Check if the controls are viable; if not, switch the controls
   * @param domain the domain containing the well manager to access individual wells
   */
  void CheckWellControlSwitch( DomainPartition * const domain );
  
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

  /// index of the rel perm constitutive model in the flow solver
  localIndex m_resRelPermIndex;
  
  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> m_resDofNumber;
  
  /// views into reservoir primary variable fields

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_resPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_deltaResPressure;

  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_resGlobalCompDensity;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_deltaResGlobalCompDensity;

  /// views into other reservoir variable fields

  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_resCompFrac;
  ElementRegionManager::ElementViewAccessor<arrayView3d<real64>> m_dResCompFrac_dCompDens;

  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_resPhaseVolFrac;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_dResPhaseVolFrac_dPres;
  ElementRegionManager::ElementViewAccessor<arrayView3d<real64>> m_dResPhaseVolFrac_dCompDens;

  /// views into reservoir material fields

  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> m_resPhaseFrac;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> m_dResPhaseFrac_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> m_dResPhaseFrac_dComp;

  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> m_resPhaseDens;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> m_dResPhaseDens_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> m_dResPhaseDens_dComp;

  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> m_resPhaseVisc;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> m_dResPhaseVisc_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> m_dResPhaseVisc_dComp;

  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> m_resPhaseCompFrac;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> m_dResPhaseCompFrac_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView5d<real64>> m_dResPhaseCompFrac_dComp;

  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_resTotalDens;

  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> m_resPhaseRelPerm;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> m_dResPhaseRelPerm_dPhaseVolFrac;

};
 
} // namespace geosx


#endif //SRC_COMPONENTS_CORE_SRC_WELLS_COMPOSITIONALMULTIPHASEWELL_HPP_
