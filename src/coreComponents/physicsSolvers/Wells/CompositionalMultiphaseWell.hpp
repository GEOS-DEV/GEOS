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

#ifndef GEOSX_PHYSICSSOLVERS_COMPOSITIONALMULTIPHASEWELLSOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_COMPOSITIONALMULTIPHASEWELLSOLVER_HPP_

#include "WellSolverBase.hpp"
#include "constitutive/RelPerm/RelativePermeabilityBase.hpp"
#include "physicsSolvers/FiniteVolume/CompositionalMultiphaseFlow.hpp"

namespace geosx
{

namespace dataRepository
{
class Group;
}
  
namespace constitutive
{
class MultiFluidBase;
}
class WellElementSubRegion;

/**
 * @class CompositionalMultiphaseWell
 *
 * A compositional multiphase well solver
 */
class CompositionalMultiphaseWell : public WellSolverBase
{
public:

  // define the column offset of the derivatives
  struct ColOffset
  {
    static constexpr integer DPRES = 0;
    static constexpr integer DCOMP = 1;
  };

  // define the row offset of the residual equations
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
                               Group * const parent );

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
  static string CatalogName() { return "CompositionalMultiphaseWell"; }
  
  virtual void RegisterDataOnMesh(Group * const meshBodies) override;
  

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/


  virtual real64
  CalculateResidualNorm( DomainPartition const * const domain,
                         DofManager const & dofManager,
                         ParallelVector const & rhs ) override;
    
  virtual bool
  CheckSystemSolution( DomainPartition const * const domain,
                       DofManager const & dofManager,
                       ParallelVector const & solution,
                       real64 const scalingFactor ) override;

  virtual void
  ApplySystemSolution( DofManager const & dofManager,
                       ParallelVector const & solution,
                       real64 const scalingFactor,
                       DomainPartition * const domain ) override;
  
  virtual void 
  ResetStateToBeginningOfStep( DomainPartition * const domain ) override;

  virtual void 
  ImplicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition * const domain ) override;

  virtual void
  SetupDofs( DomainPartition const * const domain,
             DofManager & dofManager ) const override;

  /**@}*/

  /**
   * @brief Recompute mixture densities using current values of pressure and composition
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  void UpdateMixtureDensity( WellElementSubRegion const * const subRegion );

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  virtual void UpdateState( WellElementSubRegion * const subRegion ) override;

  virtual string WellElementDofName() const override { return viewKeyStruct::dofFieldString; }

  virtual string ResElementDofName() const override { return CompositionalMultiphaseFlow::viewKeyStruct::dofFieldString; }

  virtual localIndex NumFluidComponents() const override { return m_numComponents; }

  virtual localIndex NumFluidPhases() const override { return m_numPhases; }

  /**
   * @brief assembles the flux terms for all connections between well elements
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void AssembleFluxTerms( real64 const time_n,
                                  real64 const dt,
                                  DomainPartition const * const domain,
                                  DofManager const * const dofManager,
                                  ParallelMatrix * const matrix,
                                  ParallelVector * const rhs ) override;

  /**
   * @brief assembles the volume balance terms for all well elements
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void AssembleVolumeBalanceTerms( real64 const time_n,
                                           real64 const dt,
                                           DomainPartition const * const domain,
                                           DofManager const * const dofManager,
                                           ParallelMatrix * const matrix,
                                           ParallelVector * const rhs ) override;

  /**
   * @brief assembles the pressure relations at all connections between well elements except at the well head
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void FormPressureRelations( DomainPartition const * const domain,
                                      DofManager const * const dofManager,
                                      ParallelMatrix * const matrix,
                                      ParallelVector * const rhs ) override;

  /**
   * @brief assembles the control equation for the first connection
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void FormControlEquation( DomainPartition const * const domain,
                                    DofManager const * const dofManager,
                                    ParallelMatrix * const matrix,
                                    ParallelVector * const rhs ) override;
  
  struct viewKeyStruct : WellSolverBase::viewKeyStruct
  {
    static constexpr auto dofFieldString = "compositionalWellVars";

    // inputs
    static constexpr auto temperatureString = "wellTemperature";
    static constexpr auto useMassFlagString = "useMass";

    static constexpr auto resRelPermNameString  = "wellRelPermName";
    static constexpr auto resRelPermIndexString = "elementRelPermIndex";
    
    // primary solution field
    static constexpr auto pressureString = CompositionalMultiphaseFlow::viewKeyStruct::pressureString;
    static constexpr auto deltaPressureString = CompositionalMultiphaseFlow::viewKeyStruct::deltaPressureString;
    static constexpr auto globalCompDensityString = CompositionalMultiphaseFlow::viewKeyStruct::globalCompDensityString;
    static constexpr auto deltaGlobalCompDensityString = CompositionalMultiphaseFlow::viewKeyStruct::deltaGlobalCompDensityString;
    static constexpr auto mixtureConnRateString = "wellElementMixtureConnectionRate";
    static constexpr auto deltaMixtureConnRateString = "deltaWellElementMixtureConnectionRate";

    // saturations
    static constexpr auto phaseVolumeFractionString = CompositionalMultiphaseFlow::viewKeyStruct::phaseVolumeFractionString;
    static constexpr auto dPhaseVolumeFraction_dPressureString = CompositionalMultiphaseFlow::viewKeyStruct::dPhaseVolumeFraction_dPressureString;
    static constexpr auto dPhaseVolumeFraction_dGlobalCompDensityString = CompositionalMultiphaseFlow::viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString;

    // mixture density
    static constexpr auto mixtureDensityString = "wellElementMixtureDensity";
    static constexpr auto dMixtureDensity_dPressureString = "dWellElementMixtureDensity_dPres";
    static constexpr auto dMixtureDensity_dGlobalCompDensityString = "dWellElementMixtureDensity_dComp";

    // global component fractions
    static constexpr auto globalCompFractionString = CompositionalMultiphaseFlow::viewKeyStruct::globalCompFractionString;
    static constexpr auto dGlobalCompFraction_dGlobalCompDensityString = CompositionalMultiphaseFlow::viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString;
    
    // perforation rates and derivatives
    static constexpr auto compPerforationRateString = "compPerforationRate";
    static constexpr auto dCompPerforationRate_dPresString = "dCompPerforationRate_dPres";
    static constexpr auto dCompPerforationRate_dCompString = "dCompPerforationRate_dComp";
    
    using ViewKey = dataRepository::ViewKey;

    // inputs
    ViewKey temperature = { temperatureString };
    ViewKey useMassFlag = { useMassFlagString };

    ViewKey resRelPermName  = { resRelPermNameString };
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
    
    // mixture density
    ViewKey mixtureDensity        = { mixtureDensityString };
    ViewKey dMixtureDensity_dPres = { dMixtureDensity_dPressureString };
    ViewKey dMixtureDensity_dComp = { dMixtureDensity_dGlobalCompDensityString };

    // global component fractions
    ViewKey globalComponentFrac        = { globalCompFractionString };
    ViewKey dGlobalComponentFrac_dComp = { dGlobalCompFraction_dGlobalCompDensityString };

    // perforation rates
    ViewKey compPerforationRate        = { compPerforationRateString };
    ViewKey dCompPerforationRate_dPres = { dCompPerforationRate_dPresString };
    ViewKey dCompPerforationRate_dComp = { dCompPerforationRate_dCompString };
    
  } viewKeysCompMultiphaseWell;

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysCompMultiphaseWell;

protected:

  virtual void InitializePreSubGroups( Group * const rootGroup ) override;

  virtual void InitializePostInitialConditions_PreSubGroups( Group * const rootGroup ) override;

private:
  
  /**
   * @brief Setup stored reservoir views into domain data for the current step
   */
  void ResetViews( DomainPartition * const domain ) override;

  /**
   * @brief Initialize all the primary and secondary variables in all the wells
   * @param domain the domain containing the well manager to access individual wells
   */
  void InitializeWells( DomainPartition * const domain ) override;
  
  /**
   * @brief Check if the controls are viable; if not, switch the controls
   * @param domain the domain containing the well manager to access individual wells
   */
  void CheckWellControlSwitch( DomainPartition * const domain ) override;

  /**
   * @brief Resize the allocated multidimensional fields
   * @param well the well for which the fields are resized
   */
  void ResizeFields( WellElementSubRegion * const subRegion );

  /**
   * @brief Save all the rates and pressures in the well for reporting purposes
   * @param well the well with its perforations
   */
  void RecordWellData( WellElementSubRegion const * const subRegion );
  
  /// the max number of fluid phases
  localIndex m_numPhases;

  /// the number of fluid components
  localIndex m_numComponents;
  
  /// the (uniform) temperature
  real64 m_temperature;

  /// flag indicating whether mass or molar formulation should be used
  integer m_useMass;

  /// name of the rel perm constitutive model
  string m_resRelPermName;

  /// index of the rel perm constitutive model in the flow solver
  localIndex m_resRelPermIndex;
  
  /// views into reservoir primary variable fields

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_resPressure;

  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_resGlobalCompDensity;

  /// views into other reservoir variable fields

  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_resPhaseVolFrac;

  /// views into reservoir material fields

  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> m_resPhaseDens;

};
 
} // namespace geosx


#endif //SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_WELLS_COMPOSITIONALMULTIPHASEWELL_HPP_
