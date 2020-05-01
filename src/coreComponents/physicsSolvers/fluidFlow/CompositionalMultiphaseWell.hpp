/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompositionalMultiphaseWell.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEWELLSOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEWELLSOLVER_HPP_

#include "WellSolverBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFlow.hpp"

namespace geosx
{

namespace dataRepository
{
class Group;
}

namespace constitutive
{
class ConstitutiveManager;
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
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  CompositionalMultiphaseWell( const string & name,
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

  virtual void RegisterDataOnMesh( Group * const meshBodies ) override;


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

  /**@}*/

  /**
   * @brief Recompute mixture densities using current values of pressure and composition
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  void UpdateMixtureDensity( WellElementSubRegion & subRegion, localIndex const targetIndex );

  /**
   * @brief Recompute component fractions from primary variables (component densities)
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  void UpdateComponentFraction( WellElementSubRegion & subRegion ) const;

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  void UpdateFluidModel( WellElementSubRegion & subRegion, localIndex const targetIndex );

  /**
   * @brief Recompute phase volume fractions (saturations) from constitutive and primary variables
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  void UpdatePhaseVolumeFraction( WellElementSubRegion & subRegion, localIndex const targetIndex ) const;

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  virtual void UpdateState( WellElementSubRegion & subRegion, localIndex const targetIndex ) override;

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

  arrayView1d< string const > const & relPermModelNames() const { return m_relPermModelNames; }

  struct viewKeyStruct : WellSolverBase::viewKeyStruct
  {
    static constexpr auto dofFieldString = "compositionalWellVars";

    // inputs
    static constexpr auto temperatureString = "wellTemperature";
    static constexpr auto useMassFlagString = "useMass";

    static constexpr auto relPermNamesString  = "relPermNames";

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
    static constexpr auto dPhaseVolumeFraction_dGlobalCompDensityString =
      CompositionalMultiphaseFlow::viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString;

    // mixture density
    static constexpr auto mixtureDensityString = "wellElementMixtureDensity";
    static constexpr auto dMixtureDensity_dPressureString = "dWellElementMixtureDensity_dPres";
    static constexpr auto dMixtureDensity_dGlobalCompDensityString = "dWellElementMixtureDensity_dComp";

    // global component fractions
    static constexpr auto globalCompFractionString = CompositionalMultiphaseFlow::viewKeyStruct::globalCompFractionString;
    static constexpr auto dGlobalCompFraction_dGlobalCompDensityString =
      CompositionalMultiphaseFlow::viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString;

    // perforation rates and derivatives
    static constexpr auto compPerforationRateString = "compPerforationRate";
    static constexpr auto dCompPerforationRate_dPresString = "dCompPerforationRate_dPres";
    static constexpr auto dCompPerforationRate_dCompString = "dCompPerforationRate_dComp";
  } viewKeysCompMultiphaseWell;

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {} groupKeysCompMultiphaseWell;

protected:

  virtual void PostProcessInput() override;

  virtual void InitializePreSubGroups( Group * const rootGroup ) override;

  virtual void InitializePostInitialConditions_PreSubGroups( Group * const rootGroup ) override;

  /**
   * @brief Checks constitutive models for consistency
   * @param meshLevel reference to the mesh
   * @param cm        reference to the global constitutive model manager
   *
   * This function will produce an error if one of the well constitutive models
   * (fluid, relperm) is incompatible with the corresponding models in reservoir
   * regions connected to that particular well.
   */
  void ValidateConstitutiveModels( MeshLevel const & meshLevel, constitutive::ConstitutiveManager const & cm ) const;

  /**
   * @brief Checks injection streams for validity (compositions sum to one)
   * @param meshLevel reference to the mesh
   */
  void ValidateInjectionStreams( MeshLevel const & meshLevel ) const;

private:

  /**
   * @brief Compute all the perforation rates for this well
   * @param well the well with its perforations
   */
  void ComputePerforationRates( WellElementSubRegion & subRegion, localIndex const targetIndex );

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
  void ResizeFields( WellElementSubRegion & subRegion );

  /**
   * @brief Save all the rates and pressures in the well for reporting purposes
   * @param well the well with its perforations
   */
  void RecordWellData( WellElementSubRegion const & subRegion );

  /// the max number of fluid phases
  localIndex m_numPhases;

  /// the number of fluid components
  localIndex m_numComponents;

  /// the (uniform) temperature
  real64 m_temperature;

  /// flag indicating whether mass or molar formulation should be used
  integer m_useMass;

  /// list of relative permeability model names per target region
  array1d< string > m_relPermModelNames;


  /// views into reservoir primary variable fields

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_resPressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_deltaResPressure;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_resGlobalCompDensity;

  /// views into other reservoir variable fields

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_resPhaseVolFrac;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_dResPhaseVolFrac_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_dResPhaseVolFrac_dCompDens;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_dResCompFrac_dCompDens;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_resPhaseMob;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_dResPhaseMob_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_dResPhaseMob_dCompDens;

  /// views into reservoir material fields

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_resPhaseDens;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_resPhaseVisc;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_dResPhaseVisc_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 > > m_dResPhaseVisc_dComp;

  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 > > m_resPhaseCompFrac;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 > > m_dResPhaseCompFrac_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView5d< real64 > > m_dResPhaseCompFrac_dComp;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_resPhaseRelPerm;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 > > m_dResPhaseRelPerm_dPhaseVolFrac;

};

} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEWELL_HPP_
