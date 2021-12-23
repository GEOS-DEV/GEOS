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
 * @file CompositionalMultiphaseWell.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELL_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELL_HPP_

#include "constitutive/fluid/layouts.hpp"
#include "constitutive/relativePermeability/layouts.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"

namespace geosx
{

namespace constitutive
{
class ConstitutiveManager;
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
  static string catalogName() { return "CompositionalMultiphaseWell"; }

  virtual void registerDataOnMesh( Group & meshBodies ) override;


  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/


  virtual real64
  calculateResidualNorm( DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

  virtual real64
  scalingForSystemSolution( DomainPartition const & domain,
                            DofManager const & dofManager,
                            arrayView1d< real64 const > const & localSolution ) override;

  virtual bool
  checkSystemSolution( DomainPartition const & domain,
                       DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void
  implicitStepSetup( real64 const & time,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  /**@}*/

  /**
   * @brief Recompute component fractions from primary variables (component densities)
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  void updateComponentFraction( WellElementSubRegion & subRegion ) const;

  /**
   * @brief Recompute the volumetric rates that are used in the well constraints
   * @param subRegion the well subregion containing all the primary and dependent fields
   * @param targetIndex the targetIndex of the subRegion
   */
  void updateVolRatesForConstraint( WellElementSubRegion & subRegion,
                                    localIndex const targetIndex );

  /**
   * @brief Recompute the current BHP pressure
   * @param subRegion the well subregion containing all the primary and dependent fields
   * @param targetIndex the targetIndex of the subRegion
   */
  void updateBHPForConstraint( WellElementSubRegion & subRegion,
                               localIndex const targetIndex );

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param subRegion the well subregion containing all the primary and dependent fields
   * @param targetIndex the targetIndex of the subRegion
   */
  void updateFluidModel( WellElementSubRegion & subRegion, localIndex const targetIndex );

  /**
   * @brief Recompute phase volume fractions (saturations) from constitutive and primary variables
   * @param subRegion the well subregion containing all the primary and dependent fields
   * @param targetIndex the targetIndex of the subRegion
   */
  void updatePhaseVolumeFraction( WellElementSubRegion & subRegion, localIndex const targetIndex ) const;

  /**
   * @brief Recompute total mass densities from mass density and phase volume fractions
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  void updateTotalMassDensity( WellElementSubRegion & subRegion, localIndex const targetIndex ) const;


  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param subRegion the well subregion containing all the primary and dependent fields
   * @param targetIndex the targetIndex of the subRegion
   */
  virtual void updateSubRegionState( WellElementSubRegion & subRegion, localIndex const targetIndex ) override;

  virtual string wellElementDofName() const override { return viewKeyStruct::dofFieldString(); }

  virtual string resElementDofName() const override { return CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString(); }

  virtual localIndex numFluidComponents() const override { return m_numComponents; }

  virtual localIndex numFluidPhases() const override { return m_numPhases; }

  /**
   * @brief assembles the flux terms for all connections between well elements
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void assembleFluxTerms( real64 const time_n,
                                  real64 const dt,
                                  DomainPartition const & domain,
                                  DofManager const & dofManager,
                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                  arrayView1d< real64 > const & localRhs ) override;

  /**
   * @brief assembles the accumulation term for all the well elements
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void assembleAccumulationTerms( DomainPartition const & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs ) override;

  /**
   * @brief assembles the volume balance terms for all well elements
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void assembleVolumeBalanceTerms( DomainPartition const & domain,
                                           DofManager const & dofManager,
                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                           arrayView1d< real64 > const & localRhs ) override;

  /**
   * @brief assembles the pressure relations at all connections between well elements except at the well head
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void assemblePressureRelations( DomainPartition const & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs ) override;


  /**
   * @brief Sets all the negative component densities (if any) to zero.
   * @param domain the physical domain object
   */
  void chopNegativeDensities( DomainPartition & domain );

  /**
   * @brief Backup current values of all constitutive fields that participate in the accumulation term
   * @param mesh reference to the mesh
   */
  void backupFields( MeshLevel & mesh ) const override;


  arrayView1d< string const > relPermModelNames() const { return m_relPermModelNames; }

  struct viewKeyStruct : WellSolverBase::viewKeyStruct
  {
    static constexpr char const * dofFieldString() { return "compositionalWellVars"; }

    // inputs

    static constexpr char const * temperatureString() { return "wellTemperature"; }

    static constexpr char const * deltaTemperatureString() { return "deltaWellTemperature"; }

    static constexpr char const * useMassFlagString() { return CompositionalMultiphaseBase::viewKeyStruct::useMassFlagString(); }

    static constexpr char const * relPermNamesString() { return CompositionalMultiphaseBase::viewKeyStruct::relPermNamesString(); }

    static constexpr char const * maxCompFracChangeString() { return CompositionalMultiphaseBase::viewKeyStruct::maxCompFracChangeString(); }

    static constexpr char const * maxRelativePresChangeString() { return "maxRelativePressureChange"; }

    static constexpr char const * allowLocalCompDensChoppingString() { return CompositionalMultiphaseBase::viewKeyStruct::allowLocalCompDensChoppingString(); }

    // primary solution field
    static constexpr char const * globalCompDensityString() { return extrinsicMeshData::flow::globalCompDensity::key(); }

    static constexpr char const * deltaGlobalCompDensityString() { return extrinsicMeshData::flow::deltaGlobalCompDensity::key(); }

    static constexpr char const * mixtureConnRateString() { return "wellElementMixtureConnectionRate"; }

    static constexpr char const * deltaMixtureConnRateString() { return "deltaWellElementMixtureConnectionRate"; }

    // saturations
    static constexpr char const * phaseVolumeFractionString() { return extrinsicMeshData::flow::phaseVolumeFraction::key(); }

    static constexpr char const * dPhaseVolumeFraction_dPressureString() { return extrinsicMeshData::flow::dPhaseVolumeFraction_dPressure::key(); }

    static constexpr char const * dPhaseVolumeFraction_dGlobalCompDensityString() { return extrinsicMeshData::flow::dPhaseVolumeFraction_dGlobalCompDensity::key(); }

    // global component fractions
    static constexpr char const * globalCompFractionString() { return extrinsicMeshData::flow::globalCompFraction::key(); }

    static constexpr char const * dGlobalCompFraction_dGlobalCompDensityString() { return extrinsicMeshData::flow::dGlobalCompFraction_dGlobalCompDensity::key(); }

    // total mass densities
    static constexpr char const * totalMassDensityString() { return "totalMassDensity"; }

    static constexpr char const * dTotalMassDensity_dPressureString() { return "dTotalMassDensity_dPressure"; }

    static constexpr char const * dTotalMassDensity_dGlobalCompDensityString() { return "dTotalMassDensity_dComp"; }

    // these are used to store last converged time step values
    static constexpr char const * phaseVolumeFractionOldString() { return extrinsicMeshData::flow::phaseVolumeFractionOld::key(); }

    static constexpr char const * phaseDensityOldString() { return extrinsicMeshData::flow::phaseDensityOld::key(); }

    static constexpr char const * totalDensityOldString() { return extrinsicMeshData::flow::totalDensityOld::key(); }

    static constexpr char const * phaseComponentFractionOldString() { return extrinsicMeshData::flow::phaseComponentFractionOld::key(); }


    // perforation rates and derivatives
    static constexpr char const * compPerforationRateString() { return "compPerforationRate"; }

    static constexpr char const * dCompPerforationRate_dPresString() { return "dCompPerforationRate_dPres"; }

    static constexpr char const * dCompPerforationRate_dCompString() { return "dCompPerforationRate_dComp"; }

    // control data
    static constexpr char const * currentBHPString() { return "currentBHP"; }

    static constexpr char const * dCurrentBHP_dPresString() { return "dCurrentBHP_dPres"; }

    static constexpr char const * dCurrentBHP_dCompDensString() { return "dCurrentBHP_dCompDens"; }

    static constexpr char const * currentPhaseVolRateString() { return "currentPhaseVolumetricRate"; }

    static constexpr char const * dCurrentPhaseVolRate_dPresString() { return "dCurrentPhaseVolumetricRate_dPres"; }

    static constexpr char const * dCurrentPhaseVolRate_dCompDensString() { return "dCurrentPhaseVolumetricRate_dCompDens"; }

    static constexpr char const * dCurrentPhaseVolRate_dRateString() { return "dCurrentPhaseVolumetricRate_dRate"; }

    static constexpr char const * currentTotalVolRateString() { return "currentTotalVolumetricRate"; }

    static constexpr char const * dCurrentTotalVolRate_dPresString() { return "dCurrentTotalVolumetricRate_dPres"; }

    static constexpr char const * dCurrentTotalVolRate_dCompDensString() { return "dCurrentTotalVolumetricRate_dCompDens"; }

    static constexpr char const * dCurrentTotalVolRate_dRateString() { return "dCurrentTotalVolumetricRate_dRate"; }

  } viewKeysCompMultiphaseWell;

protected:

  virtual void postProcessInput() override;

  virtual void initializePreSubGroups() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  /**
   * @brief Checks constitutive models for consistency
   * @param meshLevel reference to the mesh
   * @param cm        reference to the global constitutive model manager
   *
   * This function will produce an error if one of the well constitutive models
   * (fluid, relperm) is incompatible with the corresponding models in reservoir
   * regions connected to that particular well.
   */
  void validateConstitutiveModels( MeshLevel const & meshLevel, constitutive::ConstitutiveManager const & cm ) const;

  /**
   * @brief Checks injection streams for validity (compositions sum to one)
   * @param meshLevel reference to the mesh
   */
  void validateInjectionStreams( MeshLevel const & meshLevel ) const;

  /**
   * @brief Make sure that the well constraints are compatible
   * @param meshLevel the mesh level object (to loop over wells)
   * @param fluid the fluid model (to get the target phase index)
   */
  void validateWellConstraints( MeshLevel const & meshLevel, constitutive::MultiFluidBase const & fluid );

private:

  /**
   * @brief Compute all the perforation rates for this well
   * @param subRegion the well element region for which the fields are resized
   * @param targetIndex index of the target region
   */
  void computePerforationRates( WellElementSubRegion & subRegion, localIndex const targetIndex );

  /**
   * @brief Setup stored reservoir views into domain data for the current step
   * @param domain the domain containing the well manager to access individual wells
   */
  void resetViews( DomainPartition & domain ) override;

  /**
   * @brief Initialize all the primary and secondary variables in all the wells
   * @param domain the domain containing the well manager to access individual wells
   */
  void initializeWells( DomainPartition & domain ) override;

  /// the max number of fluid phases
  integer m_numPhases;

  /// the number of fluid components
  integer m_numComponents;

  /// flag indicating whether mass or molar formulation should be used
  integer m_useMass;

  /// list of relative permeability model names per target region
  array1d< string > m_relPermModelNames;

  /// maximum (absolute) change in a component fraction between two Newton iterations
  real64 m_maxCompFracChange;

  /// maximum (relative) change in pressure between two Newton iterations
  real64 m_maxRelativePresChange;

  /// minimum value of the scaling factor obtained by enforcing maxCompFracChange
  real64 m_minScalingFactor;

  /// flag indicating whether local (cell-wise) chopping of negative compositions is allowed
  integer m_allowCompDensChopping;

  /// index of the target phase, used to impose the phase rate constraint
  localIndex m_targetPhaseIndex;

  /// views into reservoir primary variable fields

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_resPres;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_deltaResPres;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_resTemp;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const, compflow::USD_COMP > > m_resCompDens;

  /// views into other reservoir variable fields

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const, compflow::USD_PHASE > > m_resPhaseVolFrac;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const, compflow::USD_PHASE > > m_dResPhaseVolFrac_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, compflow::USD_PHASE_DC > > m_dResPhaseVolFrac_dCompDens;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, compflow::USD_COMP_DC > > m_dResCompFrac_dCompDens;

  /// views into reservoir material fields

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > m_resPhaseDens;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > m_dResPhaseDens_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > m_dResPhaseDens_dComp;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > m_resPhaseMassDens;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > m_resPhaseVisc;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > m_dResPhaseVisc_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > m_dResPhaseVisc_dComp;

  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > m_resPhaseCompFrac;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > m_dResPhaseCompFrac_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > > m_dResPhaseCompFrac_dComp;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, constitutive::relperm::USD_RELPERM > > m_resPhaseRelPerm;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const, constitutive::relperm::USD_RELPERM_DS > > m_dResPhaseRelPerm_dPhaseVolFrac;

};

} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELL_HPP_
