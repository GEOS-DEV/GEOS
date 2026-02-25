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
 * @file CompositionalMultiphaseWell.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELL_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELL_HPP_

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/relativePermeability/Layouts.hpp"

#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"

#include "physicsSolvers/fluidFlow/wells/WellConstraintsBase.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
namespace geos
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
class CompositionalMultiphaseWell : public WellControls
{
public:

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

  /// deleted move constructor
  CompositionalMultiphaseWell( CompositionalMultiphaseWell && ) = delete;

  /// deleted assignment operator
  CompositionalMultiphaseWell & operator=( CompositionalMultiphaseWell const & ) = delete;

  /// deleted move operator
  CompositionalMultiphaseWell & operator=( CompositionalMultiphaseWell && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~CompositionalMultiphaseWell() override = default;


  virtual void registerWellDataOnMesh( WellElementSubRegion & subRegion ) override;
  /**
   * @defgroup WellManager Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   * The "Well" versions apply to individual well subRegions, whereas the others apply to all wells
   */
  /**@{*/
  /**
   *   * @brief Initialize well for the beginning of a simulation or restart
   *   @param domain the domain
   *   @param mesh the mesh level
   *   @param subRegion the well subRegion
   */
  virtual void initializeWell( DomainPartition & domain, MeshLevel & mesh, WellElementSubRegion & subRegion, real64 const & time_n ) override;

  virtual void initializeWellPostInitialConditionsPreSubGroups( WellElementSubRegion & subRegion ) override;

  virtual bool isCompositional() const override { return true; }


  /**
   * @copydoc WellControls::assembleWellAccumulationTerms()
   */
  virtual void assembleWellAccumulationTerms( real64 const & time,
                                              real64 const & dt,
                                              WellElementSubRegion & subRegion,
                                              DofManager const & dofManager,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs ) override;
  /**
   * @copydoc WellControls::assembleWellPressureRelations()
   */
  virtual void assembleWellPressureRelations( real64 const & time_n,
                                              real64 const & dt,
                                              WellElementSubRegion const & subRegion,
                                              DofManager const & dofManager,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs ) override;

  /**
   * @copydoc WellControls::assembleWellConstraintTerms()
   */
  virtual void assembleWellConstraintTerms( real64 const & time_n,
                                            real64 const & dt,
                                            WellElementSubRegion const & subRegion,
                                            DofManager const & dofManager,
                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                            arrayView1d< real64 > const & localRhs ) override;
  /**
   * @copydoc WellControls::computeWellPerforationRates()
   */
  virtual void computeWellPerforationRates( real64 const & time_n,
                                            real64 const & GEOS_UNUSED_PARAM( dt ),
                                            ElementRegionManager & elemManager,
                                            WellElementSubRegion & subRegion ) override;

  /**
   * @copydoc WellControls::assembleFluxTerms()
   */
  virtual void assembleWellFluxTerms( real64 const & time,
                                      real64 const & dt,
                                      WellElementSubRegion const & subRegion,
                                      DofManager const & dofManager,
                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                      arrayView1d< real64 > const & localRhs ) override;
  /**@}*/
  /**
   * @defgroup Well Interface Functions - required by WellManager and WellNewtonSolver
   *
   * These functions provide the primary interface that is required for derived classes
   * The "Well" versions apply to individual well subRegions
   */
  /**@{*/

  /**
   * @copydoc WellControls::calculateResidualNorm()
   */

  virtual array1d< real64 >
  calculateLocalWellResidualNorm( real64 const & time_n,
                                  real64 const & dt,
                                  NonlinearSolverParameters const & nonlinearSolverParameters,
                                  WellElementSubRegion const & subRegion,
                                  DofManager const & dofManager,
                                  arrayView1d< real64 const > const & localRhs ) override;


  virtual real64
  calculateWellResidualNorm( real64 const & time_n,
                             real64 const & dt,
                             NonlinearSolverParameters const & nonlinearSolverParameters,
                             WellElementSubRegion const & subRegion,
                             DofManager const & dofManager,
                             arrayView1d< real64 const > const & localRhs ) override;

  /**
   * @copydoc WellControls::scalingForSystemSolution()
   */
  real64 scalingForLocalSystemSolution ( WellElementSubRegion & subRegion,
                                         DofManager const & dofManager,
                                         real64 & maxDeltaPres,
                                         real64 & maxDeltaCompDens,
                                         real64 & maxDeltaTemp,
                                         real64 & minPresScalingFactor,
                                         real64 & minCompDensScalingFactor,
                                         real64 & minTempScalingFactor,
                                         arrayView1d< real64 const > const & localSolution );

  virtual real64 scalingForWellSystemSolution( WellElementSubRegion & subRegion,
                                               DofManager const & dofManager,
                                               arrayView1d< real64 const > const & localSolution ) override;

  /**
   * @copydoc WellControls::checkSystemSolution()
   */
  virtual bool
  checkWellSystemSolution( WellElementSubRegion & subRegion,
                           DofManager const & dofManager,
                           arrayView1d< real64 const > const & localSolution,
                           real64 const scalingFactor ) override;

  /**
   * @copydoc WellControls::applyWellSystemSolution()
   */

  virtual void
  applyWellSystemSolution( DofManager const & dofManager,
                           arrayView1d< real64 const > const & localSolution,
                           real64 const scalingFactor,
                           real64 const dt,
                           DomainPartition & domain,
                           MeshLevel & mesh,
                           WellElementSubRegion & subRegion ) override;

  virtual void applyWellBoundaryConditions( real64 const time_n,
                                            real64 const dt,
                                            ElementRegionManager & elemManager,
                                            WellElementSubRegion & subRegion,
                                            DofManager const & dofManager,
                                            arrayView1d< real64 > const & localRhs,
                                            CRSMatrixView< real64, globalIndex const > const & localMatrix ) override;


  virtual void resetStateToBeginningOfStep( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion ) override;

  virtual void implicitStepSetup( real64 const & time_n,
                                  real64 const & GEOS_UNUSED_PARAM( dt ),
                                  ElementRegionManager & elemManager,
                                  WellElementSubRegion & subRegion ) override;

  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        WellElementSubRegion const & subRegion ) override;

  virtual void printRates( real64 const & time_n,
                           real64 const & dt,
                           WellElementSubRegion const & subRegion ) override;

  virtual real64 updateSubRegionState( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion ) override;

  /**@}*/

  /**
   * @brief Recompute global component fractions from primary variables (component densities)
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  void updateGlobalComponentFraction( WellElementSubRegion & subRegion ) const;

  /**
   * @brief Recompute the volumetric rates that are used in the well constraints
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  void updateVolRatesForConstraint( WellElementSubRegion const & subRegion );

  /**
   * @brief Recompute the current BHP pressure
   * @param subRegion the well subregion containing all the primary and dependent fields
   * @param targetIndex the targetIndex of the subRegion
   */
  void updateBHPForConstraint( WellElementSubRegion & subRegion );

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param subRegion the well subregion containing all the primary and dependent fields
   * @param targetIndex the targetIndex of the subRegion
   */
  void updateFluidModel( WellElementSubRegion & subRegion );

  /**
   * @brief Update well separator using current values of pressure and composition at the reference
   * element
   * @param elemManager the element region manager

   */
  void updateSeparator( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion );

  /**
   * @brief  Calculate well rates at reference element
   * @param subRegion the well subregion containing all the primary and dependent fields
   * @param targetIndex the targetIndex of the subRegion
   */

  void calculateReferenceElementRates( WellElementSubRegion & subRegion );

  /**
   * @brief Recompute phase volume fractions (saturations) from constitutive and primary variables
   * @param subRegion the well subregion containing all the primary and dependent fields
   * @param targetIndex the targetIndex of the subRegion
   */
  real64 updatePhaseVolumeFraction( WellElementSubRegion & subRegion ) const;

  /**
   * @brief Recompute total mass densities from mass density and phase volume fractions
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  void updateTotalMassDensity( WellElementSubRegion & subRegion ) const;

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  virtual real64 updateWellState( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion ) override;


  virtual string wellElementDofName() const override { return viewKeyStruct::dofFieldString(); }

  virtual string resElementDofName() const override { return CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString(); }

  virtual localIndex numFluidComponents() const override { return m_numComponents; }

  virtual localIndex numFluidPhases() const override { return m_numPhases; }

  integer useTotalMassEquation() const { return m_useTotalMassEquation; }

  /**
   * @brief Sets all the negative component densities (if any) to zero.
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  void chopNegativeDensities( WellElementSubRegion & subRegion );


  struct viewKeyStruct : WellControls::viewKeyStruct
  {
    static constexpr char const * dofFieldString() { return "wellVars"; }

    // inputs

    static constexpr char const * useMassFlagString() { return CompositionalMultiphaseBase::viewKeyStruct::useMassFlagString(); }

    static constexpr char const * useTotalMassEquationString() { return CompositionalMultiphaseBase::viewKeyStruct::useTotalMassEquationString(); }

    static constexpr char const * maxCompFracChangeString() { return CompositionalMultiphaseBase::viewKeyStruct::maxCompFracChangeString(); }

    static constexpr char const * maxRelativePresChangeString() { return "maxRelativePressureChange"; }

    static constexpr char const * maxAbsolutePresChangeString() { return "maxAbsolutePressureChange"; }

    static constexpr char const * maxRelativeCompDensChangeString() { return "maxRelativeCompDensChange"; }

    static constexpr char const * maxRelativeTempChangeString() { return "maxRelativeTemperatureChange"; }

    static constexpr char const * allowLocalCompDensChoppingString() { return CompositionalMultiphaseBase::viewKeyStruct::allowLocalCompDensChoppingString(); }



  } viewKeysCompMultiphaseWell;

protected:

  virtual void postInputInitialization() override;

  virtual void initializePostSubGroups() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  void saveState( WellElementSubRegion & subRegion );
  virtual void postRestartInitialization( ) override;

  /**
   * @brief Checks fluild model compatibility and validity
   * @param[in] fluid the fluid to check
   * @param[in] referenceFluid the reference fluid model
   * @detail
   * This function will produce an error if one of the well constitutive models
   * is incompatible with the corresponding models in reservoir
   * regions connected to that particular well.
   */
  void validateFluidModel( constitutive::MultiFluidBase const & fluid, constitutive::MultiFluidBase const & referenceFluid )const;
  /**
   * @brief Make sure that the well constraints are compatible
   * @param time_n the time at the beginning of the time step
   * @param dt the time step dt
   * @param subRegion the well subRegion
   * @param elemManager the element manager
   */
  virtual void validateWellConstraints( real64 const & time_n,
                                        real64 const & dt,
                                        WellElementSubRegion const & subRegion ) override;

  /**
   * @brief Create well separator
   */
  virtual void createSeparator( WellElementSubRegion & subRegion ) override;



  virtual bool solveWHPConstraint( real64 const & time_n,
                                   real64 const & dt,
                                   integer const cycleNumber,
                                   integer const coupledIterationNumber,
                                   DomainPartition & domain,
                                   MeshLevel & mesh,
                                   ElementRegionManager & elemManager,
                                   WellElementSubRegion & subRegion )override;
private:

  virtual void setConstitutiveNames( ElementSubRegionBase & subRegion ) const override;


  /// flag indicating whether mass or molar formulation should be used
  integer m_useMass;

  /// flag indicating whether total mass equation should be used
  integer m_useTotalMassEquation;

  /// maximum (absolute) change in a component fraction between two Newton iterations
  real64 m_maxCompFracChange;

  /// maximum (relative) change in pressure between two Newton iterations
  real64 m_maxRelativePresChange;

  /// maximum (absolute) change in pressure between two Newton iterations
  real64 m_maxAbsolutePresChange;

  /// maximum (relative) change in component density between two Newton iterations
  real64 m_maxRelativeCompDensChange;

  /// maximum (relative) change in temperature in a Newton iteration
  real64 m_maxRelativeTempChange;

  /// minimum value of the scaling factor obtained by enforcing maxCompFracChange
  real64 m_minScalingFactor;

  /// flag indicating whether local (cell-wise) chopping of negative compositions is allowed
  integer m_allowCompDensChopping;

};

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELL_HPP_
