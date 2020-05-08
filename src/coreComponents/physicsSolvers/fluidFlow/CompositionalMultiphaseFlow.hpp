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
 * @file CompositionalMultiphaseFlow.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FINITEVOLUME_COMPOSITIONALMULTIPHASEFLOW_HPP_
#define GEOSX_PHYSICSSOLVERS_FINITEVOLUME_COMPOSITIONALMULTIPHASEFLOW_HPP_

#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"

namespace geosx
{

namespace dataRepository
{
class Group;

namespace keys
{
string const compositionalMultiphaseFlow = "CompositionalMultiphaseFlow";
}
}
class FieldSpecificationBase;
class FiniteElementBase;
class DomainPartition;

namespace constitutive
{
class MultiFluidBase;
}

//START_SPHINX_INCLUDE_00
/**
 * @class CompositionalMultiphaseFlow
 *
 * A compositional multiphase solver
 */
class CompositionalMultiphaseFlow : public FlowSolverBase
{
public:

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  CompositionalMultiphaseFlow( const string & name,
                               Group * const parent );

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

//START_SPHINX_INCLUDE_01

  /**
   * @brief name of the solver in the object catalog
   * @return string that contains the catalog name to generate a new object through the object catalog.
   */
  static string CatalogName() { return dataRepository::keys::compositionalMultiphaseFlow; }

  virtual void RegisterDataOnMesh( Group * const MeshBodies ) override;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual real64 SolverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition * domain ) override;

  virtual void
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition * const domain,
                     DofManager & dofManager,
                     ParallelMatrix & matrix,
                     ParallelVector & rhs,
                     ParallelVector & solution ) override;

  virtual void
  SetupDofs( DomainPartition const * const domain,
             DofManager & dofManager ) const override;

  virtual void
  AssembleSystem( real64 const time_n,
                  real64 const dt,
                  DomainPartition * const domain,
                  DofManager const & dofManager,
                  ParallelMatrix & matrix,
                  ParallelVector & rhs ) override;

  virtual void
  ApplyBoundaryConditions( real64 const time_n,
                           real64 const dt,
                           DomainPartition * const domain,
                           DofManager const & dofManager,
                           ParallelMatrix & matrix,
                           ParallelVector & rhs ) override;

  virtual real64
  CalculateResidualNorm( DomainPartition const * const domain,
                         DofManager const & dofManager,
                         ParallelVector const & rhs ) override;

  virtual void
  SolveSystem( DofManager const & dofManager,
               ParallelMatrix & matrix,
               ParallelVector & rhs,
               ParallelVector & solution ) override;

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

  /**
   * @brief Recompute component fractions from primary variables (component densities)
   * @param dataGroup the group storing the required fields
   */
  void UpdateComponentFraction( Group & dataGroup ) const;

  /**
   * @brief Recompute phase volume fractions (saturations) from constitutive and primary variables
   * @param dataGroup the group storing the required fields
   */
  void UpdatePhaseVolumeFraction( Group & dataGroup, localIndex const targetIndex ) const;

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param dataGroup the group storing the required fields
   */
  void UpdateFluidModel( Group & dataGroup, localIndex const targetIndex );

  /**
   * @brief Update all relevant solid models using current values of pressure
   * @param dataGroup the group storing the required fields
   */
  void UpdateSolidModel( Group & dataGroup, localIndex const targetIndex );

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param dataGroup the group storing the required fields
   */
  void UpdateRelPermModel( Group & dataGroup, localIndex const targetIndex );

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param dataGroup the group storing the required fields
   */
  void UpdateCapPressureModel( Group & dataGroup, localIndex const targetIndex );

  /**
   * @brief Recompute phase mobility from constitutive and primary variables
   * @param domain the domain containing the mesh and fields
   */
  void UpdatePhaseMobility( Group & dataGroup, localIndex const targetIndex ) const;

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param domain the domain containing the mesh and fields
   */
  void UpdateState( Group & dataGroup, localIndex const targetIndex );

  /**
   * @brief Get the number of fluid components (species)
   * @return the number of components
   */
  localIndex numFluidComponents() const { return m_numComponents; }

  /**
   * @brief Get the number of fluid phases
   * @return the number of phases
   */
  localIndex numFluidPhases() const { return m_numPhases; }

  /**
   * @brief assembles the accumulation terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  void AssembleAccumulationTerms( real64 const time_n,
                                  real64 const dt,
                                  DomainPartition const * const domain,
                                  DofManager const * const dofManager,
                                  ParallelMatrix * const matrix,
                                  ParallelVector * const rhs );

  /**
   * @brief assembles the flux terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  void AssembleFluxTerms( real64 const time_n,
                          real64 const dt,
                          DomainPartition const * const domain,
                          DofManager const * const dofManager,
                          ParallelMatrix * const matrix,
                          ParallelVector * const rhs );

  /**
   * @brief assembles the volume balance terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  void AssembleVolumeBalanceTerms( real64 const time_n,
                                   real64 const dt,
                                   DomainPartition const * const domain,
                                   DofManager const * const dofManager,
                                   ParallelMatrix * const matrix,
                                   ParallelVector * const rhs );

  /**@}*/

  arrayView1d< string const > const & relPermModelNames() const { return m_relPermModelNames; }

  arrayView1d< string const > const & capPresModelNames() const { return m_capPressureModelNames; }

  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {
    static constexpr auto dofFieldString = "compositionalVariables";

    // inputs
    static constexpr auto temperatureString = "temperature";
    static constexpr auto useMassFlagString = "useMass";

    static constexpr auto relPermNamesString  = "relPermNames";
    static constexpr auto capPressureNamesString  = "capPressureNames";

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

    // intermediate values for mobilities
    static constexpr auto phaseMobilityString                     = "phaseMobility";
    static constexpr auto dPhaseMobility_dPressureString          = "dPhaseMobility_dPressure";
    static constexpr auto dPhaseMobility_dGlobalCompDensityString = "dPhaseMobility_dGlobalCompDensity";

    // these are used to store last converged time step values
    static constexpr auto phaseVolumeFractionOldString     = "phaseVolumeFractionOld";
    static constexpr auto phaseDensityOldString            = "phaseDensityOld";
    static constexpr auto phaseComponentFractionOldString  = "phaseComponentFractionOld";
    static constexpr auto porosityOldString                = "porosityOld";

    // these are allocated on faces for BC application until we can get constitutive models on faces
    static constexpr auto phaseViscosityString             = "phaseViscosity";
    static constexpr auto phaseRelativePermeabilityString  = "phaseRelativePermeability";
    static constexpr auto phaseCapillaryPressureString     = "phaseCapillaryPressure";
  } viewKeysCompMultiphaseFlow;

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {} groupKeysCompMultiphaseFlow;

protected:

  virtual void PostProcessInput() override;

  virtual void InitializePreSubGroups( Group * const rootGroup ) override;

  virtual void InitializePostInitialConditions_PreSubGroups( dataRepository::Group * const rootGroup ) override;

  /**
   * @brief Checks constitutive models for consistency
   * @param cm        reference to the global constitutive model manager
   */
  void ValidateConstitutiveModels( constitutive::ConstitutiveManager const & cm ) const;

private:

  /**
   * @brief Resize the allocated multidimensional fields
   * @param domain the domain containing the mesh and fields
   *
   * Resize fields along dimensions 1 and 2 (0 is the size of containing object, i.e. element subregion)
   * once the number of phases/components is known (e.g. component fractions)
   */
  void ResizeFields( MeshLevel & meshLevel );

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
   * @brief Setup stored views into domain data for the current step
   */
  void ResetViews( DomainPartition * const domain ) override;

  /**
   * @brief Function to perform the Application of Dirichlet type BC's
   * @param time current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   */
  void ApplyDirichletBC_implicit( real64 const time,
                                  real64 const dt,
                                  DofManager const * const dofManager,
                                  DomainPartition * const domain,
                                  ParallelMatrix * const matrix,
                                  ParallelVector * const rhs );


  void ApplySourceFluxBC( real64 const time,
                          real64 const dt,
                          DofManager const * const dofManager,
                          DomainPartition * const domain,
                          ParallelMatrix * const matrix,
                          ParallelVector * const rhs );


  /// the max number of fluid phases
  localIndex m_numPhases;

  /// the number of fluid components
  localIndex m_numComponents;

  /// the (uniform) temperature
  real64 m_temperature;

  /// flag indicating whether mass or molar formulation should be used
  integer m_useMass;

  /// name of the rel perm constitutive model
  array1d< string > m_relPermModelNames;

  /// flag to determine whether or not to apply capillary pressure
  integer m_capPressureFlag;

  /// name of the cap pressure constitutive model
  array1d< string > m_capPressureModelNames;


  /// views into primary variable fields

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > >      m_pressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > >      m_deltaPressure;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > >      m_globalCompDensity;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > >      m_deltaGlobalCompDensity;

  /// views into other variable fields

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_compFrac;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_dCompFrac_dCompDens;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_phaseVolFrac;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_dPhaseVolFrac_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_dPhaseVolFrac_dCompDens;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_phaseMob;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_dPhaseMob_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_dPhaseMob_dCompDens;

  /// views into backup fields

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_porosityOld;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_phaseVolFracOld;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_phaseDensOld;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_phaseCompFracOld;

  /// views into material fields

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_pvMult;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_dPvMult_dPres;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_phaseFrac;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_dPhaseFrac_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 > > m_dPhaseFrac_dComp;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_phaseDens;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_dPhaseDens_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 > > m_dPhaseDens_dComp;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_phaseVisc;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_dPhaseVisc_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 > > m_dPhaseVisc_dComp;

  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 > > m_phaseCompFrac;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 > > m_dPhaseCompFrac_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView5d< real64 > > m_dPhaseCompFrac_dComp;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_totalDens;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_phaseRelPerm;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 > > m_dPhaseRelPerm_dPhaseVolFrac;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_phaseCapPressure;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 > > m_dPhaseCapPressure_dPhaseVolFrac;

};


} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FINITEVOLUME_COMPOSITIONALMULTIPHASEFLOW_HPP_
