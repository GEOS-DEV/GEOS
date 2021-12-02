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
 * @file CompositionalMultiphaseBase.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASE_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASE_HPP_

#include "common/DataLayouts.hpp"
#include "constitutive/fluid/layouts.hpp"
#include "constitutive/relativePermeability/layouts.hpp"
#include "constitutive/capillaryPressure/layouts.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"

namespace geosx
{

//START_SPHINX_INCLUDE_00
/**
 * @class CompositionalMultiphaseBase
 *
 * A compositional multiphase solver
 */
class CompositionalMultiphaseBase : public FlowSolverBase
{
public:

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  CompositionalMultiphaseBase( const string & name,
                               Group * const parent );

  /// deleted default constructor
  CompositionalMultiphaseBase() = delete;

  /// deleted copy constructor
  CompositionalMultiphaseBase( CompositionalMultiphaseBase const & ) = delete;

  /// default move constructor
  CompositionalMultiphaseBase( CompositionalMultiphaseBase && ) = default;

  /// deleted assignment operator
  CompositionalMultiphaseBase & operator=( CompositionalMultiphaseBase const & ) = delete;

  /// deleted move operator
  CompositionalMultiphaseBase & operator=( CompositionalMultiphaseBase && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~CompositionalMultiphaseBase() override = default;

//START_SPHINX_INCLUDE_01

  virtual void registerDataOnMesh( Group & meshBodies ) override;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual real64
  solverStep( real64 const & time_n,
              real64 const & dt,
              integer const cycleNumber,
              DomainPartition & domain ) override;

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  assembleSystem( real64 const time_n,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual void
  applyBoundaryConditions( real64 const time_n,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  virtual void
  solveSystem( DofManager const & dofManager,
               ParallelMatrix & matrix,
               ParallelVector & rhs,
               ParallelVector & solution ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  /**
   * @brief Recompute component fractions from primary variables (component densities)
   * @param dataGroup the group storing the required fields
   */
  void updateComponentFraction( Group & dataGroup ) const;

  /**
   * @brief Recompute phase volume fractions (saturations) from constitutive and primary variables
   * @param dataGroup the group storing the required fields
   */
  void updatePhaseVolumeFraction( Group & dataGroup, localIndex const targetIndex ) const;

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param dataGroup the group storing the required fields
   */
  void updateFluidModel( Group & dataGroup, localIndex const targetIndex ) const;

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param castedRelPerm the group storing the required fields
   */
  void updateRelPermModel( Group & castedRelPerm, localIndex const targetIndex ) const;

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param castedCapPres the group storing the required fields
   */
  void updateCapPressureModel( Group & castedCapPres, localIndex const targetIndex ) const;

  /**
   * @brief Recompute phase mobility from constitutive and primary variables
   * @param domain the domain containing the mesh and fields
   */
  virtual void updatePhaseMobility( Group & dataGroup, localIndex const targetIndex ) const = 0;

  void updateFluidState( Group & dataGroup, localIndex const targetIndex ) const;

  virtual void updateState( DomainPartition & domain ) override final;

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
   * @brief assembles the accumulation and volume balance terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   */
  void assembleAccumulationAndVolumeBalanceTerms( DomainPartition & domain,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief assembles the flux terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void
  assembleFluxTerms( real64 const dt,
                     DomainPartition const & domain,
                     DofManager const & dofManager,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs ) const = 0;


  /**@}*/

  arrayView1d< string const > relPermModelNames() const { return m_relPermModelNames; }

  arrayView1d< string const > capPresModelNames() const { return m_capPressureModelNames; }

  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {
    static constexpr char const * elemDofFieldString() { return "compositionalVariables"; }

    // inputs

    static constexpr char const * inputTemperatureString() { return "temperature"; }

    static constexpr char const * temperatureString() { return "temperature"; }

    static constexpr char const * deltaTemperatureString() { return "deltaTemperature"; }

    static constexpr char const * useMassFlagString() { return "useMass"; }

    static constexpr char const * isothermalFlagString()  { return "isothermal"; }

    static constexpr char const * computeCFLNumbersString() { return "computeCFLNumbers"; }

    static constexpr char const * relPermNamesString() { return "relPermNames"; }

    static constexpr char const * capPressureNamesString() { return "capPressureNames"; }

    static constexpr char const * maxCompFracChangeString() { return "maxCompFractionChange"; }

    static constexpr char const * allowLocalCompDensChoppingString() { return "allowLocalCompDensityChopping"; }

    static constexpr char const * facePressureString() { return "facePressure"; }

    static constexpr char const * bcPressureString() { return "bcPressure"; }

    static constexpr char const * globalCompDensityString() { return "globalCompDensity"; }

    static constexpr char const * deltaGlobalCompDensityString() { return "deltaGlobalCompDensity"; }

    // intermediate values for constitutive model input
    static constexpr char const * globalCompFractionString() { return "globalCompFraction"; }

    static constexpr char const * dGlobalCompFraction_dGlobalCompDensityString() { return "dGlobalCompFraction_dGlobalCompDensity"; }

    static constexpr char const * phaseVolumeFractionString() { return "phaseVolumeFraction"; }

    static constexpr char const * dPhaseVolumeFraction_dPressureString() { return "dPhaseVolumeFraction_dPressure"; }

    static constexpr char const * dPhaseVolumeFraction_dGlobalCompDensityString() { return "dPhaseVolumeFraction_dGlobalCompDensity"; }

    static constexpr char const * dPhaseVolumeFraction_dTemperatureString() { return "dPhaseVolumeFraction_dTemperature"; }

    // intermediate values for mobilities
    static constexpr char const * phaseMobilityString() { return "phaseMobility"; }

    static constexpr char const * dPhaseMobility_dPressureString() { return "dPhaseMobility_dPressure"; }

    static constexpr char const * dPhaseMobility_dGlobalCompDensityString() { return "dPhaseMobility_dGlobalCompDensity"; }

    static constexpr char const * dPhaseMobility_dTemperatureString() { return "dPhaseMobility_dTemperature"; }

    // intermediate values for CFL number computation and actual cell CFL numbers
    static constexpr char const * phaseOutfluxString() { return "phaseOutflux"; }

    static constexpr char const * componentOutfluxString() { return "componentOutflux"; }

    static constexpr char const * phaseCFLNumberString() { return "phaseCFLNumber"; }

    static constexpr char const * componentCFLNumberString() { return "componentCFLNumber"; }

    // these are used to store last converged time step values
    static constexpr char const * phaseVolumeFractionOldString() { return "phaseVolumeFractionOld"; }

    static constexpr char const * phaseDensityOldString() { return "phaseDensityOld"; }

    static constexpr char const * totalDensityOldString() { return "totalDensityOld"; }

    static constexpr char const * phaseComponentFractionOldString() { return "phaseComponentFractionOld"; }

    static constexpr char const * phaseMobilityOldString() { return "phaseMobilityOld"; }

    static constexpr char const * phaseInternalEnergyOldString() { return "phaseInternalEnergyOld"; }

    // these are allocated on faces for BC application until we can get constitutive models on faces
    static constexpr char const * phaseViscosityString() { return "phaseViscosity"; }

    static constexpr char const * phaseRelativePermeabilityString() { return "phaseRelativePermeability"; }

    static constexpr char const * phaseCapillaryPressureString() { return "phaseCapillaryPressure"; }
  };

  /**
   * @brief Initialize all variables from initial conditions
   * @param domain the domain containing the mesh and fields
   *
   * Initialize all variables from initial conditions. This calculating primary variable values
   * from prescribed intermediate values (i.e. global densities from global fractions)
   * and any applicable hydrostatic equilibration of the domain
   */
  void initializeFluidState( MeshLevel & mesh );

  /**
   * @brief Compute the hydrostatic equilibrium using the compositions and temperature input tables
   */
  void computeHydrostaticEquilibrium();

  /**
   * @brief Backup current values of all constitutive fields that participate in the accumulation term
   * @param domain the domain containing the mesh and fields
   */
  void backupFields( MeshLevel & mesh ) const;

  /**
   * @brief Function to perform the Application of Dirichlet type BC's
   * @param time current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param localMatrix local system matrix
   * @param localRhs local system right-hand side vector
   */
  void applyDirichletBC( real64 const time,
                         real64 const dt,
                         DofManager const & dofManager,
                         DomainPartition & domain,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief Apply source flux boundary conditions to the system
   * @param time current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param localMatrix local system matrix
   * @param localRhs local system right-hand side vector
   */
  void applySourceFluxBC( real64 const time,
                          real64 const dt,
                          DofManager const & dofManager,
                          DomainPartition & domain,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief Apply aquifer boundary conditions to the system
   * @param time current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param localMatrix local system matrix
   * @param localRhs local system right-hand side vector
   */
  virtual void applyAquiferBC( real64 const time,
                               real64 const dt,
                               DofManager const & dofManager,
                               DomainPartition & domain,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) const = 0;

  /**
   * @brief Sets all the negative component densities (if any) to zero.
   * @param domain the physical domain object
   */
  void chopNegativeDensities( DomainPartition & domain );

  virtual void initializePostInitialConditionsPreSubGroups() override;

protected:

  virtual void postProcessInput() override;

  virtual void initializePreSubGroups() override;

  /**
   * @brief Checks constitutive models for consistency
   * @param[in] cm reference to the global constitutive model manager
   */
  void validateConstitutiveModels( constitutive::ConstitutiveManager const & cm ) const;

  /**
   * @brief Checks aquifer boundary condition for consistency
   * @param[in] cm reference to the global constitutive model manager
   */
  void validateAquiferBC( constitutive::ConstitutiveManager const & cm ) const;

  /**
   * @brief Initialize the aquifer boundary condition (gravity vector, water phase index)
   * @param[in] cm reference to the global constitutive model manager
   */
  void initializeAquiferBC( constitutive::ConstitutiveManager const & cm ) const;

  /**
   * @brief Setup stored views into domain data for the current step
   */
  void resetViews( MeshLevel & mesh ) override;

  /// the max number of fluid phases
  integer m_numPhases;

  /// the number of fluid components
  integer m_numComponents;

  /// the input temperature
  real64 m_inputTemperature;

  /// flag indicating if the problem is isothermal
  integer m_isothermalFlag;

  /// flag indicating whether mass or molar formulation should be used
  integer m_useMass;

  /// flag indicating whether CFL numbers will be computed or not
  integer m_computeCFLNumbers;

  /// name of the rel perm constitutive model
  array1d< string > m_relPermModelNames;

  /// flag to determine whether or not to apply capillary pressure
  integer m_capPressureFlag;

  /// name of the cap pressure constitutive model
  array1d< string > m_capPressureModelNames;

  /// maximum (absolute) change in a component fraction between two Newton iterations
  real64 m_maxCompFracChange;

  /// minimum value of the scaling factor obtained by enforcing maxCompFracChange
  real64 m_minScalingFactor;

  /// flag indicating whether local (cell-wise) chopping of negative compositions is allowed
  integer m_allowCompDensChopping;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_temperature;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_deltaTemperature;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, compflow::USD_COMP_DC > > m_dCompFrac_dCompDens;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const, compflow::USD_PHASE > > m_phaseVolFrac;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const, compflow::USD_PHASE > > m_dPhaseVolFrac_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const, compflow::USD_PHASE > > m_dPhaseVolFrac_dTemp;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, compflow::USD_PHASE_DC > > m_dPhaseVolFrac_dCompDens;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const, compflow::USD_PHASE > > m_phaseMob;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const, compflow::USD_PHASE > > m_dPhaseMob_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const, compflow::USD_PHASE > > m_dPhaseMob_dTemp;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, compflow::USD_PHASE_DC > > m_dPhaseMob_dCompDens;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, constitutive::relperm::USD_RELPERM > > m_phaseRelPerm;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > m_phaseVisc;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > m_phaseDens;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > m_dPhaseDens_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > m_dPhaseDens_dTemp;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > m_dPhaseDens_dComp;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > m_phaseMassDens;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > m_dPhaseMassDens_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > m_dPhaseMassDens_dTemp;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > m_dPhaseMassDens_dComp;

  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > m_phaseCompFrac;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > m_dPhaseCompFrac_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > m_dPhaseCompFrac_dTemp;
  ElementRegionManager::ElementViewAccessor< arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > > m_dPhaseCompFrac_dComp;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > m_phaseCapPressure;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > m_dPhaseCapPressure_dPhaseVolFrac;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > m_phaseEnthalpy;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > m_dPhaseEnthalpy_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > m_dPhaseEnthalpy_dTemp;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > m_dPhaseEnthalpy_dComp;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_totalDensOld;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const, compflow::USD_PHASE > > m_phaseMobOld;

};

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASE_HPP_
