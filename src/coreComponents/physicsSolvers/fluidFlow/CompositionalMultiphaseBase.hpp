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

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASE_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASE_HPP_

#include "common/DataLayouts.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/relativePermeability/layouts.hpp"
#include "constitutive/capillaryPressure/layouts.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"

namespace geos
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
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  /**
   * @brief Recompute component fractions from primary variables (component densities)
   * @param dataGroup the group storing the required fields
   */
  void updateComponentFraction( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Recompute phase volume fractions (saturations) from constitutive and primary variables
   * @param dataGroup the group storing the required fields
   */
  real64 updatePhaseVolumeFraction( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param dataGroup the group storing the required fields
   */
  void updateFluidModel( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Update all relevant relperm models using current values of phase volume fraction
   * @param dataGroup the group storing the required fields
   */
  void updateRelPermModel( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Update all relevant capillary pressure models using current values of phase volume fraction
   * @param dataGroup the group storing the required fields
   */
  void updateCapPressureModel( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Update all relevant solid internal energy models using current values of temperature
   * @param dataGroup the group storing the required fields
   */
  void updateSolidInternalEnergyModel( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Recompute phase mobility from constitutive and primary variables
   * @param dataGroup the group storing the required field
   */
  virtual void updatePhaseMobility( ObjectManagerBase & dataGroup ) const = 0;

  real64 updateFluidState( ObjectManagerBase & dataGroup ) const;

  virtual void saveConvergedState( ElementSubRegionBase & subRegion ) const override final;

  virtual void saveIterationState( DomainPartition & domain ) const override final;

  virtual void saveIterationState( ElementSubRegionBase & subRegion ) const override final;

  virtual void updateState( DomainPartition & domain ) override final;

  /**
   * @brief Getter for the number of fluid components (species)
   * @return the number of components
   */
  integer numFluidComponents() const { return m_numComponents; }

  /**
   * @brief Getter for the number of fluid phases
   * @return the number of phases
   */
  integer numFluidPhases() const { return m_numPhases; }

  /**
   * @brief Getter for the name of the reference fluid model name
   * @return the name of the reference fluid
   */
  string referenceFluidModelName() const { return m_referenceFluidModelName; }

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

  /**
   * @brief assembles the flux terms for all cells with pressure jump stabilization
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void
  assembleStabilizedFluxTerms( real64 const dt,
                               DomainPartition const & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) const = 0;
  /**@}*/

  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {
    static constexpr char const * elemDofFieldString() { return "compositionalVariables"; }

    // inputs

    static constexpr char const * inputTemperatureString() { return "temperature"; }
    static constexpr char const * useMassFlagString() { return "useMass"; }
    static constexpr char const * relPermNamesString() { return "relPermNames"; }
    static constexpr char const * capPressureNamesString() { return "capPressureNames"; }
    static constexpr char const * thermalConductivityNamesString() { return "thermalConductivityNames"; }
    static constexpr char const * diffusionNamesString() { return "diffusionNames"; }
    static constexpr char const * dispersionNamesString() { return "dispersionNames"; }


    // time stepping controls

    static constexpr char const * solutionChangeScalingFactorString() { return "solutionChangeScalingFactor"; }
    static constexpr char const * targetRelativePresChangeString() { return "targetRelativePressureChangeInTimeStep"; }
    static constexpr char const * targetRelativeTempChangeString() { return "targetRelativeTemperatureChangeInTimeStep"; }
    static constexpr char const * targetPhaseVolFracChangeString() { return "targetPhaseVolFractionChangeInTimeStep"; }

    // nonlinear solver parameters

    static constexpr char const * maxCompFracChangeString() { return "maxCompFractionChange"; }
    static constexpr char const * maxRelativePresChangeString() { return "maxRelativePressureChange"; }
    static constexpr char const * maxRelativeTempChangeString() { return "maxRelativeTemperatureChange"; }
    static constexpr char const * allowLocalCompDensChoppingString() { return "allowLocalCompDensityChopping"; }
    static constexpr char const * useTotalMassEquationString() { return "useTotalMassEquation"; }
    static constexpr char const * useSimpleAccumulationString() { return "useSimpleAccumulation"; }
    static constexpr char const * minCompDensString() { return "minCompDens"; }

  };

  /**
   * @brief Initialize all variables from initial conditions
   * @param domain the domain containing the mesh and fields
   *
   * Initialize all variables from initial conditions. This calculating primary variable values
   * from prescribed intermediate values (i.e. global densities from global fractions)
   * and any applicable hydrostatic equilibration of the domain
   */
  void initializeFluidState( MeshLevel & mesh, arrayView1d< string const > const & regionNames );

  /**
   * @brief Compute the hydrostatic equilibrium using the compositions and temperature input tables
   */
  void computeHydrostaticEquilibrium();

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
   * @brief Utility function to keep the flow variables during a time step (used in poromechanics simulations)
   * @param[in] keepFlowVariablesConstantDuringInitStep flag to tell the solver to freeze its primary variables during a time step
   * @detail This function is meant to be called by a specific task before/after the initialization step
   */
  void keepFlowVariablesConstantDuringInitStep( bool const keepFlowVariablesConstantDuringInitStep )
  { m_keepFlowVariablesConstantDuringInitStep = keepFlowVariablesConstantDuringInitStep; }

  /**
   * @brief Function to fix the initial state during the initialization step in coupled problems
   * @param[in] time current time
   * @param[in] dt time step
   * @param[in] dofManager degree-of-freedom manager associated with the linear system
   * @param[in] domain the domain
   * @param[in] localMatrix local system matrix
   * @param[in] localRhs local system right-hand side vector
   * @detail This function is meant to be called when the flag m_keepFlowVariablesConstantDuringInitStep is on
   *         The main use case is the initialization step in coupled problems during which we solve an elastic problem for a fixed pressure
   */
  void keepFlowVariablesConstantDuringInitStep( real64 const time,
                                                real64 const dt,
                                                DofManager const & dofManager,
                                                DomainPartition & domain,
                                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                arrayView1d< real64 > const & localRhs ) const;


  /**
   * @brief Sets all the negative component densities (if any) to zero.
   * @param domain the physical domain object
   */
  void chopNegativeDensities( DomainPartition & domain );

  virtual real64 setNextDtBasedOnStateChange( real64 const & currentDt,
                                              DomainPartition & domain ) override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  integer useTotalMassEquation() const { return m_useTotalMassEquation; }

protected:

  virtual void postProcessInput() override;

  virtual void initializePreSubGroups() override;

  /**
   * @brief Utility function that checks the consistency of the constitutive models
   * @param[in] domain the domain partition
   * This function will produce an error if one of the constitutive models
   * (fluid, relperm) is incompatible with the reference fluid model.
   */
  void validateConstitutiveModels( DomainPartition const & domain ) const;

  /**
   * @brief Initialize the aquifer boundary condition (gravity vector, water phase index)
   * @param[in] cm reference to the global constitutive model manager
   */
  void initializeAquiferBC( constitutive::ConstitutiveManager const & cm ) const;

  /**
   * @brief Utility function that encapsulates the call to FieldSpecificationBase::applyFieldValue in BC application
   * @param[in] time_n the time at the beginning of the step
   * @param[in] dt the time step
   * @param[in] mesh the mesh level object
   * @param[in] logMessage the log message issued by the solver if the bc is called
   * @param[in] fieldKey the key of the field specified in the xml file
   * @param[in] boundaryFieldKey the key of the boundary field
   */
  template< typename OBJECT_TYPE >
  void applyFieldValue( real64 const & time_n,
                        real64 const & dt,
                        MeshLevel & mesh,
                        char const logMessage[],
                        string const fieldKey,
                        string const boundaryFieldKey ) const;

  /// the max number of fluid phases
  integer m_numPhases;

  /// the number of fluid components
  integer m_numComponents;

  /// the input temperature
  real64 m_inputTemperature;

  /// flag indicating whether mass or molar formulation should be used
  integer m_useMass;

  /// flag to determine whether or not to apply capillary pressure
  integer m_hasCapPressure;

  /// flag to determine whether or not to apply diffusion
  integer m_hasDiffusion;

  /// flag to determine whether or not to apply dispersion
  integer m_hasDispersion;

  /// flag to freeze the initial state during initialization in coupled problems
  integer m_keepFlowVariablesConstantDuringInitStep;

  /// maximum (absolute) change in a component fraction in a Newton iteration
  real64 m_maxCompFracChange;

  /// maximum (relative) change in pressure in a Newton iteration
  real64 m_maxRelativePresChange;

  /// maximum (relative) change in temperature in a Newton iteration
  real64 m_maxRelativeTempChange;

  /// damping factor for solution change targets
  real64 m_solutionChangeScalingFactor;

  /// target (relative) change in pressure in a time step
  real64 m_targetRelativePresChange;

  /// target (relative) change in temperature in a time step
  real64 m_targetRelativeTempChange;

  /// target (absolute) change in phase volume fraction in a time step
  real64 m_targetPhaseVolFracChange;

  /// minimum value of the scaling factor obtained by enforcing maxCompFracChange
  real64 m_minScalingFactor;

  /// flag indicating whether local (cell-wise) chopping of negative compositions is allowed
  integer m_allowCompDensChopping;

  /// flag indicating whether total mass equation is used
  integer m_useTotalMassEquation;

  /// flag indicating whether simple accumulation form is used
  integer m_useSimpleAccumulation;

  /// minimum allowed global component density
  real64 m_minCompDens;

  /// name of the fluid constitutive model used as a reference for component/phase description
  string m_referenceFluidModelName;

private:

  /**
   * @brief Utility function to validate the consistency of Dirichlet BC input
   * @param[in] domain the domain partition
   * @param[in] time the time at the end of the time step (time_n + dt)
   */
  bool validateDirichletBC( DomainPartition & domain,
                            real64 const time ) const;

  virtual void setConstitutiveNames( ElementSubRegionBase & subRegion ) const override;

};

template< typename OBJECT_TYPE >
void CompositionalMultiphaseBase::applyFieldValue( real64 const & time_n,
                                                   real64 const & dt,
                                                   MeshLevel & mesh,
                                                   char const logMessage[],
                                                   string const fieldKey,
                                                   string const boundaryFieldKey ) const
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  fsManager.apply< OBJECT_TYPE >( time_n + dt,
                                  mesh,
                                  fieldKey,
                                  [&]( FieldSpecificationBase const & fs,
                                       string const & setName,
                                       SortedArrayView< localIndex const > const & lset,
                                       OBJECT_TYPE & targetGroup,
                                       string const & )
  {
    if( fs.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
    {
      globalIndex const numTargetElems = MpiWrapper::sum< globalIndex >( lset.size() );
      GEOS_LOG_RANK_0( GEOS_FMT( logMessage,
                                 getName(), time_n+dt, FieldSpecificationBase::catalogName(),
                                 fs.getName(), setName, targetGroup.getName(), fs.getScale(), numTargetElems ) );
    }

    // Specify the bc value of the field
    fs.applyFieldValue< FieldSpecificationEqual,
                        parallelDevicePolicy<> >( lset,
                                                  time_n + dt,
                                                  targetGroup,
                                                  boundaryFieldKey );
  } );
}


} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASE_HPP_
