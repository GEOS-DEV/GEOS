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
 * @file ImmiscibleMultiphaseFlow.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_IMMISCIBLEMULTIPHASEFLOW_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_IMMISCIBLEMULTIPHASEFLOW_HPP_

#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
//#include "physicsSolvers/fluidFlow/kernels/immiscibleMultiphase/ImmiscibleMultiphaseKernels.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"  // For GravityDensityScheme
namespace geos
{

namespace immiscibleFlowUtilities
{
  
enum class ScalingType : integer
{
  Global,         // Scale the Newton update with a unique scaling factor
  Local           // Scale the Newton update locally (modifies the Newton direction)
};

enum class ScalingFactorType : integer
{
  MaxVariation,
  TrustRegion,
  TrustRegionFlux
};

ENUM_STRINGS( ScalingType,
              "Global",
              "Local" );

ENUM_STRINGS( ScalingFactorType,
              "MaxVariation",
              "TrustRegion",
              "TrustRegionFlux" );

} // namespace immiscibleFlowUtilities

using namespace immiscibleFlowUtilities;

//START_SPHINX_INCLUDE_00
/**
 * @class ImmiscibleMultiphaseFlow
 *
 * An Immiscible multiphase flow solver
 */
class ImmiscibleMultiphaseFlow : public FlowSolverBase
{
public:

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  ImmiscibleMultiphaseFlow( const string & name,
                            Group * const parent );

  /// deleted default constructor
  ImmiscibleMultiphaseFlow() = delete;

  /// deleted copy constructor
  ImmiscibleMultiphaseFlow( ImmiscibleMultiphaseFlow const & ) = delete;

  /// default move constructor
  ImmiscibleMultiphaseFlow( ImmiscibleMultiphaseFlow && ) = default;

  /// deleted assignment operator
  ImmiscibleMultiphaseFlow & operator=( ImmiscibleMultiphaseFlow const & ) = delete;

  /// deleted move operator
  ImmiscibleMultiphaseFlow & operator=( ImmiscibleMultiphaseFlow && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~ImmiscibleMultiphaseFlow() override = default;
  /**
   * @brief name of the solver in the object catalog
   * @return string that contains the catalog name to generate a new object through the object catalog.
   */
  static string catalogName() { return "ImmiscibleMultiphaseFlow"; }
  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual void registerDataOnMesh( Group & meshBodies ) override final;

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

  virtual real64
  calculateResidualNorm( real64 const & time_n,
                         real64 const & dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

  virtual real64
  scalingForSystemSolution( DomainPartition & domain,
                            DofManager const & dofManager,
                            arrayView1d< real64 > const & localSolution,
                            arrayView1d< real64 const > const & localResidual,
                            real64 const dt,
                            real64 const residualNorm,
                            integer const newtonIter ) override;

  real64
  scalingForSystemSolutionTrustRegion( DomainPartition & domain,
                                       DofManager const & dofManager,
                                       arrayView1d< real64 > const & localSolution,
                                       arrayView1d< real64 const > const & localResidual,
                                       real64 const dt,
                                       real64 const residualNorm,
                                       integer const newtonIter );

  virtual bool
  checkSystemSolution( DomainPartition & domain,
                       DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       real64 const dt,
                       DomainPartition & domain ) override;

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

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

  void updateFluidState( ElementSubRegionBase & subRegion ) const;

  void updateVolumeConstraint( ElementSubRegionBase & subRegion ) const;

  virtual void saveConvergedState( ElementSubRegionBase & subRegion ) const override final;

  virtual void updateState( DomainPartition & domain ) override final;

  void
  updateSolutionField( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       DomainPartition & domain );

  /**
   * @brief Getter for the number of fluid phases
   * @return the number of phases
   */
  integer numFluidPhases() const { return m_numPhases; }

  /**
   * @brief assembles the accumulation and volume balance terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   */
  void assembleAccumulationTerm( DomainPartition & domain,
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
                     arrayView1d< real64 > const & localRhs ) const;

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

  void applySourceFluxBC( real64 const time,
                          real64 const dt,
                          DofManager const & dofManager,
                          DomainPartition & domain,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs ) const;
  /**
   * @brief function to set the next time step size
   * @param[in] currentDt the current time step size
   * @param[in] domain the domain object
   * @return the prescribed time step size
   */

  virtual real64 setNextDtBasedOnStateChange( real64 const & currentDt,
                                              DomainPartition & domain ) override;

  virtual void initializePreSubGroups() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  virtual void initializeFluidState( MeshLevel & mesh, string_array const & regionNames ) override;

  /**
   * @brief Function to update fluid mass
   * @param subRegion subregion that contains the fields
   */
  void updatePhaseMass( ElementSubRegionBase & subRegion ) const;

  struct viewKeyStruct : public FlowSolverBase::viewKeyStruct
  {
    // inputs
    static constexpr char const * capPressureNamesString() { return "capPressureNames"; }
    static constexpr char const * relPermNamesString() { return "relPermNames"; }
    static constexpr char const * elemDofFieldString() { return "elemDofField"; }
    static constexpr char const * elemDofUpdateFieldString() { return "elemDofUpdateField"; }

    // density averaging scheme
    static constexpr char const * gravityDensitySchemeString() { return "gravityDensityScheme"; }

    // scaling scheme
    static constexpr char const * scalingTypeString() { return "scalingType"; }
    static constexpr char const * scalingFactorTypeString() { return "scalingFactorType"; }
    static constexpr char const * maxAbsolutePresChangeString() { return "maxAbsolutePressureChange"; }
    static constexpr char const * maxRelativePresChangeString() { return "maxRelativePressureChange"; }
    static constexpr char const * maxAbsoluteSatChangeString() { return "maxAbsoluteSaturationChange"; }
    static constexpr char const * maxRelativeSatChangeString() { return "maxRelativeSaturationChange"; }

    // trust region parameters
    static constexpr char const * trustRegionMinNewtonIterString() { return "trustRegionMinNewtonIter"; }
    static constexpr char const * trustRegionMinGradientString() { return "trustRegionMinGradient"; } 
    static constexpr char const * trustRegionMaxIterString() { return "trustRegionMaxIter"; }
    static constexpr char const * trustRegionMinPotentialDiffString() { return "trustRegionMinPotentialDiff"; }
    static constexpr char const * trustRegionMinDerivativeString() { return "trustRegionMinDerivative"; }
    static constexpr char const * trustRegionMinKinkFactorString() { return "trustRegionMinKinkFactor"; }
    static constexpr char const * trustRegionMinInfFactorString() { return "trustRegionMinInfFactor"; }
    static constexpr char const * trustRegionKinkFactorDeltaString() { return "trustRegionKinkFactorDelta"; }
    static constexpr char const * trustRegionRelResThresString() { return "trustRegionRelResidualThreshold"; }
    static constexpr char const * trustRegionAbsResThresString() { return "trustRegionAbsResidualThreshold"; }
    static constexpr char const * trustRegionUseAccumString() { return "trustRegionUseAccum"; }

    // time stepping controls
    static constexpr char const * solutionChangeScalingFactorString() { return "solutionChangeScalingFactor"; }
    static constexpr char const * targetRelativePresChangeString() { return "targetRelativePressureChangeInTimeStep"; }
    static constexpr char const * targetPhaseVolFracChangeString() { return "targetPhaseVolFractionChangeInTimeStep"; }

    // nonlinear solver parameters    
    static constexpr char const * useTotalMassEquationString() { return "useTotalMassEquation"; }
    static constexpr char const * allowOutOfBoundPressureString() { return "allowOutOfBoundPressure"; }
    static constexpr char const * allowOutOfBoundSatString() { return "allowOutOfBoundSaturation"; }
    static constexpr char const * allowLocalSatChoppingString() { return "allowLocalSatChopping"; }
    static constexpr char const * allowLocalPresChoppingString() { return "allowLocalPresChopping"; }
    static constexpr char const * minScalingFactorString() { return "minScalingFactor"; }
  };


private:

  virtual void postInputInitialization() override;

  /**
   * @brief Update all relevant fluid models using current values of pressure and phase volume fraction
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
   * @brief Recompute phase mobility from constitutive and primary variables
   * @param dataGroup the group storing the required field
   */
  void updatePhaseMobility( ObjectManagerBase & dataGroup ) const;

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

  /**
   * @brief Utility function to chop saturations that lie outside of physical bounds
   * @param[in] domain the domain object
   */
  void chopOutOfBoundPhaseVolFrac ( DomainPartition & domain );

  /**
   * @brief Utility function to chop pressures that lie outside of physical bounds
   * @param[in] domain the domain object
   */
  void chopOutOfBoundPressure( DomainPartition & domain );

/**
   * @brief Utility function to avoid pressure updates that lead outside of physical bounds
   * @param[in] domain the domain object
   * @param[in] dofManager the dof manager
   * @param[in] localSolution the local solution vector
   */
  void avoidOutOfBoundPressure( DomainPartition & domain, 
                                DofManager const & dofManager,
                                arrayView1d< real64 > const & localSolution );

  /**
   * @brief Utility function to avoid phase volume fraction updates that lead outside of physical bounds
   * @param[in] domain the domain object
   * @param[in] dofManager the dof manager
   * @param[in] localSolution the local solution vector
   */
  void avoidOutOfBoundPhaseVolFrac( DomainPartition & domain, 
                                    DofManager const & dofManager,
                                    arrayView1d< real64 > const & localSolution );

   /**
   * @brief Utility function to reset the local scaling factors
   * @param[in] domain the domain object
   */
  void resetLocalScalingFactors( DomainPartition & domain );

  /**
   * @brief Utility function to check the maximum gradient of the update vector
   * @param[in] domain the domain partition
   * @param[in] dofManager the dof manager
   * @param[in] localSolution the local solution vector
   */
  real64 checkMaxGradient( DomainPartition & domain, 
                           DofManager const & dofManager,
                           arrayView1d< real64 const > const & localSolution );

  /// the max number of fluid phases
  integer m_numPhases;

  /// flag to determine whether or not to apply capillary pressure
  integer m_hasCapPressure;

  /// flag to determine whether or not to use total velocity formulation
  integer m_useTotalMassEquation;

  /// flag to determine whether to use the flux or the residual inflection algorithm
  integer m_fluxInflection;

  /// flag to determine whether to allow negative pressures
  integer m_allowOutOfBoundPressure;

  /// flag to determine whether to chop negative pressures
  integer m_allowPresChopping;

  /// flag to determine whether to allow out of bounds saturation
  integer m_allowOutOfBoundSaturation;

  /// flag to determine whether to chop saturations that lie outside of physical bounds
  integer m_allowSatChopping;

  /// maximum (absolute) pressure change in a Newton iteration
  real64 m_maxAbsolutePresChange;

  /// maximum (absolute) saturation change in a Newton iteration
  real64 m_maxAbsoluteSatChange;

  /// maximum (relative) change in pressure in a Newton iteration
  real64 m_maxRelativePresChange;

  /// maximum (relative) change in saturation in a Newton iteration
  real64 m_maxRelativeSatChange;

  /// scheme for density treatment in gravity
  GravityDensityScheme m_gravityDensityScheme;

  /// scaling type for solution update  
  ScalingType m_scalingType;

  /// scaling type for solution update  
  ScalingFactorType m_scalingFactorType;

  /// minimum value of the scaling factor
  real64 m_minScalingFactor;  

  /// target (relative) change in pressure in a time step
  real64 m_targetRelativePresChange;

  /// target (absolute) change in phase volume fraction in a time step
  real64 m_targetPhaseVolFracChange;

  /// damping factor for solution change targets
  real64 m_solutionChangeScalingFactor;

  /// type of scaling to applying in current Newton iteration
  ScalingType m_currentScaling;

  /// flag to indicate stagnation and turn off physical chopping
  bool m_stagnation;

  /// previous residual norms
  real64 m_prevResidualNorm = 0.0;
  real64 m_prevResidualNorm2 = 0.0;

public:

  /// trust region parameters
  struct TrustRegionParameters
  {
    /// minimum number of Newton iterations before applying trust region
    integer minNewtonIter;

    /// minimum gradient to enable trust region solver
    real64 minGradient;

    /// maximum number of nonlinear iterations
    integer maxIter;

    /// minimum potential difference to apply a damping factor
    real64 dPhiMin;

    /// minimum directional second derivative to apply a damping factor
    real64 d2RMin;
   
    /// minimum discontinuity damping factor
    real64 minKinkFactor;

    /// minimum inflection damping factor
    real64 minInfFactor;

    /// stretching factor applied to damping factor to allow for a small crossing of discontinuities
    real64 kinkFactorDelta;

    /// minimum relative residual threshold for applying damping factor
    real64 relResThres;

    /// minimum absolute residual threshold for applying damping factor
    real64 absResThres;

    /// flag on whether to use accumulation term
    integer useAccum;

  } m_trustRegionParams;

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
void ImmiscibleMultiphaseFlow::applyFieldValue( real64 const & time_n,
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
                                 getName(), time_n+dt, fs.getCatalogName(), fs.getName(),
                                 setName, targetGroup.getName(), fs.getScale(), numTargetElems ) );
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

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_IMMISCIBLEMULTIPHASEFLOW_HPP_
