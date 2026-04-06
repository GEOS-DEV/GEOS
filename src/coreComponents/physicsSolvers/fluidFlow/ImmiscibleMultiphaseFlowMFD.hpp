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
 * @file ImmiscibleMultiphaseFlowMFD.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_IMMISCIBLEMULTIPHASEFLOWMFD_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_IMMISCIBLEMULTIPHASEFLOWMFD_HPP_

#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/fluidFlow/ImmiscibleMultiphaseFlow.hpp"
#include "physicsSolvers/fluidFlow/kernels/immiscibleMultiphase/ImmiscibleMultiphaseMFDKernels.hpp"

namespace geos
{

/**
 * @class ImmiscibleMultiphaseFlowMFD
 *
 * Immiscible two-phase flow solver using a mixed hybrid (mimetic) discretization
 * for the pressure equation (cell + face pressures) while reusing the standard
 * cell-centered transport treatment for phase volume fractions.
 *
 * This is an initial structural implementation that assembles a hybrid flux
 * contribution only into the pressure equation and the face constraints.
 */
class ImmiscibleMultiphaseFlowMFD : public FlowSolverBase
{
public:
  ImmiscibleMultiphaseFlowMFD( const string & name, Group * const parent );
  ImmiscibleMultiphaseFlowMFD() = delete;
  ImmiscibleMultiphaseFlowMFD( ImmiscibleMultiphaseFlowMFD const & ) = delete;
  ImmiscibleMultiphaseFlowMFD( ImmiscibleMultiphaseFlowMFD && ) = default;
  ImmiscibleMultiphaseFlowMFD & operator=( ImmiscibleMultiphaseFlowMFD const & ) = delete;
  ImmiscibleMultiphaseFlowMFD & operator=( ImmiscibleMultiphaseFlowMFD && ) = delete;
  virtual ~ImmiscibleMultiphaseFlowMFD() override = default;

  static string catalogName() { return "ImmiscibleMultiphaseFlowMFD"; }
  string getCatalogName() const override { return catalogName(); }

  // data registration
  void registerDataOnMesh( Group & meshBodies ) override;

  // dofs (add face pressure dofs)
  void setupDofs( DomainPartition const & domain, DofManager & dofManager ) const override;

  // time step hooks
  void implicitStepSetup( real64 const & time_n, real64 const & dt, DomainPartition & domain ) override;
  void implicitStepComplete( real64 const & time, real64 const & dt, DomainPartition & domain ) override;
  void resetStateToBeginningOfStep( DomainPartition & domain ) override;

  // assembly
  void assembleSystem( real64 const time_n,
                       real64 const dt,
                       DomainPartition & domain,
                       DofManager const & dofManager,
                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                       arrayView1d< real64 > const & localRhs ) override;

  void assembleFluxTermsHybrid( real64 const dt,
                                DomainPartition const & domain,
                                DofManager const & dofManager,
                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                arrayView1d< real64 > const & localRhs );

  // boundary conditions (extend pressure BCs to face pressures)
  void applyBoundaryConditions( real64 const time_n,
                                real64 const dt,
                                DomainPartition & domain,
                                DofManager const & dofManager,
                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                arrayView1d< real64 > const & localRhs ) override;

  // residual norm (add face equation contribution)
  real64 calculateResidualNorm( real64 const & time_n,
                                real64 const & dt,
                                DomainPartition const & domain,
                                DofManager const & dofManager,
                                arrayView1d< real64 const > const & localRhs ) override;

  // solution application (add face pressures)
  void applySystemSolution( DofManager const & dofManager,
                            arrayView1d< real64 const > const & localSolution,
                            real64 const scalingFactor,
                            real64 const dt,
                            DomainPartition & domain ) override;

  // state update hook (override base to avoid error)
  virtual void updateState( DomainPartition & domain ) override;

  // multiphase specific public hooks
  integer numFluidPhases() const { return m_numPhases; }

  struct viewKeyStruct : public FlowSolverBase::viewKeyStruct
  {
    static constexpr char const * elemDofFieldString() { return "elemDofField"; }
    static constexpr char const * capPressureNamesString() { return "capillary_pressure"; }
    static constexpr char const * relPermNamesString() { return "relative_permeability"; }
    static constexpr char const * useTotalMassEquationString() { return "useTotalMassEquation"; }
    static constexpr char const * dependentPhaseIndexString() { return "dependentPhaseIndex"; }
  };

  void initializePreSubGroups() override;
  void initializePostInitialConditionsPreSubGroups() override;

private:
  // Add override to register constitutive model names (fluid, relperm) and call base for solid/permeability
  virtual void setConstitutiveNames( ElementSubRegionBase & subRegion ) const override;

  // helper methods migrated from ImmiscibleMultiphaseFlow
  void updateFluidState( ElementSubRegionBase & subRegion ) const;
  void updatePhaseMass( ElementSubRegionBase & subRegion ) const;
  void updatePhaseMobility( ObjectManagerBase & dataGroup ) const;
  void updateVolumeConstraint( ElementSubRegionBase & subRegion ) const;
  void updateRelPermModel( ObjectManagerBase & dataGroup ) const;
  void updateCapPressureModel( ObjectManagerBase & dataGroup ) const;
  void updateFluidModel( ObjectManagerBase & dataGroup ) const;
  void applyDirichletBC( real64 const time_n,
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
  void assembleAccumulationTerm( DomainPartition & domain,
                                 DofManager const & dofManager,
                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                 arrayView1d< real64 > const & localRhs ) const;
  void assembleFluxTerms( real64 const dt,
                          DomainPartition const & domain,
                          DofManager const & dofManager,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs ) const;

  // Dedup helpers (no behavior change)
  void copyPressureToPrev( ElementSubRegionBase & subRegion ) const;
  void copyPhaseStateToPrev( ElementSubRegionBase & subRegion ) const;
  void restoreStateFromPrev( ElementSubRegionBase & subRegion ) const;
  void copyFacePressureToPrev( MeshLevel & mesh ) const;
  void restoreFacePressureFromPrev( MeshLevel & mesh ) const;

  // multiphase state
  integer m_numPhases;
  bool m_hasCapPressure;
  integer m_useTotalMassEquation;
  integer m_dependentPhaseIndex; // 0 or 1, indicates which phase saturation is dependent (s_dep = 1 - s_ind)
  // time stepping targets (simplified defaults)
  real64 m_targetRelativePresChange;
  real64 m_targetPhaseVolFracChange;
  real64 m_solutionChangeScalingFactor;
  GravityDensityScheme m_gravityDensityScheme;
  // region filter & tolerance already present
  SortedArray< localIndex > m_regionFilter;
  real64 m_areaRelTol;

  // Small helper: independent phase index (0 or 1) when one phase is dependent
  GEOS_HOST_DEVICE inline integer indepPhaseIndex() const { return 1 - m_dependentPhaseIndex; }
};

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_IMMISCIBLEMULTIPHASEFLOWMFD_HPP_
