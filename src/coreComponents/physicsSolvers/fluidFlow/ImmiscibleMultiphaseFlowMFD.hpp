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
class ImmiscibleMultiphaseFlowMFD : public ImmiscibleMultiphaseFlow
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

  // update mobility to provide total mobility field expected by hybrid kernels
  void updatePhaseMobility( ObjectManagerBase & dataGroup ) const override;

  void initializePreSubGroups() override;
  void initializePostInitialConditionsPreSubGroups() override;

private:
  void computeTotalMobility( ObjectManagerBase & dataGroup ) const;

  // region filter used in hybrid flux assembly
  SortedArray< localIndex > m_regionFilter;
  real64 m_areaRelTol; // tolerance for transmissibility (like single-phase)
};

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_IMMISCIBLEMULTIPHASEFLOWMFD_HPP_
