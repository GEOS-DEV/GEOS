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
 * @file SinglePhaseMixedMFD.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEMIXEDMFD_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEMIXEDMFD_HPP_

#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"

namespace geos
{

/**
 * @class SinglePhaseMixedMFD
 *
 * Single-phase flow solver based on the mixed mimetic finite difference formulation:
 * the primary unknowns are the cell pressures and the face mass fluxes, coupled in a
 * saddle-point system. The cell-wise inner product is selected adaptively (TPFA/MFD)
 * through residual-based Global Adaptation indicators.
 */
class SinglePhaseMixedMFD : public SinglePhaseBase
{
public:

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  SinglePhaseMixedMFD( const string & name,
                       Group * const parent );

  /// deleted default constructor
  SinglePhaseMixedMFD() = delete;

  /// deleted copy constructor
  SinglePhaseMixedMFD( SinglePhaseMixedMFD const & ) = delete;

  /// default move constructor
  SinglePhaseMixedMFD( SinglePhaseMixedMFD && ) = default;

  /// deleted assignment operator
  SinglePhaseMixedMFD & operator=( SinglePhaseMixedMFD const & ) = delete;

  /// deleted move operator
  SinglePhaseMixedMFD & operator=( SinglePhaseMixedMFD && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseMixedMFD() override = default;

  /**
   * @brief name of the solver in the object catalog
   * @return string that contains the catalog name to generate a new object through the object catalog.
   */
  static string catalogName()
  { return "SinglePhaseMixedMFD"; }

  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual void registerDataOnMesh( Group & meshBodies ) override;

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  setupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               ParallelVector & rhs,
               ParallelVector & solution,
               bool const setSparsity = true ) override;

  virtual void
  applyBoundaryConditions( real64 const time_n,
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

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       real64 const dt,
                       DomainPartition & domain ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  virtual void
  assembleFluxTerms( real64 const dt,
                     DomainPartition const & domain,
                     DofManager const & dofManager,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs ) override;

  virtual void
  assembleStabilizedFluxTerms( real64 const dt,
                               DomainPartition const & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

  virtual void
  assembleEDFMFluxTerms( real64 const time_n,
                         real64 const dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs,
                         string const & jumpDofKey ) override final;

  virtual void
  applyAquiferBC( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) const override final;

  virtual void
  saveAquiferConvergedState( real64 const & time,
                             real64 const & dt,
                             DomainPartition & domain ) override;

  virtual void initializePreSubGroups() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

private:

  /**
   * @brief Populate the boundary face pressure values from the field specifications.
   * @param[in] time the time at which the boundary values are evaluated
   * @param[in] domain the domain
   */
  void applyFacePressureBCValues( real64 const time,
                                  DomainPartition & domain );

  /**
   * @brief Run the residual-based Global Adaptation pipeline and mark the cells.
   * @param[in] domain the domain
   */
  void computeGlobalAdaptationIndicators( DomainPartition & domain );

  /**
   * @brief Build the per-dof labels used by the stencilFlag-guided three-level MGR strategy:
   *        0 = face flux with exactly-diagonal row (all adjacent cells TPFA-compatible),
   *        1 = face flux adjacent to at least one MFD-compatible cell,
   *        2 = cell pressure.
   * @param[in] domain the domain
   * @param[in] dofManager the dof manager (dof numbers must be finalized)
   */
  void computeMgrPointMarkers( DomainPartition const & domain,
                               DofManager const & dofManager );

  /// relative tolerance used in the mass matrix computations
  real64 m_areaRelTol;

  /// region filter used in flux assembly
  SortedArray< localIndex > m_regionFilter;

};

} /* namespace geos */

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEMIXEDMFD_HPP_
