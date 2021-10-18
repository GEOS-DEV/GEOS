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
 * @file CompositionalMultiphaseFVM.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFVM_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFVM_HPP_

#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"

namespace geosx
{

/**
 * @class CompositionalMultiphaseFVM
 *
 * A compositional multiphase solver
 * using only cell-centered variables
 * works with both TPFA and MPFA
 */
//START_SPHINX_INCLUDE_00
class CompositionalMultiphaseFVM : public CompositionalMultiphaseBase
{
//END_SPHINX_INCLUDE_00
public:

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  CompositionalMultiphaseFVM( const string & name,
                              Group * const parent );

  /// deleted default constructor
  CompositionalMultiphaseFVM() = delete;

  /// deleted copy constructor
  CompositionalMultiphaseFVM( CompositionalMultiphaseFVM const & ) = delete;

  /// default move constructor
  CompositionalMultiphaseFVM( CompositionalMultiphaseFVM && ) = default;

  /// deleted assignment operator
  CompositionalMultiphaseFVM & operator=( CompositionalMultiphaseFVM const & ) = delete;

  /// deleted move operator
  CompositionalMultiphaseFVM & operator=( CompositionalMultiphaseFVM && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~CompositionalMultiphaseFVM() override = default;

//START_SPHINX_INCLUDE_01
  /**
   * @brief name of the solver in the object catalog
   * @return string that contains the catalog name to generate a new object through the object catalog.
   */
  static string catalogName() { return "CompositionalMultiphaseFVM"; }
//END_SPHINX_INCLUDE_01

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

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
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;


  /**@}*/

  virtual void
  assembleFluxTerms( real64 const dt,
                     DomainPartition const & domain,
                     DofManager const & dofManager,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs ) const override;


  virtual void
  updatePhaseMobility( Group & dataGroup, localIndex const targetIndex ) const override;

  virtual void
  applyAquiferBC( real64 const time,
                  real64 const dt,
                  DofManager const & dofManager,
                  DomainPartition & domain,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) const override;


  /**
   * @brief Compute the largest CFL number in the domain
   * @param dt the time step size
   * @param domain the domain containing the mesh and fields
   */
  void
  computeCFLNumbers( real64 const & dt, DomainPartition & domain );


  virtual void initializePreSubGroups() override;

private:

  // no data needed here, see CompositionalMultiphaseBase

};


} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFVM_HPP_
