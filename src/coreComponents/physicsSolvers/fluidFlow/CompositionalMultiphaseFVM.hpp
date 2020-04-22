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
 * @file CompositionalMultiphaseFVM.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFVM_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFVM_HPP_

#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"

namespace geosx
{

namespace dataRepository
{
class Group;

namespace keys
{
string const compositionalMultiphaseFVM = "CompositionalMultiphaseFVM";
}
}

namespace constitutive
{
class MultiFluidBase;
}

/**
 * @class CompositionalMultiphaseFVM
 *
 * A compositional multiphase solver
 * using only cell-centered variables
 * works with both TPFA and MPFA
 */
class CompositionalMultiphaseFVM : public CompositionalMultiphaseBase
{
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

  /**
   * @brief name of the solver in the object catalog
   * @return string that contains the catalog name to generate a new object through the object catalog.
   */
  static string CatalogName() { return dataRepository::keys::compositionalMultiphaseFVM; }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void
  SetupDofs( DomainPartition const * const domain,
             DofManager & dofManager ) const override;

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

  /**
   * @brief assembles the flux terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual
  void AssembleFluxTerms( real64 const time_n,
                          real64 const dt,
                          DomainPartition const * const domain,
                          DofManager const * const dofManager,
                          ParallelMatrix * const matrix,
                          ParallelVector * const rhs ) override;

  /**@}*/

  struct viewKeyStruct : CompositionalMultiphaseBase::viewKeyStruct
  {} viewKeysCompMultiphaseFVM;

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {} groupKeysCompMultiphaseFVM;

private:

  /**
   * @brief Function to perform the Application of Dirichlet type BC's
   * @param time current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  void ApplyDirichletBC_implicit( real64 const time,
                                  real64 const dt,
                                  DofManager const * const dofManager,
                                  DomainPartition * const domain,
                                  ParallelMatrix * const matrix,
                                  ParallelVector * const rhs );

  /**
   * @brief Function to perform the application of a flux boundary condition
   * @param time current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  void ApplySourceFluxBC( real64 const time,
                          real64 const dt,
                          DofManager const * const dofManager,
                          DomainPartition * const domain,
                          ParallelMatrix * const matrix,
                          ParallelVector * const rhs );

  // no data needed here, see CompositionalMultiphaseBase

};


} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFVM_HPP_
