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
 * @file SinglePhaseCellCentered.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASECELLCENTERED_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASECELLCENTERED_HPP_

#include "physicsSolvers/fluidFlow/SinglePhaseFlowBase.hpp"

namespace geosx
{

namespace dataRepository
{
class Group;
}
class FieldSpecificationBase;

class FiniteElementBase;

class DomainPartition;

/**
 * @class SinglePhaseCellCentered
 *
 * class to perform a single phase finite volume solve 
 * using only cell-centered variables 
 * works with both TPFA and MPFA
 */
class SinglePhaseCellCentered : public SinglePhaseFlowBase
{
public:
  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  SinglePhaseCellCentered( const std::string & name,
                           Group * const parent );


  /// deleted default constructor
  SinglePhaseCellCentered() = delete;

  /// deleted copy constructor
  SinglePhaseCellCentered( SinglePhaseCellCentered const & ) = delete;

  /// default move constructor
  SinglePhaseCellCentered( SinglePhaseCellCentered && ) = default;

  /// deleted assignment operator
  SinglePhaseCellCentered & operator=( SinglePhaseCellCentered const & ) = delete;

  /// deleted move operator
  SinglePhaseCellCentered & operator=( SinglePhaseCellCentered && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseCellCentered() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName()
  { return "SinglePhaseCellCentered"; }

  /**
   * @defgroup Solver Interface Function
   *
   * This function provides the primary interface that is required for derived classes
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
  virtual void 
  AssembleFluxTerms( real64 const time_n,
                     real64 const dt,
                     DomainPartition const * const domain,
                     DofManager const * const dofManager,
                     ParallelMatrix * const matrix,
                     ParallelVector * const rhs ) override;

  /**@}*/

  struct viewKeyStruct : SinglePhaseFlowBase::viewKeyStruct
  {
  } viewKeysSinglePhaseCellCentered;

  viewKeyStruct & viewKeys()
  { return viewKeysSinglePhaseCellCentered; }

  viewKeyStruct const & viewKeys() const
  { return viewKeysSinglePhaseCellCentered; }

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysSinglePhaseCellCentered;

  groupKeyStruct & groupKeys()
  { return groupKeysSinglePhaseCellCentered; }

  groupKeyStruct const & groupKeys() const
  { return groupKeysSinglePhaseCellCentered; }

private:

  /**
   * @brief Function to perform the application of Dirichlet BCs on faces
   * @param time_n current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  void ApplyFaceDirichletBC_implicit( real64 const time_n,
                                      real64 const dt,
                                      DofManager const * const dofManager,
                                      DomainPartition * const domain,
                                      ParallelMatrix * const matrix,
                                      ParallelVector * const rhs );

  // no data needed here, see SinglePhaseFlowBase

};


} /* namespace geosx */

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASECELLCENTERED_HPP_
