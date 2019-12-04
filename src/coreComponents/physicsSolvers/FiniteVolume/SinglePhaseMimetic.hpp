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
 * @file SinglePhaseMimetic.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FINITEVOLUME_SINGLEPHASEMIMETIC_HPP_
#define GEOSX_PHYSICSSOLVERS_FINITEVOLUME_SINGLEPHASEMIMETIC_HPP_

#include "physicsSolvers/FiniteVolume/SinglePhaseFlowBase.hpp"

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
 * @class SinglePhaseMimetic
 *
 * class to assemble the single-phase flow equations based
 * on a mixed hybrid formulation known as the mimetic method
 */
class SinglePhaseMimetic : public SinglePhaseFlowBase
{
public:
  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  SinglePhaseMimetic( const std::string & name,
                   Group * const parent );


  /// deleted default constructor
  SinglePhaseMimetic() = delete;

  /// deleted copy constructor
  SinglePhaseMimetic( SinglePhaseMimetic const & ) = delete;

  /// default move constructor
  SinglePhaseMimetic( SinglePhaseMimetic && ) = default;

  /// deleted assignment operator
  SinglePhaseMimetic & operator=( SinglePhaseMimetic const & ) = delete;

  /// deleted move operator
  SinglePhaseMimetic & operator=( SinglePhaseMimetic && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseMimetic() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName()
  { return "SinglePhaseMimetic"; }

  virtual void RegisterDataOnMesh( Group * const MeshBodies ) override;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

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

  virtual void
  ResetStateToBeginningOfStep( DomainPartition * const domain ) override;

  virtual void
  ImplicitStepComplete( real64 const & time,
                        real64 const & dt,
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
                     ParallelVector * const rhs );

  /**@}*/


  struct viewKeyStruct : SinglePhaseFlowBase::viewKeyStruct
  {
    // primary face-based field
    static constexpr auto facePressureString      = "facePressure";
    static constexpr auto deltaFacePressureString = "deltaFacePressure";

    // intermediate half-face-based field
    static constexpr auto oneSidedFluxString      = "oneSidedFlux";

  } viewKeysSinglePhaseMimetic;

  viewKeyStruct & viewKeys()
  { return viewKeysSinglePhaseMimetic; }

  viewKeyStruct const & viewKeys() const
  { return viewKeysSinglePhaseMimetic; }

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysSinglePhaseMimetic;

  groupKeyStruct & groupKeys()
  { return groupKeysSinglePhaseMimetic; }

  groupKeyStruct const & groupKeys() const
  { return groupKeysSinglePhaseMimetic; }

  /**
   * @brief Setup stored views into domain data for the current step
   */
  void ResetViews( DomainPartition * const domain ) override;

protected:

  virtual void InitializePostInitialConditions_PreSubGroups( dataRepository::Group * const rootGroup ) override;

private:

  // For each local face, we compute the one-sided volumetric flux
  void ComputeLocalVolumetricOneSidedFluxes( DomainPartition const * const domain, 
                                             localIndex ei );

  // Multiply the one-sided fluxes by density and add mass conservation residual
  void AssembleLocalMassFluxTerms( real64 const time_n,
                                   real64 const dt,
                                   DomainPartition const * const domain,
                                   DofManager const * const dofManager,
                                   ParallelMatrix * const matrix,
                                   ParallelVector * const rhs, 
                                   localIndex ei );

  // Add the volumetric one-sided fluxes to the flux continuity constraint residual 
  void AssembleLocalContinuityConstraint( DomainPartition const * const domain,
                                          DofManager const * const dofManager,
                                          ParallelMatrix * const matrix,
                                          ParallelVector * const rhs );

};


} /* namespace geosx */

#endif //GEOSX_PHYSICSSOLVERS_FINITEVOLUME_SINGLEPHASEMIMETIC_HPP_
