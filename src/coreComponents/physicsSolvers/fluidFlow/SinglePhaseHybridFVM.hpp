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
 * @file SinglePhaseHybridFVM.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVM_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVM_HPP_

#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseHybridFVMKernels.hpp"

namespace geosx
{


/**
 * @class SinglePhaseHybridFVM
 *
 * class to assemble the single-phase flow equations based
 * on a mixed hybrid formulation known as the mimetic method
 */
class SinglePhaseHybridFVM : public SinglePhaseBase
{
public:

  struct EquationType
  {
    static constexpr integer MASS_CONS  = 0;
    static constexpr integer CONSTRAINT = 1;
  };

  struct InnerProductType
  {
    static constexpr integer TPFA = 0;
    static constexpr integer QUASI_TPFA = 1;
  };


  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  SinglePhaseHybridFVM( const std::string & name,
                        Group * const parent );


  /// deleted default constructor
  SinglePhaseHybridFVM() = delete;

  /// deleted copy constructor
  SinglePhaseHybridFVM( SinglePhaseHybridFVM const & ) = delete;

  /// default move constructor
  SinglePhaseHybridFVM( SinglePhaseHybridFVM && ) = default;

  /// deleted assignment operator
  SinglePhaseHybridFVM & operator=( SinglePhaseHybridFVM const & ) = delete;

  /// deleted move operator
  SinglePhaseHybridFVM & operator=( SinglePhaseHybridFVM && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseHybridFVM() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName()
  { return "SinglePhaseHybridFVM"; }

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
                     ParallelVector * const rhs ) override;


  /**@}*/


  struct viewKeyStruct : SinglePhaseBase::viewKeyStruct
  {
    // primary face-based field
    static constexpr auto facePressureString      = "facePressure";
    static constexpr auto deltaFacePressureString = "deltaFacePressure";

  } viewKeysSinglePhaseHybridFVM;

  viewKeyStruct & viewKeys()
  { return viewKeysSinglePhaseHybridFVM; }

  viewKeyStruct const & viewKeys() const
  { return viewKeysSinglePhaseHybridFVM; }

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {} groupKeysSinglePhaseHybridFVM;

  groupKeyStruct & groupKeys()
  { return groupKeysSinglePhaseHybridFVM; }

  groupKeyStruct const & groupKeys() const
  { return groupKeysSinglePhaseHybridFVM; }


private:

  /**
   * @brief Assemble the mass conservation equations and face constraints in the face subregion
   * @param[in] er index of this element's region
   * @param[in] esr index of this element's subregion
   * @param[in] subRegion pointer to the cell element subregion
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] mesh the mesh object (single level only)
   * @param[in] nodePosition position of the nodes
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] faceToNodes map from face to nodes
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] dFacePres the accumulated pressure updates at the mesh face
   * @param[in] faceGravCoef the depth at the mesh faces
   * @param[in] lengthTolerance tolerance used in the transmissibility calculations
   * @param[in] dt time step size
   * @param[in] dofManager the dof manager
   * @param[inout] matrix the system matrix
   * @param[inout] rhs the system right-hand side vector
   */
  void FluxLaunch( localIndex er,
                   localIndex esr,
                   FaceElementSubRegion const * const subRegion,
                   SortedArrayView< localIndex const > const & regionFilter,
                   MeshLevel const * const mesh,
                   arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                   arrayView2d< localIndex const > const & elemRegionList,
                   arrayView2d< localIndex const > const & elemSubRegionList,
                   arrayView2d< localIndex const > const & elemList,
                   ArrayOfArraysView< localIndex const > const & faceToNodes,
                   arrayView1d< globalIndex const > const & faceDofNumber,
                   arrayView1d< real64 const > const & facePres,
                   arrayView1d< real64 const > const & dFacePres,
                   arrayView1d< real64 const > const & faceGravCoef,
                   real64 const lengthTolerance,
                   real64 const dt,
                   DofManager const * const dofManager,
                   ParallelMatrix * const matrix,
                   ParallelVector * const rhs );

  /**
   * @brief Assemble the mass conservation equations and face constraints in the cell subregion
   * @param[in] er index of this element's region
   * @param[in] esr index of this element's subregion
   * @param[in] subRegion pointer to the cell element subregion
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] mesh the mesh object (single level only)
   * @param[in] nodePosition position of the nodes
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] faceToNodes map from face to nodes
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] dFacePres the accumulated pressure updates at the mesh face
   * @param[in] faceGravCoef the depth at the mesh faces
   * @param[in] lengthTolerance tolerance used in the transmissibility calculations
   * @param[in] dt time step size
   * @param[in] dofManager the dof manager
   * @param[inout] matrix the system matrix
   * @param[inout] rhs the system right-hand side vector
   */
  void FluxLaunch( localIndex er,
                   localIndex esr,
                   CellElementSubRegion const * const subRegion,
                   SortedArrayView< localIndex const > const & regionFilter,
                   MeshLevel const * const mesh,
                   arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                   arrayView2d< localIndex const > const & elemRegionList,
                   arrayView2d< localIndex const > const & elemSubRegionList,
                   arrayView2d< localIndex const > const & elemList,
                   ArrayOfArraysView< localIndex const > const & faceToNodes,
                   arrayView1d< globalIndex const > const & faceDofNumber,
                   arrayView1d< real64 const > const & facePres,
                   arrayView1d< real64 const > const & dFacePres,
                   arrayView1d< real64 const > const & faceGravCoef,
                   real64 const lengthTolerance,
                   real64 const dt,
                   DofManager const * const dofManager,
                   ParallelMatrix * const matrix,
                   ParallelVector * const rhs );

  /**
   * @brief In a given element, recompute the transmissibility matrix
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemVolume the volume of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] orthonormalizeWithSVD flag to indicate whether SVD is used (if not, Gram-Schmidt is used)
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] transMatrix
   *
   * This function is in this class until we find a better place for it
   *
   */
  void ComputeTransmissibilityMatrix( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                                      ArrayOfArraysView< localIndex const > const & faceToNodes,
                                      arraySlice1d< localIndex const > const elemToFaces,
                                      R1Tensor const & elemCenter,
                                      real64 const & elemVolume,
                                      R1Tensor const & elemPerm,
                                      real64 const & lengthTolerance,
                                      bool const & orthonormalizeWithSVD,
                                      stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                    *HybridFVMInnerProduct::MAX_NUM_FACES > const & transMatrix ) const;


  /// Dof key for the member functions that do not have access to the coupled Dof manager
  string m_faceDofKey;

  /// relative tolerance (redundant with FluxApproximationBase)
  real64 m_areaRelTol;

  /// type of inner product for the mimetic method
  /// This is only const for now
  integer const m_ipType;

  /// flag to decide we orthonormalize with SVD or with MGS
  /// This is only const for now
  bool const m_orthonormalizeWithSVD;

};

} /* namespace geosx */

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVM_HPP_
