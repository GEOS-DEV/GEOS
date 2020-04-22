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
 * @file TwoPhaseHybridFVM.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_TWOPHASEHYBRIDFVM_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_TWOPHASEHYBRIDFVM_HPP_

#include "physicsSolvers/fluidFlow/TwoPhaseBase.hpp"
#include "finiteVolume/HybridFVMInnerProduct.hpp"

namespace geosx
{

/**
 * @class TwoPhaseHybridFVM
 *
 * class to assemble the single-phase flow equations based
 * on a mixed hybrid finite-volume formulation
 */
class TwoPhaseHybridFVM : public TwoPhaseBase
{
public:

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
  TwoPhaseHybridFVM( const std::string & name,
                     Group * const parent );


  /// deleted default constructor
  TwoPhaseHybridFVM() = delete;

  /// deleted copy constructor
  TwoPhaseHybridFVM( TwoPhaseHybridFVM const & ) = delete;

  /// default move constructor
  TwoPhaseHybridFVM( TwoPhaseHybridFVM && ) = default;

  /// deleted assignment operator
  TwoPhaseHybridFVM & operator=( TwoPhaseHybridFVM const & ) = delete;

  /// deleted move operator
  TwoPhaseHybridFVM & operator=( TwoPhaseHybridFVM && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~TwoPhaseHybridFVM() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName()
  { return "TwoPhaseHybridFVM"; }

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
   * @brief assembles the flux terms for all elems
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


  struct viewKeyStruct : TwoPhaseBase::viewKeyStruct
  {
    static constexpr auto faceDofFieldString = "faceCenteredVariables";

    // primary face-based field
    static constexpr auto facePhasePotentialString      = "facePhasePotential";
    static constexpr auto deltaFacePhasePotentialString = "deltaFacePhasePotential";

  } viewKeysTwoPhaseHybridFVM;

  viewKeyStruct & viewKeys()
  { return viewKeysTwoPhaseHybridFVM; }

  viewKeyStruct const & viewKeys() const
  { return viewKeysTwoPhaseHybridFVM; }

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {} groupKeysTwoPhaseHybridFVM;

  groupKeyStruct & groupKeys()
  { return groupKeysTwoPhaseHybridFVM; }

  groupKeyStruct const & groupKeys() const
  { return groupKeysTwoPhaseHybridFVM; }

protected:

  virtual void ResizeFields( MeshLevel & meshLevel ) override;

private:

  /**
   * @brief Assemble the mass conservation equations and face constraints in the face subregion
   * @param[in] er index of this element's region
   * @param[in] esr index of this element's subregion
   * @param[in] subRegion reference to the face element subregion
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] mesh the mesh object (single level only)
   * @param[in] nodePosition position of the nodes
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] faceToNodes map from face to nodes
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] facePotential the pressure at the mesh faces at the beginning of the time step
   * @param[in] dFacePotential the accumulated pressure updates at the mesh face
   * @param[in] faceGravCoef the depth at the mesh faces
   * @param[in] lengthTolerance tolerance used in the transmissibility calculations
   * @param[in] dt time step size
   * @param[in] dofManager the dof manager
   * @param[inout] matrix the system matrix
   * @param[inout] rhs the system right-hand side vector
   */
  void FluxLaunch( localIndex er,
                   localIndex esr,
                   FaceElementSubRegion const & subRegion,
                   SortedArray< localIndex > regionFilter,
                   MeshLevel const & mesh,
                   arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                   array2d< localIndex > const & elemRegionList,
                   array2d< localIndex > const & elemSubRegionList,
                   array2d< localIndex > const & elemList,
                   ArrayOfArraysView< localIndex const > const & faceToNodes,
                   arrayView1d< globalIndex const > const & faceDofNumber,
                   arrayView2d< real64 const > const & facePotential,
                   arrayView2d< real64 const > const & dFacePotential,
                   arrayView1d< real64 const > const & faceGravCoef,
                   real64 const lengthTolerance,
                   real64 const dt,
                   DofManager const * const dofManager,
                   ParallelMatrix * const matrix,
                   ParallelVector * const rhs );

  /**
   * @brief Assemble the mass conservation equations and face constraints in the face subregion
   * @param[in] er index of this element's region
   * @param[in] esr index of this element's subregion
   * @param[in] subRegion reference to the cell element subregion
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] mesh the mesh object (single level only)
   * @param[in] nodePosition position of the nodes
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] faceToNodes map from face to nodes
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] facePotential the pressure at the mesh faces at the beginning of the time step
   * @param[in] dFacePotential the accumulated pressure updates at the mesh face
   * @param[in] faceGravCoef the depth at the mesh faces
   * @param[in] lengthTolerance tolerance used in the transmissibility calculations
   * @param[in] dt time step size
   * @param[in] dofManager the dof manager
   * @param[inout] matrix the system matrix
   * @param[inout] rhs the system right-hand side vector
   */
  void FluxLaunch( localIndex er,
                   localIndex esr,
                   CellElementSubRegion const & subRegion,
                   SortedArray< localIndex > regionFilter,
                   MeshLevel const & mesh,
                   arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                   array2d< localIndex > const & elemRegionList,
                   array2d< localIndex > const & elemSubRegionList,
                   array2d< localIndex > const & elemList,
                   ArrayOfArraysView< localIndex const > const & faceToNodes,
                   arrayView1d< globalIndex const > const & faceDofNumber,
                   arrayView2d< real64 const > const & facePotential,
                   arrayView2d< real64 const > const & dFacePotential,
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
   * @param[in] elemToFaces elem-to-faces maps
   * @param[in] elemCenter the center of the element
   * @param[in] elemVolume the volume of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[in] orthonormalizeWithSVD flag to indicate whether SVD is used (if not, Gram-Schmidt is used)
   * @param[inout] transMatrix
   *
   * I will probably move this function somewhere else (to HybridFVMInnerProduct) at some point
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
                                      stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES *HybridFVMInnerProduct::MAX_NUM_FACES > & transMatrix ) const;


  /// Dof key for the member functions that do not have access to the coupled Dof manager
  string m_faceDofKey;

  /// Dof key for the member functions that do not have access to the coupled Dof manager
  string m_elemDofKey;

  /// relative tolerance (redundant with FluxApproximationBase)
  real64 const m_areaRelTol;

  /// type of inner product for the mimetic method
  /// This is const (hard-coded) for now
  integer const m_ipType;

  /// flag to decide we orthonormalize with SVD or with MGS
  /// This is const (hard-coded) for now
  bool const m_orthonormalizeWithSVD;

};

} /* namespace geosx */

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_TWOPHASEHYBRIDFVM_HPP_
