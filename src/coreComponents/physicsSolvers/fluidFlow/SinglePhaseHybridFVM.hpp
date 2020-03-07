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

  static constexpr localIndex MAX_NUM_FACES = 15;

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
  {
  } groupKeysSinglePhaseHybridFVM;

  groupKeyStruct & groupKeys()
  { return groupKeysSinglePhaseHybridFVM; }

  groupKeyStruct const & groupKeys() const
  { return groupKeysSinglePhaseHybridFVM; }


private:

  /**
   * @brief In a given element, compute the one-sided volumetric fluxes at this element's faces
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] dFacePres the accumulated pressure updates at the mesh face 
   * @param[in] faceGravDepth the depth at the mesh faces
   * @param[in] elemToFaces the map from one-sided face to face 
   * @param[in] elemPres the pressure at this element's center
   * @param[in] dElemPres the accumulated pressure updates at this element's center
   * @param[in] elemGravDepth the depth at this element's center
   * @param[in] elemDens the density at this elenent's center  
   * @param[in] dElemDens_dp the derivative of the density at this element's center 
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[out] oneSidedVolFlux the volumetric fluxes at this element's faces  
   * @param[out] dOneSidedVolFlux_dp the derivatives of the vol fluxes wrt to this element's cell centered pressure 
   * @param[out] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   */
  void ComputeOneSidedVolFluxes( arrayView1d<real64 const> const & facePres,
                                 arrayView1d<real64 const> const & dFacePres,
                                 arrayView1d<real64 const> const & faceGravDepth,
                                 arraySlice1d<localIndex const> const elemToFaces,
                                 real64 const & elemPres,
                                 real64 const & dElemPres,
                                 real64 const & elemGravDepth,
                                 real64 const & elemDens,
                                 real64 const & dElemDens_dp,
                                 stackArray2d<real64, MAX_NUM_FACES
                                                     *MAX_NUM_FACES> const & transMatrix,
                                 stackArray1d<real64, MAX_NUM_FACES> & oneSidedVolFlux,
                                 stackArray1d<real64, MAX_NUM_FACES> & dOneSidedVolFlux_dp,
                                 stackArray1d<real64, MAX_NUM_FACES> & dOneSidedVolFlux_dfp ) const;

  /**
   * @brief In a given element, collect the upwinded mobilities at this element's faces 
   * @param[in] mesh the mesh object (single level only)
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] elemToFaces elem-to-faces maps
   * @param[in] mob the mobilities in the domain (non-local)
   * @param[in] dMob_dp the derivatives of the mobilities in the domain wrt cell-centered pressure (non-local)
   * @param[in] er index of this element's region
   * @param[in] esr index of this element's subregion  
   * @param[in] ei index of this element 
   * @param[in] elemDofNumber the dof number of this element's cell centered pressure
   * @param[in] elemDofKey 
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces  
   * @param[inout] upwMobility the upwinded mobilities at this element's faces  
   * @param[inout] dUpwMobility_dp the derivatives of the upwinded mobilities wrt the cell-centered pressures (local or neighbor)  
   * @param[inout] upwDofNumber  the dof number of the upwind pressure 
   *
   * Note: because of the upwinding, this function requires non-local information
   */
  void UpdateUpwindedCoefficients( MeshLevel const * const mesh,
                                   array2d<localIndex> const & elemRegionList,
                                   array2d<localIndex> const & elemSubRegionList,
                                   array2d<localIndex> const & elemList,
                                   SortedArray<localIndex> const & regionFilter,
                                   arraySlice1d<localIndex const> const elemToFaces,
                                   ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & mob,
                                   ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & dMob_dp,
                                   localIndex const er,
                                   localIndex const esr,
                                   localIndex const ei,
                                   globalIndex const elemDofNumber,
                                   string const elemDofKey,    
                                   stackArray1d<real64, MAX_NUM_FACES> const & oneSidedVolFlux,
                                   stackArray1d<real64, MAX_NUM_FACES> & upwMobility,
                                   stackArray1d<real64, MAX_NUM_FACES> & dUpwMobility_dp,
                                   stackArray1d<globalIndex, MAX_NUM_FACES> & upwDofNumber ) const;

  /**
   * @brief In a given element, assemble the mass conservation equation
   * @param[in] dt the time step size 
   * @param[in] faceDofNumber the dof numbers of the face pressures 
   * @param[in] elemToFaces the map from one-sided face to face to access face Dof numbers
   * @param[in] elemDofNumber the dof number of this element's cell centered pressure
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces  
   * @param[in] dOneSidedVolFlux_dp the derivatives of the vol fluxes wrt to this element's cell centered pressure 
   * @param[in] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[in] upwMobility the upwinded mobilities at this element's faces  
   * @param[in] dUpwMobility_dp the derivatives of the upwinded mobilities wrt the cell-centered pressures (local or neighbor)  
   * @param[in] upwDofNumber  the dof number of the upwind pressure 
   * @param[inout] matrix the jacobian matrix
   * @param[inout] rhs the residual
   */
  void AssembleOneSidedMassFluxes( real64 const & dt,
                                   arrayView1d<globalIndex const> const & faceDofNumber,
                                   arraySlice1d<localIndex const> const elemToFaces,
                                   globalIndex const elemDofNumber,
                                   stackArray1d<real64, MAX_NUM_FACES> const & oneSidedVolFlux,
                                   stackArray1d<real64, MAX_NUM_FACES> const & dOneSidedVolFlux_dp,
                                   stackArray1d<real64, MAX_NUM_FACES> const & dOneSidedVolFlux_dfp,
                                   stackArray1d<real64, MAX_NUM_FACES> const & upwMobility,
                                   stackArray1d<real64, MAX_NUM_FACES> const & dUpwMobility_dp,
                                   stackArray1d<globalIndex, MAX_NUM_FACES> const & upwDofNumber,
                                   ParallelMatrix * const matrix,
                                   ParallelVector * const rhs ) const;


  /**
   * @brief In a given element, assemble the constraints at this element's faces
   * @param[in] faceDofNumber the dof numbers of the face pressures 
   * @param[in] elemToFaces the map from one-sided face to face to access face Dof numbers
   * @param[in] elemDofNumber the dof number of this element's cell centered pressure
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces  
   * @param[in] dOneSidedVolFlux_dp the derivatives of the vol fluxes wrt to this element's cell centered pressure 
   * @param[in] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[inout] matrix the jacobian matrix
   * @param[inout] rhs the residual
   */
  void AssembleConstraints( arrayView1d<globalIndex const> const & faceDofNumber,
                            arraySlice1d<localIndex const> const elemToFaces,
                            globalIndex const elemDofNumber,
                            stackArray1d<real64, MAX_NUM_FACES> const & oneSidedVolFlux,
                            stackArray1d<real64, MAX_NUM_FACES> const & dOneSidedVolFlux_dp,
                            stackArray1d<real64, MAX_NUM_FACES> const & dOneSidedVolFlux_dfp,
                            ParallelMatrix * const matrix,
                            ParallelVector * const rhs ) const; 


  /**
   * @brief In a given element, recompute the transmissibility matrix
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemVolume the volume of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] transMatrix
   *
   * This function is in this class until we find a better place for it
   * 
   */
  void ComputeTransmissibilityMatrix( arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & nodePosition, 
                                      ArrayOfArraysView<localIndex const> const & faceToNodes, 
                                      arraySlice1d<localIndex const> const elemToFaces,
                                      R1Tensor const & elemCenter,
                                      real64   const & elemVolume,
                                      R1Tensor const & elemPerm,
                                      real64   const & lengthTolerance,
                                      stackArray2d<real64, MAX_NUM_FACES
                                                          *MAX_NUM_FACES> const & transMatrix ) const; 

  
  /**
   * @brief In a given element, recompute the transmissibility matrix using TPFA
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] transMatrix
   *
   * This function is in this class until we find a better place for it
   * 
   */
  void ComputeTPFAInnerProduct( arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & nodePosition, 
                                ArrayOfArraysView<localIndex const> const & faceToNodes, 
                                arraySlice1d<localIndex const> const elemToFaces,
                                R1Tensor const & elemCenter,
                                R1Tensor const & elemPerm,
                                real64   const & lengthTolerance,
                                stackArray2d<real64, MAX_NUM_FACES
                                                    *MAX_NUM_FACES> const & transMatrix ) const; 

  /**
   * @brief In a given element, recompute the transmissibility matrix using a consistent inner product
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] tParam parameter used in the transmissibility matrix computations 
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] transMatrix
   *   
   * When tParam = 2, we obtain a scheme that reduces to TPFA
   * on orthogonal meshes, but remains consistent on non-orthogonal meshes
   *
   * This function is in this class until we find a better place for it
   * 
   */
  void ComputeQFamilyInnerProduct( arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & nodePosition, 
                                   ArrayOfArraysView<localIndex const> const & faceToNodes, 
                                   arraySlice1d<localIndex const> const elemToFaces,
                                   R1Tensor const & elemCenter,
                                   real64   const & elemVolume,
                                   R1Tensor const & elemPerm,
                                   real64   const & tParam, 
                                   real64   const & lengthTolerance,
                                   stackArray2d<real64, MAX_NUM_FACES
                                                       *MAX_NUM_FACES> const & transMatrix ) const; 


  
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
