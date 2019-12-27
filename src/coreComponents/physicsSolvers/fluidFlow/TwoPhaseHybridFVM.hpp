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

namespace geosx
{

  static constexpr localIndex MAX_NUM_FACES_IN_ELEM = 15;
  
/**
 * @class TwoPhaseHybridFVM
 *
 * class to assemble the single-phase flow equations based
 * on a mixed hybrid formulation known as the mimetic method
 */
class TwoPhaseHybridFVM : public TwoPhaseBase
{
public:

  static constexpr localIndex MAX_NUM_FACES_IN_ELEM = 15;

  
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


  struct viewKeyStruct : TwoPhaseBase::viewKeyStruct
  {
    static constexpr auto faceDofFieldString = "faceCenteredVariables";
    
    // primary face-based field
    static constexpr auto facePressureString      = "facePressure";
    static constexpr auto deltaFacePressureString = "deltaFacePressure";
  
  } viewKeysTwoPhaseHybridFVM;

  viewKeyStruct & viewKeys()
  { return viewKeysTwoPhaseHybridFVM; }

  viewKeyStruct const & viewKeys() const
  { return viewKeysTwoPhaseHybridFVM; }

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysTwoPhaseHybridFVM;

  groupKeyStruct & groupKeys()
  { return groupKeysTwoPhaseHybridFVM; }

  groupKeyStruct const & groupKeys() const
  { return groupKeysTwoPhaseHybridFVM; }


private:

  /**
   * @brief In a given element, compute the one-sided volumetric fluxes at this element's faces
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] dFacePres the accumulated pressure updates at the mesh face 
   * @param[in] faceGravDepth the depth at the mesh faces
   * @param[in] oneSidedFaceToFace the map from one-sided face to face 
   * @param[in] elemPres the pressure at this element's center
   * @param[in] dElemPres the accumulated pressure updates at this element's center
   * @param[in] elemGravDepth the depth at this element's center
   * @param[in] elemPhaseDens the density at this elenent's center  
   * @param[in] dElemPhaseDens_dp the derivative of the density at this element's center 
   * @param[in] elemOffset the offset of this element in the map oneSidedFaceToFace
   * @param[in] numFacesInElem the number of faces in this element
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[out] oneSidedVolFlux the volumetric fluxes at this element's faces  
   * @param[out] dOneSidedVolFlux_dp the derivatives of the vol fluxes wrt to this element's cell centered pressure 
   * @param[out] dOneSidedVolFlux_dS the derivatives of the vol fluxes wrt to this element's cell centered saturation
   * @param[out] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   */
  void ComputeOneSidedVolFluxes( arrayView1d<real64 const> const & facePres,
                                 arrayView1d<real64 const> const & dFacePres,
                                 arrayView1d<real64 const> const & faceGravDepth,
                                 arraySlice1d<localIndex const> const elemToFaces,
                                 real64 const & elemPres,
                                 real64 const & dElemPres,
                                 real64 const & elemGravDepth,
                                 arraySlice1d<real64 const> const elemPhaseMob,
                                 arraySlice1d<real64 const> const dElemPhaseMob_dp,
                                 arraySlice1d<real64 const> const dElemPhaseMob_dS,
                                 arraySlice1d<real64 const> const elemPhaseDens,
                                 arraySlice1d<real64 const> const dElemPhaseDens_dp,
                                 stackArray2d<real64, MAX_NUM_FACES_IN_ELEM
                                                     *MAX_NUM_FACES_IN_ELEM> const & transMatrix,
                                 stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> & oneSidedVolFlux,
                                 stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> & dOneSidedVolFlux_dp,
                                 stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> & dOneSidedVolFlux_dS,
                                 stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> & dOneSidedVolFlux_dfp ) const;

  /**
   * @brief In a given element, collect the upwinded mobilities at this element's faces 
   * @param[in] neighborRegionId the region index of the neighbor element of the one-sided faces (non-local)
   * @param[in] neighborSubRegionId the subregion index of the neighbor element of the one-sided faces (non-local)
   * @param[in] neighborElemId the elem index of the neighbor element of the one-sided faces (non-local)
   * @param[in] neighborDofNumber the dof number of the neighbor element of the one-sided faces (non-local)
   * @param[in] domainMobility the mobilities in the domain (non-local)
   * @param[in] dDomainMobility_dp the derivatives of the mobilities in the domain wrt cell-centered pressure (non-local)
   * @param[in] elemMobility the mobility in this element
   * @param[in] dElemMobility_dp the derivative of the mobility wrt this element cell-centered pressure 
   * @param[in] dElemMobility_dS the derivative of the mobility wrt this element cell-centered saturation 
   * @param[in] elemDofNumber the dof number of this element's cell centered pressure
   * @param[in] elemOffset the offset of this element in the map oneSidedFaceToFace
   * @param[in] numFacesInElem the number of faces in this element
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces  
   * @param[inout] upwMobility the upwinded mobilities at this element's faces  
   * @param[inout] dUpwMobility_dp the derivatives of the upwinded mobilities wrt the cell-centered pressures (local or neighbor)  
   * @param[inout] dUpwMobility_dS the derivatives of the upwinded mobilities wrt the cell-centered saturations (local or neighbor)  
   * @param[inout] upwDofNumber  the dof number of the upwind pressure 
   *
   * Note: because of the upwinding, this function requires non-local information
   */
   void UpdateUpwindedCoefficients( MeshLevel const * const  mesh,
                                    array2d<localIndex> const & elemRegionList,
                                    array2d<localIndex> const & elemSubRegionList,
                                    array2d<localIndex> const & elemList,
                                    set<localIndex> const & regionFilter,
                                    arraySlice1d<localIndex const> const elemToFaces,
                                    ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & domainMobility,
                                    ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & dDomainMobility_dp,
                                    ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & dDomainMobility_dS,
                                    localIndex const er,
                                    localIndex const esr,
                                    localIndex const ei,
                                    globalIndex const elemDofNumber,
                                    string const elemDofKey,
                                    stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & oneSidedVolFlux,
                                    stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> & viscousMobilityRatio,
                                    stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> & dViscousMobilityRatio_dp,
                                    stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> & dViscousMobilityRatio_dS,
                                    stackArray1d<globalIndex, MAX_NUM_FACES_IN_ELEM> & viscousDofNumber,
                                    stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> & buoyancyMobilityRatio,
                                    stackArray2d<real64, 2*MAX_NUM_FACES_IN_ELEM> & dBuoyancyMobilityRatio_dp,
                                    stackArray2d<real64, 2*MAX_NUM_FACES_IN_ELEM> & dBuoyancyMobilityRatio_dS,
                                    stackArray2d<globalIndex, 2*MAX_NUM_FACES_IN_ELEM> & buoyancyDofNumber ) const;

  /**
   * @brief In a given element, assemble the mass conservation equation
   * @param[in] dt the time step size 
   * @param[in] faceDofNumber the dof numbers of the face pressures 
   * @param[in] oneSidedFaceToFace the map from one-sided face to face to access face Dof numbers
   * @param[in] elemDofNumber the dof number of this element's cell centered pressure
   * @param[in] elemOffset the offset of this element in the map oneSidedFaceToFace
   * @param[in] numFacesInElem the number of faces in this element
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces  
   * @param[in] dOneSidedVolFlux_dp the derivatives of the vol fluxes wrt to this element's cell centered pressure 
   * @param[in] dOneSidedVolFlux_dS the derivatives of the vol fluxes wrt to this element's cell centered saturation
   * @param[in] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[in] upwPhaseMobility the upwinded mobilities at this element's faces  
   * @param[in] dUpwPhaseMobility_dp the derivatives of the upwinded mobilities wrt the cell-centered pressures (local or neighbor)  
   * @param[in] dUpwPhaseMobility_dS the derivatives of the upwinded mobilities wrt the cell-centered saturations (local or neighbor)  
   * @param[in] upwDofNumber  the dof number of the upwind pressure 
   * @param[inout] matrix the jacobian matrix
   * @param[inout] rhs the residual
   */
  void AssembleOneSidedMassFluxes( real64 const & dt,
                                   arrayView1d<globalIndex const> const & faceDofNumber,
                                   arraySlice1d<localIndex const> const elemToFaces,
                                   globalIndex const elemDofNumber,
                                   stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & oneSidedVolFlux,
                                   stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dOneSidedVolFlux_dp,
                                   stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dOneSidedVolFlux_dS,
                                   stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dOneSidedVolFlux_dfp,
                                   stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & viscousMobilityRatio,
                                   stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dViscousMobilityRatio_dp,
                                   stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dViscousMobilityRatio_dS,
                                   stackArray1d<globalIndex, MAX_NUM_FACES_IN_ELEM> const & viscousDofNumber,
                                   stackArray1d<real64, MAX_NUM_FACES_IN_ELEM>   const & buoyancyMobilityRatio,
                                   stackArray2d<real64, 2*MAX_NUM_FACES_IN_ELEM> const & dBuoyancyMobilityRatio_dp,
                                   stackArray2d<real64, 2*MAX_NUM_FACES_IN_ELEM> const & dBuoyancyMobilityRatio_dS,
                                   stackArray2d<globalIndex, 2*MAX_NUM_FACES_IN_ELEM> const & buoyancyDofNumber,
                                   ParallelMatrix * const matrix,
                                   ParallelVector * const rhs ) const;

  /**
   * @brief In a given element, assemble the constraints at this element's faces
   * @param[in] faceDofNumber the dof numbers of the face pressures 
   * @param[in] oneSidedFaceToFace the map from one-sided face to face to access face Dof numbers
   * @param[in] elemDofNumber the dof number of this element's cell centered pressure
   * @param[in] elemOffset the offset of this element in the map oneSidedFaceToFace
   * @param[in] numFacesInElem the number of faces in this element
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces  
   * @param[in] dOneSidedVolFlux_dp the derivatives of the vol fluxes wrt to this element's cell centered pressure 
   * @param[in] dOneSidedVolFlux_dS the derivatives of the vol fluxes wrt to this element's cell centered saturation
   * @param[in] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[inout] matrix the jacobian matrix
   * @param[inout] rhs the residual
   */
  void AssembleConstraints( arrayView1d<globalIndex const> const & faceDofNumber,
                            arraySlice1d<localIndex const> const elemToFaces,
                            globalIndex const elemDofNumber,
                            stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & oneSidedVolFlux,
                            stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dOneSidedVolFlux_dp,
                            stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dOneSidedVolFlux_dS,                        
                            stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dOneSidedVolFlux_dfp,
                            ParallelMatrix * const matrix,
                            ParallelVector * const rhs ) const; 


  /**
   * @brief In a given element, recompute the transmissibility matrix
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] oneSidedFaceToFace the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] permeability the permeability in the element
   * @param[in] elemOffset the position of the element in one-sided face based maps
   * @param[in] numFacesInElem the number of faces in this element
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] transMatrix
   *
   * This function is in this class until we find a better place for it
   * 
   */
  void ComputeTransmissibilityMatrix( arrayView1d<R1Tensor const> const & nodePosition, 
                                      ArrayOfArraysView<localIndex const> const & faceToNodes, 
                                      arraySlice1d<localIndex const> const elemToFaces,
                                      R1Tensor const & elemCenter,
                                      real64   const & elemVolume,
                                      R1Tensor const & elemPerm,
                                      real64   const & lengthTolerance,
                                      stackArray2d<real64, MAX_NUM_FACES_IN_ELEM
                                                          *MAX_NUM_FACES_IN_ELEM> & transMatrix ) const; 
  

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
  void ComputeTPFAInnerProduct( arrayView1d<R1Tensor const> const & nodePosition, 
                                ArrayOfArraysView<localIndex const> const & faceToNodes, 
                                arraySlice1d<localIndex const> const elemToFaces,
                                R1Tensor const & elemCenter,
                                R1Tensor const & elemPerm,
                                real64   const & lengthTolerance,
                                stackArray2d<real64, MAX_NUM_FACES_IN_ELEM
                                                    *MAX_NUM_FACES_IN_ELEM> const & transMatrix ) const;
  
  
  /// Dof key for the member functions that do not have access to the coupled Dof manager
  string m_faceDofKey; 
  
  /// relative tolerance (redundant with FluxApproximationBase)
  real64 m_areaRelTol;
  
};

} /* namespace geosx */

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVM_HPP_
