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

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEMIMETIC_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEMIMETIC_HPP_

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
 * @class SinglePhaseMimetic
 *
 * class to assemble the single-phase flow equations based
 * on a mixed hybrid formulation known as the mimetic method
 */
class SinglePhaseMimetic : public SinglePhaseFlowBase
{
public:

  static constexpr localIndex MAX_NUM_FACES_IN_ELEM = 15;

  struct EqType
  {
    static constexpr integer MASS_CONS  = 0;
    static constexpr integer CONSTRAINT = 1;
  };
  
  
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
                     ParallelVector * const rhs ) override;

  /**@}*/


  struct viewKeyStruct : SinglePhaseFlowBase::viewKeyStruct
  {
    // primary face-based field
    static constexpr auto facePressureString      = "facePressure";
    static constexpr auto deltaFacePressureString = "deltaFacePressure";
  
    // one-sided face-based connectivity maps
    static constexpr auto oneSidedFaceToFaceString  = "oneSidedFaceToFace"; 
    static constexpr auto neighborRegionIdString    = "neighborRegionIndex";
    static constexpr auto neighborSubRegionIdString = "neighborSubRegionIndex";
    static constexpr auto neighborElemIdString      = "neighborElemIndex";
    static constexpr auto neighborDofNumberString   = "neighborDofNumber";

    // elem-based map to access the one-sided face vars from the elements
    static constexpr auto elemOffsetString         = "elemOffsetString";

    
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


private:

    /**
   * @brief Register the one-sided face maps on the mesh entities
   * @param[in,out] MeshBody the group of MeshBody objects to register data on
   */
  void RegisterOneSidedFaceData( Group * const MeshBodies );


  /**
   * @brief Set the size of the one-sided face based arrays
   * @param[in,out] domain the domain containing the mesh and fields
   *
   * This function is used to
   *  1. Count the number of one-sided faces in the mesh and set m_numOneSidedFaces
   *  2. Resize the one-sided face based arrays using m_numOneSidedFaces 
   */
  void ResizeOneSidedFaceFields( DomainPartition * const domain );

  
  /**
   * @brief Construct the connectivity maps that are used in the mimetic method
   * @param[in,out] domain the domain containing the mesh and fields
   *
   * In this class, we assemble the fluxes element-by-element. In each
   * element, we iterate over the one-sided faces of the element. This function
   * constructs the arrays storing the connectivity information necessary to 
   * conveniently do that.
   */
  void ConstructOneSidedFaceMaps( DomainPartition * const domain,
                                  DofManager const & dofManager );

  /**
   * @brief In a given element, compute the one-sided volumetric fluxes at this element's faces
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] dFacePres the accumulated pressure updates at the mesh face 
   * @param[in] faceGravDepth the depth at the mesh faces
   * @param[in] oneSidedFaceToFace the map from one-sided face to face 
   * @param[in] elemPres the pressure at this element's center
   * @param[in] dElemPres the accumulated pressure updates at this element's center
   * @param[in] elemGravDepth the depth at this element's center
   * @param[in] elemDens the density at this elenent's center  
   * @param[in] dElemDens_dp the derivative of the density at this element's center 
   * @param[in] elemOffset the offset of this element in the map oneSidedFaceToFace
   * @param[in] numFacesInElem the number of faces in this element
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[out] oneSidedVolFlux the volumetric fluxes at this element's faces  
   * @param[out] dOneSidedVolFlux_dp the derivatives of the vol fluxes wrt to this element's cell centered pressure 
   * @param[out] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   */
  void ComputeOneSidedVolFluxes( arrayView1d<real64 const> const & facePres,
                                 arrayView1d<real64 const> const & dFacePres,
                                 arrayView1d<real64 const> const & faceGravDepth,
                                 arrayView1d<localIndex const> const & oneSidedFaceToFace,
                                 real64 const & elemPres,
                                 real64 const & dElemPres,
                                 real64 const & elemGravDepth,
                                 real64 const & elemDens,
                                 real64 const & dElemDens_dp,
                                 localIndex const elemOffset,
                                 localIndex const numFacesInElem,
                                 stackArray2d<real64, MAX_NUM_FACES_IN_ELEM
                                                     *MAX_NUM_FACES_IN_ELEM> const & transMatrix,
                                 stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> & oneSidedVolFlux,
                                 stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> & dOneSidedVolFlux_dp,
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
   * @param[in] elemDofNumber the dof number of this element's cell centered pressure
   * @param[in] elemOffset the offset of this element in the map oneSidedFaceToFace
   * @param[in] numFacesInElem the number of faces in this element
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces  
   * @param[inout] upwMobility the upwinded mobilities at this element's faces  
   * @param[inout] dUpwMobility_dp the derivatives of the upwinded mobilities wrt the cell-centered pressures (local or neighbor)  
   * @param[inout] upwDofNumber  the dof number of the upwind pressure 
   *
   * Note: because of the upwinding, this function requires non-local information
   */
  void UpdateUpwindedCoefficients( arrayView1d<localIndex const> const & neighborRegionId,
                                   arrayView1d<localIndex const> const & neighborSubRegionId,
                                   arrayView1d<localIndex const> const & neighborElemId,
                                   arrayView1d<globalIndex const> const & neighborDofNumber,
                                   ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & domainMobility,
                                   ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & dDomainMobility_dp,
                                   real64 const & elemMobility,
                                   real64 const & dElemMobility_dp,
                                   globalIndex const elemDofNumber,
                                   localIndex const elemOffset,
                                   localIndex const numFacesInElem,
                                   stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & oneSidedVolFlux,
                                   stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> & upwMobility,
                                   stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> & dUpwMobility_dp,
                                   stackArray1d<globalIndex, MAX_NUM_FACES_IN_ELEM> & upwDofNumber ) const;

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
   * @param[in] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[in] upwMobility the upwinded mobilities at this element's faces  
   * @param[in] dUpwMobility_dp the derivatives of the upwinded mobilities wrt the cell-centered pressures (local or neighbor)  
   * @param[in] upwDofNumber  the dof number of the upwind pressure 
   * @param[inout] matrix the jacobian matrix
   * @param[inout] rhs the residual
   */
  void AssembleOneSidedMassFluxes( real64 const & dt,
                                   arrayView1d<globalIndex const> const & faceDofNumber,
                                   arrayView1d<localIndex const> const & oneSidedFaceToFace,
                                   globalIndex const elemDofNumber,
                                   localIndex const elemOffset,
                                   localIndex const numFacesInElem,
                                   stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & oneSidedVolFlux,
                                   stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dOneSidedVolFlux_dp,
                                   stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dOneSidedVolFlux_dfp,
                                   stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & upwMobility,
                                   stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dUpwMobility_dp,
                                   stackArray1d<globalIndex, MAX_NUM_FACES_IN_ELEM> const & upwDofNumber,
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
   * @param[in] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[inout] matrix the jacobian matrix
   * @param[inout] rhs the residual
   */
  void AssembleConstraints( arrayView1d<globalIndex const> const & faceDofNumber,
                            arrayView1d<localIndex const> const & oneSidedFaceToFace,
                            globalIndex const elemDofNumber,
                            localIndex const elemOffset,
                            localIndex const numFacesInElem,
                            stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & oneSidedVolFlux,
                            stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dOneSidedVolFlux_dp,
                            stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dOneSidedVolFlux_dfp,
                            ParallelMatrix * const matrix,
                            ParallelVector * const rhs ) const;


  /**
   * @brief In a given element, recompute the transmissibility matrix
   * @param[in] X the position of the nodes
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
                                      arrayView1d<localIndex const> const & oneSidedFaceToFace, 
                                      R1Tensor const & elemCenter, 
                                      R1Tensor const & elemPerm,
                                      real64 const   & elemOffset,
                                      real64 const   & numFacesInElem,
                                      real64 const   & lengthTolerance,
                                      stackArray2d<real64, MAX_NUM_FACES_IN_ELEM
                                                          *MAX_NUM_FACES_IN_ELEM> & transMatrix ) const; 
  

  /// Number of one-sided faces on this MPI rank
  localIndex m_numOneSidedFaces;

  /// relative tolerance (redundant with FluxApproximationBase)
  real64 m_areaRelTol;
  
};

} /* namespace geosx */

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEMIMETIC_HPP_
