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

  struct CellPos
  {
    static constexpr integer LOCAL    = 0;
    static constexpr integer NEIGHBOR = 1;
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

    // auxiliary face-based maps
    static constexpr auto numAdjacentElementsString = "numAdjacentElements";
    static constexpr auto faceToOneSidedFaceString  = "faceToOneSidedFace";
    static constexpr auto faceToElemDofNumberString = "faceToElemDofNumberString";

    // auxiliary one-sided face-based vars
    static constexpr auto oneSidedVolFluxString                = "oneSidedVolumetricFlux";
    static constexpr auto dOneSidedVolFlux_dPressureString     = "dOneSidedVolumetricFlux_dPressure";
    static constexpr auto dOneSidedVolFlux_dFacePressureString = "dOneSidedVolumetricFlux_dFacePressure";

    static constexpr auto upwMobilityString            = "upwindedMobility";
    static constexpr auto dUpwMobility_dPressureString = "dUpwindedMobility_dPressure";
    
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

  /**
   * @brief Setup stored views into domain data for the current step
   */
  void ResetViews( DomainPartition * const domain ) override;

protected:

  virtual void InitializePostInitialConditions_PreSubGroups( dataRepository::Group * const rootGroup ) override;

private:

    /**
   * @brief Register the one-sided face maps on the mesh entities
   * @param[in,out] MeshBody the group of MeshBody objects to register data on
   */
  void RegisterOneSidedFaceData( Group * const MeshBodies );

  
  /**
   * @brief Set the size of the arrays needed to implement the solver
   * @param[in,out] domain the domain containing the mesh and fields
   */
  void ResizeFaceFields( DomainPartition * const domain );
  
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
   * @brief Reset the face maps before calling ConstructFaceMaps  
   * @param[in,out] domain the domain containing the mesh and fields
   */
  void ResetFaceMaps( DomainPartition * const domain );
  
  /**
   * @brief Construct the connectivity maps that are used in the mimetic method
   * @param[in,out] domain the domain containing the mesh and fields
   *
   * In the mimetic method, we assemble the fluxes element-by-element. In each
   * element, we iterate over the one-sided faces of the element. This function
   * constructs the arrays storing the connectivity information necessary to 
   * conveniently do that.
   */
  void ConstructFaceMaps( DomainPartition * const domain );

  /**
   * @brief Compute the one-sided volumetric fluxes for all the elements of the domain
   * @param[in] domain the domain containing the mesh and fields
   *
   * For each element, this function loops over the one-sided faces. For each
   * one-sided face, it computes and stores the volumetric flux.  
   */
  void ComputeOneSidedVolFluxes( DomainPartition const * const domain );

  /**
   * @brief Collect the upwind mobility for all the one-sided faces of the domain
   * @param[in] domain the domain containing the mesh and fields
   *
   * For each element, this function loops over the one-sided faces. For each
   * one-sided face, it uses the sign of the one-sided volumetric flux to 
   * detect the upwind element. The upwind mobilities are then stored
   */
  void UpdateUpwindedTransportCoefficients( DomainPartition const * const domain );

  /**
   * @brief Assemble the upwinded one-sided mass fluxes for all the elements of the domain
   * @param time_n time at the beginning of the step
   * @param dt the perscribed timestep
   * @param domain the domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   * 
   * For each element, this function loops over the one-sided faces. For 
   * each one-sided face, it multiplies the upwinded mobility by the one-sided
   * volumetric flux to obtain the upwinded one-sided mass flux. Then this  
   * function performs the assembly of the flux terms in the residual and Jacobian
   * matrix
   */
  void AssembleUpwindedOneSidedMassFluxes( real64 const time_n,
                                           real64 const dt,
                                           DomainPartition const * const domain,
                                           DofManager const * const dofManager,
                                           ParallelMatrix * const matrix,
                                           ParallelVector * const rhs );

  /**
   * @brief Assemble the constraints enforcing flux continuity at the faces
   * @param time_n time at the beginning of the step
   * @param dt the prescribed timestep
   * @param domain the domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   * 
   * For each face, we have to make sure that the one-sided fluxes coming from the 
   * two adjacent elements match. This is done by imposing a constraint at each face.
   * To assemble these constraints, we loop over the elements. For each element, we 
   * loop over the one-sided faces. For each one-sided face, we add to the contribution
   * of the one-sided flux to the constraint at this face.
   */
  void AssembleConstraints( real64 const time_n,
                            real64 const dt,
                            DomainPartition const * const domain,
                            DofManager const * const dofManager,
                            ParallelMatrix * const matrix,
                            ParallelVector * const rhs );


  /**
   * @brief At a given one-sided face, compute the mass flux
   * @param dt time step
   * @param upwMobility 
   * @param dUpwMobility_dp 
   * @param dUpwMobility_dp_neighbor
   * @param oneSidedVolFlux
   * @param dOneSidedVolFlux_dp
   * @param dOneSidedVolFlux_dfp
   * @param sumOneSidedMassFluxes
   * @param dSumOneSidedMassFluxes_dp
   * @param dSumOneSidedMassFluxes_dp_neighbor
   * @param dSumOneSidedMassFluxes_dfp
   */
  void IncrementLocalMassFluxSum( real64 const & dt,
                                  real64 const & upwMobility,
                                  real64 const & dUpwMobility_dp,
                                  real64 const & dUpwMobility_dp_neighbor,
                                  real64 const & oneSidedVolFlux,
                                  real64 const & dOneSidedVolFlux_dp,
                                  real64 const & dOneSidedVolFlux_dfp,
                                  real64       & sumOneSidedMassFluxes,
                                  real64       & dSumOneSidedMassFluxes_dp,
                                  real64       & dSumOneSidedMassFluxes_dp_neighbor,
                                  real64       & dSumOneSidedMassFluxes_dfp ) const;


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
   * @param[out] oneSidedTrans 
   *
   * This function is in this class until we find a better place for it
   * 
   */
  void RecomputeOneSidedTransmissibilities( arrayView1d<R1Tensor const> const & X, 
 	      	      		            ArrayOfArraysView<localIndex const> const & faceToNodes, 
					    arrayView1d<localIndex const> const & oneSidedFaceToFace, 
					    R1Tensor const & elemCenter, 
 					    R1Tensor const & permeability,
					    real64 const   & elemOffset,
					    real64 const   & numFacesInElem,
					    real64 const   & lengthTolerance,
					    stackArray2d<real64, MAX_NUM_FACES_IN_ELEM
 					                        *MAX_NUM_FACES_IN_ELEM> & oneSidedTrans ) const; 
  

  /// Number of one-sided faces on this MPI rank
  localIndex m_numOneSidedFaces;

  /// relative tolerance (redundant with FluxApproximationBase)
  real64 m_areaRelTol;
  
};

} /* namespace geosx */

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEMIMETIC_HPP_
