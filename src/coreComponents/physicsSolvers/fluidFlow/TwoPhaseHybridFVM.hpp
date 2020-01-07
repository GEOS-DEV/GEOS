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

/**
 * @class TwoPhaseHybridFVM
 *
 * class to assemble the single-phase flow equations based
 * on a mixed hybrid formulation known as the mimetic method
 */
class TwoPhaseHybridFVM : public TwoPhaseBase
{
public:

  static constexpr localIndex MAX_NUM_FACES = 15;

  
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

    // weighted gravity for the transport gravity term
    static constexpr auto weightedGravityDepthString = "weightedGravityDepth";
    static constexpr auto sumTransmissibilityString  = "sumTransmissibility"; 
    
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


protected:

  void InitializePostInitialConditions_PreSubGroups( Group * const rootGroup ) override;
  
private:

  /**
   * @brief In a given element, compute the one-sided total volumetric fluxes at this element's faces
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] dFacePres the accumulated pressure updates at the mesh face 
   * @param[in] faceGravDepth the depth at the mesh faces
   * @param[in] elemToFaces elem-to-faces maps
   * @param[in] elemPres the pressure at this element's center
   * @param[in] dElemPres the accumulated pressure updates at this element's center
   * @param[in] elemGravDepth the depth at this element's center
   * @param[in] elemMob the phase mobilities at this elenent's center  
   * @param[in] dElemMob_dp the derivative of the mobilities at this element's center 
   * @param[in] dElemMob_dS the derivative of the mobilities at this element's center 
   * @param[in] elemDens the phase densities at this elenent's center  
   * @param[in] dElemDens_dp the derivative of the density at this element's center 
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[out] totalVolFlux the volumetric fluxes at this element's faces  
   * @param[out] dTotalVolFlux_dp the derivatives of the vol fluxes wrt to this element's elem centered pressure 
   * @param[out] dTotalVolFlux_dS the derivatives of the vol fluxes wrt to this element's elem centered saturation
   * @param[out] dTotalVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   *
   * For each face of the element, we compute the one sided total volumetric flux as: 
   * \sum_p T \lambda_p ( \nabla p - \rho_p g \nabla d)
   */
  void ComputeTotalVolFluxes( arrayView1d<real64 const> const & facePres,
                              arrayView1d<real64 const> const & dFacePres,
                              arrayView1d<real64 const> const & faceGravDepth,
                              arraySlice1d<localIndex const> const elemToFaces,
                              real64 const & elemPres,
                              real64 const & dElemPres,
                              real64 const & elemGravDepth,
                              arraySlice1d<real64 const> const elemMob,
                              arraySlice1d<real64 const> const dElemMob_dp,
                              arraySlice1d<real64 const> const dElemMob_dS,
                              arraySlice1d<real64 const> const elemDens,
                              arraySlice1d<real64 const> const dElemDens_dp,
                              stackArray2d<real64, MAX_NUM_FACES*MAX_NUM_FACES> const & transMatrix,
                              stackArray1d<real64, MAX_NUM_FACES> & totalVolFlux,
                              stackArray1d<real64, MAX_NUM_FACES> & dTotalVolFlux_dp,
                              stackArray1d<real64, MAX_NUM_FACES> & dTotalVolFlux_dS,
                              stackArray1d<real64, MAX_NUM_FACES> & dTotalVolFlux_dfp ) const;

  /**
   * @brief In a given element, compute the difference between phase gravity heads at this element's faces
   * @param[in] weightedFaceGravDepth the trans-weighted depth at the mesh faces
   * @param[in] elemToFaces elem-to-faces maps
   * @param[in] dens the densities in the domain (non-local)
   * @param[in] dDens_dp the derivatives of the densities in the domain wrt elem-centered pressure (non-local)
   * @param[in] elemIds the region, subregion and index of the local element
   * @param[in] neighborIds the region, subregion and index of the neighbor element
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[out] difGravHead the volumetric fluxes at this element's faces  
   * @param[out] dDifGravHead_dp the derivatives of the vol fluxes wrt to this element's centered pressure 
   *
   * For each face of the element, we compute the difference between phase gravity heads at this element's faces as: 
   * T ( \rho_p - \rho_m ) g \nabla d
   *
   * Note: because of the averaging of the densities across the face, this function requires non-local information
   */
  void ComputeGravityHead( arrayView1d<real64 const> const & weightedFaceGravDepth,
                           arraySlice1d<localIndex const> const elemToFaces,
                           real64 const & elemGravDepth,
                           ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dens,
                           ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dDens_dp,
                           stackArray1d<localIndex, 3>                       const & elemIds,
                           stackArray2d<localIndex, 3*MAX_NUM_FACES>         const & neighborIds,
                           stackArray2d<real64, MAX_NUM_FACES*MAX_NUM_FACES> const & transMatrix,
                           stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES>          & difGravHead,
                           stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES>        & dDifGravHead_dp ) const;

  /**
   * @brief In a given element, collect the upwinded mobility ratios at this element's faces 
   * @param[in] elemToFaces elem-to-faces maps
   * @param[in] mob the mobilities in the domain (non-local)
   * @param[in] dMob_dp the derivatives of the mobilities in the domain wrt elem-centered pressure (non-local)
   * @param[in] dMob_dS the derivatives of the mobilities in the domain wrt elem-centered saturation (non-local)
   * @param[in] dens the densities in the domain (non-local)
   * @param[in] dDens_dp the derivatives of the densities in the domain wrt elem-centered pressure (non-local)
   * @param[in] elemIds the region, subregion and index of the local element
   * @param[in] neighborIds the region, subregion and index of the neighbor element
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces  
   * @param[inout] viscousCoef the upwinded viscous mobility ratio at this element's faces  
   * @param[inout] dViscousCoef_dp the derivatives of the viscous mobilities ratio wrt the elem-centered pressures 
   * @param[inout] dViscousCoef_dS the derivatives of the viscous mobilities ratio wrt the elem-centered saturations 
   * @param[inout] gravCoef the upwinded buoyancy mobility ratio at this element's faces  
   * @param[inout] dGravCoef_dp the derivatives of the buoyancy mobilities ratio wrt the elem-centered pressures 
   * @param[inout] dGravCoef_dS the derivatives of the buoyancy mobilities ratio wrt the elem-centered saturations 
   *
   * Note: because of the upwinding, this function requires non-local information
   */
  void UpdateUpwindedCoefficients( arraySlice1d<localIndex const> const elemToFaces,
                                   ElementRegionManager::ElementViewAccessor<arrayView2d<real64>>  const & mob,
                                   ElementRegionManager::ElementViewAccessor<arrayView2d<real64>>  const & dMob_dp,
                                   ElementRegionManager::ElementViewAccessor<arrayView2d<real64>>  const & dMob_dS,
                                   ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dens,
                                   ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dDens_dp,
                                   stackArray1d<localIndex, 3>                    const & elemIds,
                                   stackArray2d<localIndex, 3*MAX_NUM_FACES>      const & neighborIds,
                                   stackArray1d<real64, MAX_NUM_FACES>            const & totalVolFlux,
                                   stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES> const & difGravHead,
                                   stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES>       & viscousCoef,
                                   stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES>     & dViscousCoef_dp,
                                   stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES>     & dViscousCoef_dS,
                                   stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES>       & gravCoef,
                                   stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES>     & dGravCoef_dp,
                                   stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES>     & dGravCoef_dS ) const;

  /**
   * @brief For a given one-sided face, collect the viscous upwinded mobility ratios
   * @param[in] mob the mobilities in the domain (non-local)
   * @param[in] dMob_dp the derivatives of the mobilities in the domain wrt elem-centered pressure (non-local)
   * @param[in] dMob_dS the derivatives of the mobilities in the domain wrt elem-centered saturation (non-local)
   * @param[in] dens the densities in the domain (non-local)
   * @param[in] dDens_dp the derivatives of the densities in the domain wrt elem-centered pressure (non-local)
   * @param[in] elemIds the region, subregion and index of the local element
   * @param[in] neighborIds the region, subregion and index of the neighbor element
   * @param[in] totalVolFlux total volumetric flux at the one-sided face 
   * @param[inout] viscousCoef the upwinded viscous mobility ratio at this face  
   * @param[inout] dViscousCoef_dp the derivatives of the viscous mobility ratios wrt the elem-centered pressures 
   * @param[inout] dViscousCoef_dS the derivatives of the viscous mobility ratios wrt the elem-centered saturations 
   *
   * We compute the viscous coefficient at this face as
   * \rho_p \lambda_p / \lambda_T
   * 
   */
  void UpdateLocalViscousCoefficients( ElementRegionManager::ElementViewAccessor<arrayView2d<real64>>  const & mob,
                                       ElementRegionManager::ElementViewAccessor<arrayView2d<real64>>  const & dMob_dp,
                                       ElementRegionManager::ElementViewAccessor<arrayView2d<real64>>  const & dMob_dS,
                                       ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dens,
                                       ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dDens_dp,
                                       stackArray1d<localIndex, 3>    const & elemIds,
                                       arraySlice1d<localIndex const> const neighborIds,
                                       real64 const & totalVolFlux,
                                       arraySlice1d<real64> const viscousCoef,
                                       arraySlice2d<real64> const dViscousCoef_dp,
                                       arraySlice2d<real64> const dViscousCoef_dS ) const;
  
  /**
   * @brief For a given one-sided face, collect the buoyancy upwinded mobility ratios
   * @param[in] mob the mobilities in the domain (non-local)
   * @param[in] dMob_dp the derivatives of the mobilities in the domain wrt elem-centered pressure (non-local)
   * @param[in] dMob_dS the derivatives of the mobilities in the domain wrt elem-centered saturation (non-local)
   * @param[in] dens the densities in the domain (non-local)
   * @param[in] dDens_dp the derivatives of the densities in the domain wrt elem-centered pressure (non-local)
   * @param[in] elemIds the region, subregion and index of the local element
   * @param[in] neighborIds the region, subregion and index of the neighbor element
   * @param[in] difGravCoef difference between the phase gravity heads T (\rho_p - \rho_m) g \nabla z
   * @param[inout] gravCoef the upwinded gravity mobility ratio at this face  
   * @param[inout] dGravCoef_dp the derivatives of the gravity mobility ratios wrt the elem-centered pressures 
   * @param[inout] dGravCoef_dS the derivatives of the gravity mobility ratios wrt the elem-centered saturations 
   *
   * We compute the viscous coefficient at this face as
   * \rho_p \lambda_p \lambda_m / \lambda_T
   * 
   */
  void UpdateLocalGravCoefficients( ElementRegionManager::ElementViewAccessor<arrayView2d<real64>>  const & mob,
                                    ElementRegionManager::ElementViewAccessor<arrayView2d<real64>>  const & dMob_dp,
                                    ElementRegionManager::ElementViewAccessor<arrayView2d<real64>>  const & dMob_dS,
                                    ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dens,
                                    ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dDens_dp,
                                    stackArray1d<localIndex, 3>    const & elemIds,
                                    arraySlice1d<localIndex const> const neighborIds,
                                    arraySlice1d<real64 const> const difGravCoef,                                                 
                                    arraySlice1d<real64> const gravCoef,
                                    arraySlice2d<real64> const dGravCoef_dp,
                                    arraySlice2d<real64> const dGravCoef_dS ) const;
  
  /**
   * @brief In a given element, assemble the mass conservation equations
   * @param[in] dt the time step size 
   * @param[in] faceDofNumber the dof numbers of the face pressures 
   * @param[in] elemToFaces elem-to-faces maps
   * @param[in] elemDofNumber the dof number of this element
   * @param[in] neighborDofNumbers the dof number of this element's neighbors
   * @param[in] totalVolFlux the total volumetric fluxes at this element's faces  
   * @param[in] dTotalVolFlux_dp the derivatives of the total vol fluxes wrt to this element's elem centered pressure
   * @param[in] dTotalVolFlux_dS the derivatives of the total vol fluxes wrt to this element's elem centered saturation
   * @param[in] dTotalVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[in] difGravHead the difference between the phase gravity heads at this element's faces  
   * @param[in] dDifGravHead_dp the derivatives of the difference between the phase gravity heads wrt to this element's elem centered pressure
   * @param[inout] viscousCoef the upwinded viscous coefficient (density times mobility ratio) at this element's faces  
   * @param[inout] dViscousCoef_dp the derivatives of the viscous coefficient wrt the elem-centered pressures 
   * @param[inout] dViscousCoef_dS the derivatives of the viscous coefficient wrt the elem-centered saturations 
   * @param[inout] gravCoef the upwinded gravity mobility ratio at this element's faces  
   * @param[inout] dGravCoef_dp the derivatives of the gravity mobility ratios wrt the elem-centered pressures 
   * @param[inout] dGravCoef_dS the derivatives of the gravity mobility ratios wrt the elem-centered saturations 
   * @param[inout] matrix the jacobian matrix
   * @param[inout] rhs the residual
   */
  void AssembleOneSidedMassFluxes( real64 const & dt,
                                   arrayView1d<globalIndex const> const & faceDofNumber,
                                   arraySlice1d<localIndex const> const elemToFaces,
                                   globalIndex const elemDofNumber,
                                   stackArray1d<globalIndex, MAX_NUM_FACES>         const & neighborDofNumbers,
                                   stackArray1d<real64, MAX_NUM_FACES>              const & totalVolFlux,
                                   stackArray1d<real64, MAX_NUM_FACES>              const & dTotalVolFlux_dp,
                                   stackArray1d<real64, MAX_NUM_FACES>              const & dTotalVolFlux_dS,
                                   stackArray1d<real64, MAX_NUM_FACES>              const & dTotalVolFlux_dfp,
                                   stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES>   const & difGravHead,
                                   stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES> const & dDifGravHead_dp,
                                   stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES>   const & viscousCoef,
                                   stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES> const & dViscousCoef_dp,
                                   stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES> const & dViscousCoef_dS,
                                   stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES>   const & gravCoef,
                                   stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES> const & dGravCoef_dp,
                                   stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES> const & dGravCoef_dS,
                                   ParallelMatrix * const matrix,
                                   ParallelVector * const rhs ) const;
  /**
   * @brief In a given element, assemble the constraints at this element's faces
   * @param[in] faceDofNumber the dof numbers of the face pressures 
   * @param[in] elemToFaces elem-to-faces maps
   * @param[in] elemDofNumber the dof number of this element's elem centered pressure
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces  
   * @param[in] dOneSidedVolFlux_dp the derivatives of the vol fluxes wrt to this element's elem centered pressure 
   * @param[in] dOneSidedVolFlux_dS the derivatives of the vol fluxes wrt to this element's elem centered saturation 
   * @param[in] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[inout] matrix the jacobian matrix
   * @param[inout] rhs the residual
   */
  void AssembleConstraints( arrayView1d<globalIndex const> const & faceDofNumber,
                            arraySlice1d<localIndex const> const elemToFaces,
                            globalIndex const elemDofNumber,
                            stackArray1d<real64, MAX_NUM_FACES> const & totalVolFlux,
                            stackArray1d<real64, MAX_NUM_FACES> const & dTotalVolFlux_dp,
                            stackArray1d<real64, MAX_NUM_FACES> const & dTotalVolFlux_dS,                        
                            stackArray1d<real64, MAX_NUM_FACES> const & dTotalVolFlux_dfp,
                            ParallelMatrix * const matrix,
                            ParallelVector * const rhs ) const; 


  /**
   * @brief For a given one-sided face, collect the buoyancy upwinded mobility ratios
   * @param[in] mesh the mesh object (single level only)
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] elemToFaces elem-to-faces maps
   * @param[in] ifaceLoc index of the one-sided face
   * @param[in] er index of the region of the local elem
   * @param[in] esr index of the subregion of the local elem
   * @param[in] ei index of the local elem
   * @param[inout] erNeighbor index of the region of the neighbor elem
   * @param[inout] esrNeighbor index of the subregion of the neighbor elem
   * @param[inout] eiNeighbor index of the neighbor elem
   * @param[inout] neighborDofNumber the dof number of the neighbor element pressure
   *
   */
  void FindAllNeighborsInTarget( MeshLevel const * const mesh,
                                 array2d<localIndex> const & elemRegionList,
                                 array2d<localIndex> const & elemSubRegionList,
                                 array2d<localIndex> const & elemList,
                                 set<localIndex>     const & regionFilter,
                                 arraySlice1d<localIndex const> const elemToFaces,
                                 stackArray1d<localIndex, 3> const & elemIds,
                                 globalIndex const elemDofNumber,   
                                 stackArray2d<localIndex, 3*MAX_NUM_FACES> & neighborIds,
                                 stackArray1d<globalIndex, MAX_NUM_FACES>  & neighborDofNumber ) const;


  /**
   * @brief This function generates various discretization information for later use.
   * @param domain the domain parition
   */
  void PrecomputeData( DomainPartition * const domain );
  
  /**
   * @brief In a given element, recompute the transmissibility matrix
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces elem-to-faces maps
   * @param[in] elemCenter the center of the element
   * @param[in] elemVolume the volume of the element
   * @param[in] elemPerm the permeability in the element
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
                                      stackArray2d<real64, MAX_NUM_FACES*MAX_NUM_FACES> & transMatrix ) const; 
  

  /**
   * @brief In a given element, recompute the transmissibility matrix using TPFA
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces elem-to-faces maps
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
                                stackArray2d<real64, MAX_NUM_FACES*MAX_NUM_FACES> const & transMatrix ) const;
  
  
  /// Dof key for the member functions that do not have access to the coupled Dof manager
  string m_faceDofKey; 

  /// Dof key for the member functions that do not have access to the coupled Dof manager
  string m_elemDofKey; 
  
  /// relative tolerance (redundant with FluxApproximationBase)
  real64 const m_areaRelTol;

  /// minimum value of the total mobility
  real64 const m_minTotalMob;
};

} /* namespace geosx */

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVM_HPP_
