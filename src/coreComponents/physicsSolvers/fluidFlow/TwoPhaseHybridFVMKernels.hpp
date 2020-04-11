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
 * @file TwoPhaseHybridFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_TWOPHASEHYBRIDFVMKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_TWOPHASEHYBRIDFVMKERNELS_HPP

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "physicsSolvers/fluidFlow/TwoPhaseBase.hpp"
#include "finiteVolume/HybridFVMInnerProduct.hpp"

namespace geosx
{

namespace TwoPhaseHybridFVMKernels
{

/******************************** FluxKernelHelper ********************************/

struct FluxKernelHelper
{

  /**
   * @brief In a given element, compute the one-sided volumetric fluxes at this element's faces
   * @param[in] facePotential the phase potential at the mesh faces at the beginning of the time step
   * @param[in] dFacePotential the accumulated phase potential updates at the mesh face
   * @param[in] faceGravCoef the depth at the mesh faces
   * @param[in] elemToFaces elem-to-faces maps
   * @param[in] elemPres the pressure at this element's center
   * @param[in] dElemPres the accumulated pressure updates at this element's center
   * @param[in] elemGravCoef the depth at this element's center
   * @param[in] elemDens the phase densities at this elenent's center
   * @param[in] dElemDens_dp the derivative of the density at this element's center
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[out] volFlux the volumetric fluxes at this element's faces
   * @param[out] dVolFlux_dp the derivatives of the vol fluxes wrt to this element's elem centered pressure
   * @param[out] dVolFlux_dS the derivatives of the vol fluxes wrt to this element's elem centered saturation
   * @param[out] dVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   *
   * For each face of the element, we compute the one sided volumetric flux as:
   * T ( \nabla p_\ell - \rho_p g \nabla d)
   */
  static
  void ComputeOneSidedVolFluxes( arrayView2d< real64 const > const & facePotential,
                                 arrayView2d< real64 const > const & dFacePotential,
                                 arrayView1d< real64 const > const & faceGravCoef,
                                 arraySlice1d< localIndex const > const elemToFaces,
                                 real64 const & elemPres,
                                 real64 const & dElemPres,
                                 real64 const & elemGravCoef,
                                 arraySlice1d< real64 const > const elemDens,
                                 arraySlice1d< real64 const > const dElemDens_dp,
                                 stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                               *HybridFVMInnerProduct::MAX_NUM_FACES > const & transMatrix,
                                 stackArray2d< real64, TwoPhaseBase::NUM_PHASES *HybridFVMInnerProduct::MAX_NUM_FACES > & volFlux,
                                 stackArray2d< real64, TwoPhaseBase::NUM_PHASES *HybridFVMInnerProduct::MAX_NUM_FACES > & dVlFlux_dp,
                                 stackArray2d< real64, TwoPhaseBase::NUM_PHASES *HybridFVMInnerProduct::MAX_NUM_FACES > & dVolFlux_dS,
                                 stackArray3d< real64,
                                               TwoPhaseBase::NUM_PHASES *HybridFVMInnerProduct::MAX_NUM_FACES *HybridFVMInnerProduct::MAX_NUM_FACES > & dTotalVolFlux_dfp );

  /**
   * @brief In a given element, collect the upwinded mobility ratios at this element's faces
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] elemToFaces elem-to-faces maps
   * @param[in] mob the mobilities in the domain (non-local)
   * @param[in] dMob_dp the derivatives of the mobilities in the domain wrt elem-centered pressure (non-local)
   * @param[in] dMob_dS the derivatives of the mobilities in the domain wrt elem-centered saturation (non-local)
   * @param[in] elemDofNumber the element dof numbers in the domain (non-local)
   * @param[in] elemIds the region, subregion and index of the local element
   * @param[in] volFlux the volumetric fluxes at this element's faces
   * @param[inout] upwMobility the upwinded phase mobilities at this element's faces
   * @param[inout] dUpwMobility_dp the derivatives of the upwinded phase mobilities wrt the elem-centered pressures
   * @param[inout] dUpwMobility_dS the derivatives of the upwinded phase mobilities wrt the elem-centered saturations
   * @param[inout] upwDofNumber the dof numbers of the upwind elements for each phase
   * @param[inout] neighborDofNumber the dof numbers of the neighbor elements
   *
   * Note: because of the upwinding, this function requires non-local information
   */
  static
  void UpdateUpwindedCoefficients( array2d< localIndex > const & elemRegionList,
                                   array2d< localIndex > const & elemSubRegionList,
                                   array2d< localIndex > const & elemList,
                                   SortedArray< localIndex > const & regionFilter,
                                   arraySlice1d< localIndex const > const elemToFaces,
                                   ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > >::ViewTypeConst const & mob,
                                   ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > >::ViewTypeConst const & dMob_dp,
                                   ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > >::ViewTypeConst const & dMob_dS,
                                   ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >::ViewTypeConst const & elemDofNumber,
                                   stackArray1d< localIndex, 3 > const & elemIds,
                                   stackArray2d< real64, TwoPhaseBase::NUM_PHASES *HybridFVMInnerProduct::MAX_NUM_FACES > const & volFlux,
                                   stackArray2d< real64, TwoPhaseBase::NUM_PHASES *HybridFVMInnerProduct::MAX_NUM_FACES > & upwMobility,
                                   stackArray2d< real64, TwoPhaseBase::NUM_PHASES *HybridFVMInnerProduct::MAX_NUM_FACES > & dUpwMobility_dp,
                                   stackArray2d< real64, TwoPhaseBase::NUM_PHASES *HybridFVMInnerProduct::MAX_NUM_FACES > & dUpwMobility_dS,
                                   stackArray2d< globalIndex, TwoPhaseBase::NUM_PHASES *HybridFVMInnerProduct::MAX_NUM_FACES > & upwDofNumber,
                                   stackArray1d< globalIndex, HybridFVMInnerProduct::MAX_NUM_FACES > & neighborDofNumber );

  /**
   * @brief In a given element, assemble the mass conservation equations
   * @param[in] dt the time step size
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] elemToFaces elem-to-faces maps
   * @param[in] elemDofNumber the dof number of this element
   * @param[in] phaseToRow map from phase id to row id
   * @param[in] volFlux the volumetric fluxes at this element's faces
   * @param[in] dVolFlux_dp the derivatives of the vol fluxes wrt to this element's elem centered pressure
   * @param[in] dVolFlux_dS the derivatives of the vol fluxes wrt to this element's elem centered saturation
   * @param[in] dVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[in] upwMobility the upwinded phase mobilities at this element's faces
   * @param[in] dUpwMobility_dp the derivatives of the upwinded phase mobilities wrt the elem-centered pressures
   * @param[in] dUpwMobility_dS the derivatives of the upwinded phase mobilities wrt the elem-centered saturations
   * @param[in] upwDofNumber the dof numbers of the upwind elements for each phase
   * @param[in] neighborDofNumber the dof numbers of the neighbor elements
   * @param[inout] matrix the jacobian matrix
   * @param[inout] rhs the residual
   */
  static
  void AssembleOneSidedMassFluxes( real64 const & dt,
                                   arrayView1d< globalIndex const > const & faceDofNumber,
                                   arraySlice1d< localIndex const > const elemToFaces,
                                   globalIndex const elemDofNumber,
                                   array1d< localIndex > const & phaseToRow,
                                   stackArray2d< real64, TwoPhaseBase::NUM_PHASES *HybridFVMInnerProduct::MAX_NUM_FACES > const & volFlux,
                                   stackArray2d< real64, TwoPhaseBase::NUM_PHASES *HybridFVMInnerProduct::MAX_NUM_FACES > const & dVolFlux_dp,
                                   stackArray2d< real64, TwoPhaseBase::NUM_PHASES *HybridFVMInnerProduct::MAX_NUM_FACES > const & dVolFlux_dS,
                                   stackArray3d< real64,
                                                 TwoPhaseBase::NUM_PHASES *HybridFVMInnerProduct::MAX_NUM_FACES *HybridFVMInnerProduct::MAX_NUM_FACES > const & dVolFlux_dfp,
                                   stackArray2d< real64, TwoPhaseBase::NUM_PHASES *HybridFVMInnerProduct::MAX_NUM_FACES > const & upwMobility,
                                   stackArray2d< real64, TwoPhaseBase::NUM_PHASES *HybridFVMInnerProduct::MAX_NUM_FACES > const & dUpwMobility_dp,
                                   stackArray2d< real64, TwoPhaseBase::NUM_PHASES *HybridFVMInnerProduct::MAX_NUM_FACES > const & dUpwMobility_dS,
                                   stackArray2d< globalIndex, TwoPhaseBase::NUM_PHASES *HybridFVMInnerProduct::MAX_NUM_FACES > const & upwDofNumber,
                                   stackArray1d< globalIndex, HybridFVMInnerProduct::MAX_NUM_FACES > const & neighborDofNumber,
                                   ParallelMatrix * const matrix,
                                   ParallelVector * const rhs );
  /**
   * @brief In a given element, assemble the constraints at this element's faces
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] elemToFaces elem-to-faces maps
   * @param[in] elemDofNumber the dof number of this element's elem centered pressure
   * @param[in] volFlux the volumetric fluxes at this element's faces
   * @param[in] dVolFlux_dp the derivatives of the vol fluxes wrt to this element's elem centered pressure
   * @param[in] dVolFlux_dS the derivatives of the vol fluxes wrt to this element's elem centered saturation
   * @param[in] dVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[inout] matrix the jacobian matrix
   * @param[inout] rhs the residual
   */
  static
  void AssembleConstraints( arrayView1d< globalIndex const > const & faceDofNumber,
                            arraySlice1d< localIndex const > const elemToFaces,
                            globalIndex const elemDofNumber,
                            stackArray2d< real64, TwoPhaseBase::NUM_PHASES *HybridFVMInnerProduct::MAX_NUM_FACES > const & volFlux,
                            stackArray2d< real64, TwoPhaseBase::NUM_PHASES *HybridFVMInnerProduct::MAX_NUM_FACES > const & dVolFlux_dp,
                            stackArray2d< real64, TwoPhaseBase::NUM_PHASES *HybridFVMInnerProduct::MAX_NUM_FACES > const & dVolFlux_dS,
                            stackArray3d< real64,
                                          TwoPhaseBase::NUM_PHASES *HybridFVMInnerProduct::MAX_NUM_FACES *HybridFVMInnerProduct::MAX_NUM_FACES > const & dVolFlux_dfp,
                            ParallelMatrix * const matrix,
                            ParallelVector * const rhs );

  /**
   * @brief For a given one-sided face, collect the buoyancy upwinded mobility ratios
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] elemToFaces elem-to-faces maps
   * @param[in] elemIds indices of the region, subregion, element
   * @param[in] elemDofNumber the element dof numbers in the domain (non-local)
   * @param[inout] neighborIds indices of the region, subregion, element of the neighbor
   * @param[inout] neighborDofNumber the dof number of the neighbor element pressure
   *
   */
  static
  void FindNeighborsInTarget( array2d< localIndex > const & elemRegionList,
                              array2d< localIndex > const & elemSubRegionList,
                              array2d< localIndex > const & elemList,
                              SortedArray< localIndex > const & regionFilter,
                              arraySlice1d< localIndex const > const elemToFaces,
                              stackArray1d< localIndex, 3 > const & elemIds,
                              ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >::ViewTypeConst const & elemDofNumber,
                              stackArray2d< localIndex, 3*HybridFVMInnerProduct::MAX_NUM_FACES > & neighborIds,
                              stackArray1d< globalIndex, HybridFVMInnerProduct::MAX_NUM_FACES > & neighborDofNumber );

};

/******************************** FluxKernel ********************************/

template< typename REGIONTYPE >
struct FluxKernel
{};

template<>
struct FluxKernel< CellElementSubRegion >
{
  /**
   * @brief In a given element, assemble the mass conservation equation and the contribution of this element to the face
   * constraints
   * @param[in] er index of this element's region
   * @param[in] esr index of this element's subregion
   * @param[in] ei index of this element
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] facePotential the phase potentials at the mesh faces at the beginning of the time step
   * @param[in] dFacePotential the accumulated potential updates at the mesh face
   * @param[in] faceGravDepth the depth at the mesh faces
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] elemPres the pressure at this element's center
   * @param[in] dElemPres the accumulated pressure updates at this element's center
   * @param[in] elemGravDepth the depth at this element's center
   * @param[in] fluidIndex the fluid index
   * @param[in] elemDens the density at this elenent's center
   * @param[in] dElemDens_dp the derivative of the density at this element's center
   * @param[in] mobility the mobilities in the domain (non-local)
   * @param[in] dMobility_dPres the derivatives of the mobilities in the domain wrt cell-centered pressure (non-local)
   * @param[in] dMobility_dSat the derivatives of the mobilities in the domain wrt cell-centered pressure (non-local)
   * @param[in] elemDofNumber the dof number of the cell centered pressures (non-local)
   * @param[in] phaseToRow map from phase id to row id
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[in] dt time step size
   * @param[inout] matrix the system matrix
   * @param[inout] rhs the system right-hand side vector
   */
  inline static void
  Compute( localIndex const er,
           localIndex const esr,
           localIndex const ei,
           SortedArray< localIndex > const & regionFilter,
           array2d< localIndex > const & elemRegionList,
           array2d< localIndex > const & elemSubRegionList,
           array2d< localIndex > const & elemList,
           arrayView1d< globalIndex const > const & faceDofNumber,
           arrayView2d< real64 const > const & facePotential,
           arrayView2d< real64 const > const & dFacePotential,
           arrayView1d< real64 const > const & faceGravCoef,
           arraySlice1d< localIndex const > const elemToFaces,
           real64 const & elemPres,
           real64 const & dElemPres,
           real64 const & elemGravCoef,
           localIndex const fluidIndex,
           ElementRegionManager::MaterialViewAccessor< arrayView3d< real64 const > >::ViewTypeConst const & elemPhaseDens,
           ElementRegionManager::MaterialViewAccessor< arrayView3d< real64 const > >::ViewTypeConst const & dElemPhaseDens_dp,
           ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > >::ViewTypeConst const & elemPhaseMob,
           ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > >::ViewTypeConst const & dElemPhaseMob_dp,
           ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > >::ViewTypeConst const & dElemPhaseMob_dS,
           ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >::ViewTypeConst const & elemDofNumber,
           array1d< localIndex > const & phaseToRow,
           stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                         *HybridFVMInnerProduct::MAX_NUM_FACES > const & transMatrix,
           real64 const & dt,
           ParallelMatrix * const matrix,
           ParallelVector * const rhs )
  {

    // max number of faces allowed in an element
    localIndex constexpr maxNumFaces = HybridFVMInnerProduct::MAX_NUM_FACES;
    // number of phases
    localIndex constexpr numPhases = TwoPhaseBase::NUM_PHASES;

    localIndex const numFacesInElem = elemToFaces.size();

    stackArray1d< globalIndex, maxNumFaces > neighborDofNumber( numFacesInElem );
    stackArray1d< localIndex, 3 > elemIds( 3 );
    elemIds[0] = er;
    elemIds[1] = esr;
    elemIds[2] = ei;
    neighborDofNumber = -1;

    // one sided phase flux: T ( \nabla p_\ell - \rho_\ell g \nabla d)
    stackArray2d< real64, numPhases *maxNumFaces > volFlux( numFacesInElem, numPhases );
    stackArray2d< real64, numPhases *maxNumFaces > dVolFlux_dp( numFacesInElem, numPhases );
    stackArray2d< real64, numPhases *maxNumFaces > dVolFlux_dS( numFacesInElem, numPhases );
    stackArray3d< real64, numPhases *maxNumFaces *maxNumFaces > dVolFlux_dfp( numFacesInElem, numFacesInElem, numPhases );

    // upwinded mobility
    stackArray2d< real64, numPhases *maxNumFaces > upwMobility( numFacesInElem, numPhases );
    stackArray2d< real64, numPhases *maxNumFaces > dUpwMobility_dp( numFacesInElem, numPhases );
    stackArray2d< real64, numPhases *maxNumFaces > dUpwMobility_dS( numFacesInElem, numPhases );
    stackArray2d< globalIndex, numPhases *maxNumFaces > upwDofNumber( numFacesInElem, numPhases );

    /*
     * In this function, we want to assemble the phase fluxes
     * To do that, we first compute at all the faces of this element
     * the one-sided phase "volumetric" fluxes
     * Once we know these quantities, we can upwind the mobilities
     */

    // For each one-sided face of the elem,
    // compute the phase volumetric fluxes using transMatrix
    FluxKernelHelper::ComputeOneSidedVolFluxes( facePotential,
                                                dFacePotential,
                                                faceGravCoef,
                                                elemToFaces,
                                                elemPres,
                                                dElemPres,
                                                elemGravCoef,
                                                elemPhaseDens[er][esr][fluidIndex][ei][0],
                                                dElemPhaseDens_dp[er][esr][fluidIndex][ei][0],
                                                transMatrix,
                                                volFlux,
                                                dVolFlux_dp,
                                                dVolFlux_dS,
                                                dVolFlux_dfp );


    // At this point, we know the local flow direction in the element
    // So we can upwind the transport coefficients (mobilities) at the one sided faces
    // ** this function needs non-local information because of the upwinding **
    FluxKernelHelper::UpdateUpwindedCoefficients( elemRegionList,
                                                  elemSubRegionList,
                                                  elemList,
                                                  regionFilter,
                                                  elemToFaces,
                                                  elemPhaseMob,
                                                  dElemPhaseMob_dp,
                                                  dElemPhaseMob_dS,
                                                  elemDofNumber,
                                                  elemIds,
                                                  volFlux,
                                                  upwMobility,
                                                  dUpwMobility_dp,
                                                  dUpwMobility_dS,
                                                  upwDofNumber,
                                                  neighborDofNumber );

    /*
     * At this point, we have computed all the fluxes in the element
     * We perform the assembly in this element in two steps:
     * 1) mass conservation equations
     * 2) face constraints
     */

    // Use the computed one sided phase vol fluxes and the upwinded mobilities
    // to assemble the upwinded mass fluxes in the mass conservation eqns of the elem
    FluxKernelHelper::AssembleOneSidedMassFluxes( dt,
                                                  faceDofNumber,
                                                  elemToFaces,
                                                  elemDofNumber[er][esr][ei],
                                                  phaseToRow,
                                                  volFlux,
                                                  dVolFlux_dp,
                                                  dVolFlux_dS,
                                                  dVolFlux_dfp,
                                                  upwMobility,
                                                  dUpwMobility_dp,
                                                  dUpwMobility_dS,
                                                  upwDofNumber,
                                                  neighborDofNumber,
                                                  matrix,
                                                  rhs );

    // Use the computed one sided vol fluxes to assemble the constraints
    // enforcing phase vol flux continuity at this element's faces
    FluxKernelHelper::AssembleConstraints( faceDofNumber,
                                           elemToFaces,
                                           elemDofNumber[er][esr][ei],
                                           volFlux,
                                           dVolFlux_dp,
                                           dVolFlux_dS,
                                           dVolFlux_dfp,
                                           matrix,
                                           rhs );
  }

};

} // namespace TwoPhaseHybridFVMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_TWOPHASEHYBRIDFVMKERNELS_HPP
