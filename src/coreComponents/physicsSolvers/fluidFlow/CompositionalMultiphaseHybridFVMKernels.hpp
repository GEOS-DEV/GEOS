/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompositionalMultiphaseHybridFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMKERNELS_HPP

#include "common/DataTypes.hpp"
#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"
#include "finiteVolume/mimeticInnerProducts/TPFAInnerProduct.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mesh/MeshLevel.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVMAssemblerHelperKernels.hpp"

namespace geosx
{

namespace CompositionalMultiphaseHybridFVMKernels
{

/******************************** AssemblerKernel ********************************/

struct AssemblerKernel
{

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  /**
   * @brief In a given element, assemble the mass conservation equations and the contribution of this element to the face
   * constraints
   * @param[in] er index of this element's region
   * @param[in] esr index of this element's subregion
   * @param[in] ei index of this element
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] faceGhostRank ghost rank of each face
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] dFacePres the accumulated pressure updates at the mesh face
   * @param[in] faceGravCoef the depth at the mesh faces
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] elemPres the pressure at this element's center
   * @param[in] dElemPres the accumulated pressure updates at this element's center
   * @param[in] elemGravCoef the depth at this element's center
   * @param[in] phaseDens the phase densities in the domain (non-local)
   * @param[in] dPhaseDens_dPres the derivatives of the phase densities in the domain wrt pressure (non-local)
   * @param[in] dPhaseDens_dCompFrac the derivatives of the phase densities in the domain wrt component fraction (non-local)
   * @param[in] phaseMob the phase mobilities in the domain (non-local)
   * @param[in] dPhaseMob_dPres the derivatives of the phase mobilities in the domain wrt pressure (non-local)
   * @param[in] dPhaseMob_dCompDens the derivatives of the phase mobilities in the domain wrt component fraction (non-local)
   * @param[in] dCompFrac_dCompDens the derivatives of the component fraction in the domain wrt component density (non-local)
   * @param[in] phaseCompFrac the phase component fractions in the domain (non-local)
   * @param[in] dPhaseCompFrac_dPres the derivatives of the phase component fractions in the domain wrt pressure (non-local)
   * @param[in] dPhaseCompFrac_dCompFrac the derivatives of the phase component fractions in the domain wrt component fraction (non-local)
   * @param[in] elemDofNumber the dof number of the cell centered pressures (non-local)
   * @param[in] elemGhostRank the ghost rank of the element in which we assemble the fluxes
   * @param[in] rankOffset the offset of this rank
   * @param[in] lengthTolerance tolerance used in the transmissibility matrix computation
   * @param[in] dt time step size
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[inout] matrix the system matrix
   * @param[inout] rhs the system right-hand side vector
   */
  template< localIndex NF, localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
  compute( localIndex const er, localIndex const esr, localIndex const ei,
           SortedArrayView< localIndex const > const & regionFilter,
           arrayView2d< localIndex const > const & elemRegionList,
           arrayView2d< localIndex const > const & elemSubRegionList,
           arrayView2d< localIndex const > const & elemList,
           arrayView1d< globalIndex const > const & faceDofNumber,
           arrayView1d< integer const > const & faceGhostRank,
           arrayView1d< real64 const > const & facePres,
           arrayView1d< real64 const > const & dFacePres,
           arrayView1d< real64 const > const & faceGravCoef,
           arrayView1d< real64 const > const & mimFaceGravCoef,
           arraySlice1d< localIndex const > const & elemToFaces,
           real64 const & elemPres,
           real64 const & dElemPres,
           real64 const & elemGravCoef,
           ElementViewConst< arrayView3d< real64 const > > const & phaseDens,
           ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres,
           ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dCompFrac,
           ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
           ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
           ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dCompFrac,
           ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
           ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
           ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dCompDens,
           ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
           ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac,
           ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
           ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dCompFrac,
           ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
           integer const elemGhostRank,
           globalIndex const rankOffset,
           real64 const & dt,
           arraySlice2d< real64 const > const & transMatrix,
           arraySlice2d< real64 const > const & transMatrixGrav,
           CRSMatrixView< real64, globalIndex const > const & localMatrix,
           arrayView1d< real64 > const & localRhs )
  {
    // one sided flux
    real64 oneSidedVolFlux[ NF ]{};
    real64 dOneSidedVolFlux_dPres[ NF ]{};
    real64 dOneSidedVolFlux_dFacePres[ NF ][ NF ]{};
    real64 dOneSidedVolFlux_dCompDens[ NF ][ NC ]{};

    localIndex const localIds[3] = { er, esr, ei };

    /*
     * compute auxiliary quantities at the one sided faces of this element:
     * 1) One-sided volumetric fluxes
     * 2) Upwinded mobilities
     */

    // for each one-sided face of the elem,
    // compute the volumetric flux using transMatrix
    AssemblerKernelHelper::applyGradient< NF, NC, NP >( facePres,
                                                        dFacePres,
                                                        faceGravCoef,
                                                        elemToFaces,
                                                        elemPres,
                                                        dElemPres,
                                                        elemGravCoef,
                                                        phaseMassDens[er][esr][ei][0],
                                                        dPhaseMassDens_dPres[er][esr][ei][0],
                                                        dPhaseMassDens_dCompFrac[er][esr][ei][0],
                                                        phaseMob[er][esr][ei],
                                                        dPhaseMob_dPres[er][esr][ei],
                                                        dPhaseMob_dCompDens[er][esr][ei],
                                                        dCompFrac_dCompDens[er][esr][ei],
                                                        transMatrix,
                                                        oneSidedVolFlux,
                                                        dOneSidedVolFlux_dPres,
                                                        dOneSidedVolFlux_dFacePres,
                                                        dOneSidedVolFlux_dCompDens );

    // at this point, we know the local flow direction in the element
    // so we can upwind the transport coefficients (mobilities) at the one sided faces
    // ** this function needs non-local information **
    if( elemGhostRank < 0 )
    {
      /*
       * perform assembly in this element in two steps:
       * 1) mass conservation equations
       * 2) face constraints
       */

      // use the computed one sided vol fluxes and the upwinded mobilities
      // to assemble the upwinded mass fluxes in the mass conservation eqn of the elem
      AssemblerKernelHelper::assembleFluxDivergence< NF, NC, NP >( localIds,
                                                                   rankOffset,
                                                                   elemRegionList,
                                                                   elemSubRegionList,
                                                                   elemList,
                                                                   regionFilter,
                                                                   faceDofNumber,
                                                                   mimFaceGravCoef,
                                                                   elemToFaces,
                                                                   elemGravCoef,
                                                                   phaseDens,
                                                                   dPhaseDens_dPres,
                                                                   dPhaseDens_dCompFrac,
                                                                   phaseMassDens,
                                                                   dPhaseMassDens_dPres,
                                                                   dPhaseMassDens_dCompFrac,
                                                                   phaseMob,
                                                                   dPhaseMob_dPres,
                                                                   dPhaseMob_dCompDens,
                                                                   dCompFrac_dCompDens,
                                                                   phaseCompFrac,
                                                                   dPhaseCompFrac_dPres,
                                                                   dPhaseCompFrac_dCompFrac,
                                                                   elemDofNumber,
                                                                   transMatrixGrav,
                                                                   oneSidedVolFlux,
                                                                   dOneSidedVolFlux_dPres,
                                                                   dOneSidedVolFlux_dFacePres,
                                                                   dOneSidedVolFlux_dCompDens,
                                                                   dt,
                                                                   localMatrix,
                                                                   localRhs );
    }

    // use the computed one sided vol fluxes to assemble the constraints
    // enforcing flux continuity at this element's faces
    AssemblerKernelHelper::assembleFaceConstraints< NF, NC, NP >( faceDofNumber,
                                                                  faceGhostRank,
                                                                  elemToFaces,
                                                                  elemDofNumber[er][esr][ei],
                                                                  rankOffset,
                                                                  oneSidedVolFlux,
                                                                  dOneSidedVolFlux_dPres,
                                                                  dOneSidedVolFlux_dFacePres,
                                                                  dOneSidedVolFlux_dCompDens,
                                                                  localMatrix,
                                                                  localRhs );

  }


};


/******************************** FluxKernel ********************************/

struct FluxKernel
{

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  /**
   * @brief In a given subRegion, assemble the mass conservation equations and the contribution of the elements of this subRegion  to the
   * face constraints.
   * @param[in] er index of this element's region
   * @param[in] esr index of this element's subregion
   * @param[in] subRegion the subRegion in which we are going to assemble the fluxes
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] nodePosition position of the nodes
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] faceToNodes map from face to nodes
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] faceGhostRank  the ghost ranks of the face pressures
   * @param[in] facePres the pressures at the mesh faces at the beginning of the time step
   * @param[in] dFacePres the accumulated pressure updates at the mesh face
   * @param[in] faceGravDepth the depth at the mesh faces
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] elemPres the pressure at this element's center
   * @param[in] dElemPres the accumulated pressure updates at this element's center
   * @param[in] elemGravDepth the depth at this element's center
   * @param[in] phaseDens the phase densities in the domain (non-local)
   * @param[in] dPhaseDens_dPres the derivatives of the phase densities in the domain wrt pressure (non-local)
   * @param[in] dPhaseDens_dCompFrac the derivatives of the phase densities in the domain wrt component fraction (non-local)
   * @param[in] phaseMob the phase mobilities in the domain (non-local)
   * @param[in] dPhaseMob_dPres the derivatives of the phase mobilities in the domain wrt pressure (non-local)
   * @param[in] dPhaseMob_dCompDens the derivatives of the phase mobilities in the domain wrt component fraction (non-local)
   * @param[in] dCompFrac_dCompDens the derivatives of the component fraction in the domain wrt component density (non-local)
   * @param[in] phaseCompFrac the phase component fractions in the domain (non-local)
   * @param[in] dPhaseCompFrac_dPres the derivatives of the phase component fractions in the domain wrt pressure (non-local)
   * @param[in] dPhaseCompFrac_dCompFrac the derivatives of the phase component fractions in the domain wrt component fraction (non-local)
   * @param[in] elemDofNumber the dof number of the cell centered pressures (non-local)
   * @param[in] rankOffset the offset of this rank
   * @param[in] lengthTolerance tolerance used in the transmissibility matrix computation
   * @param[in] dt time step size
   * @param[inout] matrix the system matrix
   * @param[inout] rhs the system right-hand side vector
   */
  template< typename IP_TYPE, localIndex NF, localIndex NC, localIndex NP >
  static void
  launch( localIndex er, localIndex esr,
          CellElementSubRegion const & subRegion,
          SortedArrayView< localIndex const > const & regionFilter,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
          arrayView2d< localIndex const > const & elemRegionList,
          arrayView2d< localIndex const > const & elemSubRegionList,
          arrayView2d< localIndex const > const & elemList,
          ArrayOfArraysView< localIndex const > const & faceToNodes,
          arrayView1d< globalIndex const > const & faceDofNumber,
          arrayView1d< integer const > const & faceGhostRank,
          arrayView1d< real64 const > const & facePres,
          arrayView1d< real64 const > const & dFacePres,
          arrayView1d< real64 const > const & faceGravCoef,
          arrayView1d< real64 const > const & mimFaceGravCoef,
          arrayView1d< real64 const > const & transMultiplier,
          ElementViewConst< arrayView3d< real64 const > > const & phaseDens,
          ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres,
          ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dCompFrac,
          ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
          ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
          ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dCompFrac,
          ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
          ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dCompDens,
          ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
          ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac,
          ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
          ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dCompFrac,
          ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
          globalIndex const rankOffset,
          real64 const lengthTolerance,
          real64 const dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
  {
    // get the cell-centered DOF numbers and ghost rank for the assembly
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

    // get the map from elem to faces
    arrayView2d< localIndex const > const & elemToFaces = subRegion.faceList();

    // get the cell-centered pressures
    arrayView1d< real64 const > const & elemPres  =
      subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseBase::viewKeyStruct::pressureString );
    arrayView1d< real64 const > const & dElemPres =
      subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseBase::viewKeyStruct::deltaPressureString );

    // get the element data needed for transmissibility computation
    arrayView2d< real64 const > const & elemCenter =
      subRegion.getReference< array2d< real64 > >( CellBlock::viewKeyStruct::elementCenterString );
    arrayView1d< real64 const > const & elemVolume =
      subRegion.getReference< array1d< real64 > >( CellBlock::viewKeyStruct::elementVolumeString );
    arrayView2d< real64 const > const & elemPerm =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseBase::viewKeyStruct::permeabilityString );

    // get the cell-centered depth
    arrayView1d< real64 const > const & elemGravCoef =
      subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseBase::viewKeyStruct::gravityCoefString );

    // assemble the residual and Jacobian element by element
    // in this loop we assemble both equation types: mass conservation in the elements and constraints at the faces
    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_DEVICE ( localIndex const ei )
    {

      // transmissibility matrix
      stackArray2d< real64, NF *NF > transMatrix( NF, NF );
      stackArray2d< real64, NF *NF > transMatrixGrav( NF, NF );

      real64 const perm[ 3 ] = { elemPerm[ei][0], elemPerm[ei][1], elemPerm[ei][2] };

      // recompute the local transmissibility matrix at each iteration
      // we can decide later to precompute transMatrix if needed
      IP_TYPE::template compute< NF >( nodePosition,
                                       transMultiplier,
                                       faceToNodes,
                                       elemToFaces[ei],
                                       elemCenter[ei],
                                       elemVolume[ei],
                                       perm,
                                       lengthTolerance,
                                       transMatrix );

      // currently the gravity term in the transport scheme is treated as in MRST, that is, always with TPFA
      // this is why below we have to recompute the TPFA transmissibility in addition to the transmissibility matrix above
      // TODO: treat the gravity term with a consistent inner product
      mimeticInnerProduct::TPFAInnerProduct::compute< NF >( nodePosition,
                                                            transMultiplier,
                                                            faceToNodes,
                                                            elemToFaces[ei],
                                                            elemCenter[ei],
                                                            elemVolume[ei],
                                                            perm,
                                                            lengthTolerance,
                                                            transMatrixGrav );

      // perform flux assembly in this element
      CompositionalMultiphaseHybridFVMKernels::AssemblerKernel::compute< NF, NC, NP >( er, esr, ei,
                                                                                       regionFilter,
                                                                                       elemRegionList,
                                                                                       elemSubRegionList,
                                                                                       elemList,
                                                                                       faceDofNumber,
                                                                                       faceGhostRank,
                                                                                       facePres,
                                                                                       dFacePres,
                                                                                       faceGravCoef,
                                                                                       mimFaceGravCoef,
                                                                                       elemToFaces[ei],
                                                                                       elemPres[ei],
                                                                                       dElemPres[ei],
                                                                                       elemGravCoef[ei],
                                                                                       phaseDens,
                                                                                       dPhaseDens_dPres,
                                                                                       dPhaseDens_dCompFrac,
                                                                                       phaseMassDens,
                                                                                       dPhaseMassDens_dPres,
                                                                                       dPhaseMassDens_dCompFrac,
                                                                                       phaseMob,
                                                                                       dPhaseMob_dPres,
                                                                                       dPhaseMob_dCompDens,
                                                                                       dCompFrac_dCompDens,
                                                                                       phaseCompFrac,
                                                                                       dPhaseCompFrac_dPres,
                                                                                       dPhaseCompFrac_dCompFrac,
                                                                                       elemDofNumber,
                                                                                       elemGhostRank[ei],
                                                                                       rankOffset,
                                                                                       dt,
                                                                                       transMatrix,
                                                                                       transMatrixGrav,
                                                                                       localMatrix,
                                                                                       localRhs );
    } );
  }

};


/******************************** PhaseMobilityKernel ********************************/

/**
 * @brief Functions to compute phase mobilities and derivatives from viscosity and relperm
 */
struct PhaseMobilityKernel
{
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  compute( arraySlice2d< real64 const > const & dCompFrac_dCompDens,
           arraySlice1d< real64 const > const & phaseVisc,
           arraySlice1d< real64 const > const & dPhaseVisc_dPres,
           arraySlice2d< real64 const > const & dPhaseVisc_dComp,
           arraySlice1d< real64 const > const & phaseRelPerm,
           arraySlice2d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
           arraySlice1d< real64 const > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64 const > const & dPhaseVolFrac_dComp,
           arraySlice1d< real64 > const & phaseMob,
           arraySlice1d< real64 > const & dPhaseMob_dPres,
           arraySlice2d< real64 > const & dPhaseMob_dComp );

  template< localIndex NC, localIndex NP >
  static void
  launch( localIndex const size,
          arrayView3d< real64 const > const & dCompFrac_dCompDens,
          arrayView3d< real64 const > const & phaseVisc,
          arrayView3d< real64 const > const & dPhaseVisc_dPres,
          arrayView4d< real64 const > const & dPhaseVisc_dComp,
          arrayView3d< real64 const > const & phaseRelPerm,
          arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const > const & dPhaseVolFrac_dComp,
          arrayView2d< real64 > const & phaseMob,
          arrayView2d< real64 > const & dPhaseMob_dPres,
          arrayView3d< real64 > const & dPhaseMob_dComp );

  template< localIndex NC, localIndex NP >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          arrayView3d< real64 const > const & dCompFrac_dCompDens,
          arrayView3d< real64 const > const & phaseVisc,
          arrayView3d< real64 const > const & dPhaseVisc_dPres,
          arrayView4d< real64 const > const & dPhaseVisc_dComp,
          arrayView3d< real64 const > const & phaseRelPerm,
          arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const > const & dPhaseVolFrac_dComp,
          arrayView2d< real64 > const & phaseMob,
          arrayView2d< real64 > const & dPhaseMob_dPres,
          arrayView3d< real64 > const & dPhaseMob_dComp );
};


/******************************** ResidualNormKernel ********************************/

struct ResidualNormKernel
{
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  template< typename POLICY, typename REDUCE_POLICY >
  static void
  launch( arrayView1d< real64 const > const & localResidual,
          globalIndex const rankOffset,
          localIndex const numPhases,
          arrayView1d< globalIndex const > const & facePresDofNumber,
          arrayView1d< integer const > const & faceGhostRank,
          SortedArrayView< localIndex const > const & regionFilter,
          arrayView2d< localIndex const > const & elemRegionList,
          arrayView2d< localIndex const > const & elemSubRegionList,
          arrayView2d< localIndex const > const & elemList,
          ElementViewConst< arrayView1d< real64 const > > const & elemVolume,
          ElementViewConst< arrayView2d< real64 const > > const & phaseMobOld,
          real64 & localResidualNorm )
  {
    RAJA::ReduceSum< REDUCE_POLICY, real64 > sumScaled( 0.0 );

    forAll< POLICY >( facePresDofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iface )
    {
      // if not ghost face and if adjacent to target region, increment the residual norm
      if( faceGhostRank[iface] < 0 && facePresDofNumber[iface] >= 0 )
      {
        real64 normalizer = 0;
        localIndex elemCounter = 0;
        for( localIndex k=0; k<elemRegionList.size( 1 ); ++k )
        {
          localIndex const er  = elemRegionList[iface][k];
          localIndex const esr = elemSubRegionList[iface][k];
          localIndex const ei  = elemList[iface][k];

          bool const onBoundary = (er == -1 || esr == -1 || ei == -1);
          bool const isInTarget = regionFilter.contains( er );

          if( !onBoundary && isInTarget )
          {
            // compute a normalizer to obtain a dimensionless norm
            real64 sumMobOld = 0.0;
            for( localIndex ip = 0; ip < numPhases; ++ip )
            {
              sumMobOld += phaseMobOld[er][esr][ei][ip];
            }
            real64 const totalMobOld = ( sumMobOld < 1e-3 ) ? 1e-3 : sumMobOld;
            normalizer += elemVolume[er][esr][ei] / totalMobOld;
            elemCounter++;
          }
        }
        normalizer /= elemCounter;

        localIndex const lid = LvArray::integerConversion< localIndex >( facePresDofNumber[iface] - rankOffset );
        // note: unit of localResidual[lid] * totalMobOld: m^3, so this is dimensionless
        real64 const val = localResidual[lid] / normalizer;
        sumScaled += val * val;
      }
    } );

    localResidualNorm = localResidualNorm + sumScaled.get();
  }

};

/******************************** SolutionCheckKernel ********************************/

struct SolutionCheckKernel
{

  template< typename POLICY, typename REDUCE_POLICY >
  static localIndex
  launch( arrayView1d< real64 const > const & localSolution,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & ghostRank,
          arrayView1d< real64 const > const & facePres,
          arrayView1d< real64 const > const & dFacePres,
          real64 const scalingFactor )
  {
    RAJA::ReduceMin< REDUCE_POLICY, integer > check( 1 );

    forAll< POLICY >( dofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iface )
    {
      if( ghostRank[iface] < 0 && dofNumber[iface] >= 0 )
      {
        localIndex const localRow = LvArray::integerConversion< localIndex >( dofNumber[iface] - rankOffset );
        {
          real64 const newFacePres = facePres[iface] + dFacePres[iface] + scalingFactor * localSolution[localRow];
          check.min( newFacePres >= 0.0 );
        }
      }
    } );
    return check.get();
  }

};

/******************************** PrecomputeKernel ********************************/

struct PrecomputeKernel
{

  template< typename IP_TYPE, localIndex NF >
  static void
  launch( localIndex const subRegionSize,
          localIndex const faceManagerSize,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
          ArrayOfArraysView< localIndex const > const & faceToNodes,
          arrayView2d< real64 const > const & elemCenter,
          arrayView1d< real64 const > const & elemVolume,
          arrayView2d< real64 const > const & elemPerm,
          arrayView1d< real64 const > const & elemGravCoef,
          arrayView2d< localIndex const > const & elemToFaces,
          arrayView1d< real64 const > const & transMultiplier,
          real64 const & lengthTolerance,
          arrayView1d< RAJA::ReduceSum< serialReduce, real64 > > const & mimFaceGravCoefNumerator,
          arrayView1d< RAJA::ReduceSum< serialReduce, real64 > > const & mimFaceGravCoefDenominator,
          arrayView1d< real64 > const & mimFaceGravCoef )
  {
    forAll< serialPolicy >( subRegionSize, [=] ( localIndex const ei )
    {
      stackArray2d< real64, NF *NF > transMatrix( NF, NF );

      real64 const perm[ 3 ] = { elemPerm[ei][0], elemPerm[ei][1], elemPerm[ei][2] };

      IP_TYPE::template compute< NF >( nodePosition,
                                       transMultiplier,
                                       faceToNodes,
                                       elemToFaces[ei],
                                       elemCenter[ei],
                                       elemVolume[ei],
                                       perm,
                                       lengthTolerance,
                                       transMatrix );

      for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
      {
        mimFaceGravCoefNumerator[elemToFaces[ei][ifaceLoc]] += elemGravCoef[ei] * transMatrix[ifaceLoc][ifaceLoc];
        mimFaceGravCoefDenominator[elemToFaces[ei][ifaceLoc]] += transMatrix[ifaceLoc][ifaceLoc];
      }
    } );

    forAll< serialPolicy >( faceManagerSize, [=] GEOSX_HOST_DEVICE ( localIndex const iface )
    {
      if( !isZero( mimFaceGravCoefDenominator[iface].get() ) )
      {
        mimFaceGravCoef[iface] = mimFaceGravCoefNumerator[iface].get() / mimFaceGravCoefDenominator[iface].get();
      }
    } );
  }
};

/******************************** Kernel switches ********************************/

namespace internal
{

template< typename T, typename LAMBDA >
void KernelLaunchSelectorFaceSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "KernelLaunchSelectorFaceSwitch: type should be integral" );

  switch( value )
  {
    case 4:
    { lambda( std::integral_constant< T, 4 >() ); return;}
    case 5:
    { lambda( std::integral_constant< T, 5 >() ); return;}
    case 6:
    { lambda( std::integral_constant< T, 6 >() ); return;}
    default: GEOSX_ERROR( "Unknown numFacesInElem value: " << value );
  }
}

} // namespace internal

template< typename IP_TYPE, typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector( localIndex numFacesInElem, localIndex numComps, localIndex numPhases, ARGS && ... args )
{
  internal::KernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
  {
    if( numPhases == 2 )
    {
      if( numComps == 2 )
      {
        KERNELWRAPPER::template launch< IP_TYPE, NF(), 2, 2 >( std::forward< ARGS >( args )... );
      }
      else if( numComps == 3 )
      {
        KERNELWRAPPER::template launch< IP_TYPE, NF(), 3, 2 >( std::forward< ARGS >( args )... );
      }
      else if( numComps == 4 )
      {
        KERNELWRAPPER::template launch< IP_TYPE, NF(), 4, 2 >( std::forward< ARGS >( args )... );
      }
      else if( numComps == 5 )
      {
        KERNELWRAPPER::template launch< IP_TYPE, NF(), 5, 2 >( std::forward< ARGS >( args )... );
      }
      else
      {
        GEOSX_ERROR( "Unsupported number of components: " << numComps );
      }
    }
    else if( numPhases == 3 )
    {
      if( numComps == 2 )
      {
        KERNELWRAPPER::template launch< IP_TYPE, NF(), 2, 3 >( std::forward< ARGS >( args )... );
      }
      else if( numComps == 3 )
      {
        KERNELWRAPPER::template launch< IP_TYPE, NF(), 3, 3 >( std::forward< ARGS >( args )... );
      }
      else if( numComps == 4 )
      {
        KERNELWRAPPER::template launch< IP_TYPE, NF(), 4, 3 >( std::forward< ARGS >( args )... );
      }
      else if( numComps == 5 )
      {
        KERNELWRAPPER::template launch< IP_TYPE, NF(), 5, 3 >( std::forward< ARGS >( args )... );
      }
      else
      {
        GEOSX_ERROR( "Unsupported number of components: " << numComps );
      }
    }
    else
    {
      GEOSX_ERROR( "Unsupported number of phases: " << numPhases );
    }
  } );
}

} // namespace CompositionalMultiphaseHybridFVMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMKERNELS_HPP
