/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file AssemblerKernelHelper.hpp
 */

#ifndef GEOSX_ASSEMBLERKERNELHELPER_HPP
#define GEOSX_ASSEMBLERKERNELHELPER_HPP

#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "mesh/ElementRegionManager.hpp"


namespace geos
{

namespace compositionalMultiphaseHybridFVMKernels
{

using namespace constitutive;

/******************************** AssemblerKernelHelper ********************************/

struct AssemblerKernelHelper
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
   * @brief In a given element, compute the transmissibility-weighted pressure gradients in the cell
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] faceGravCoef the depth at the mesh facesb
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] elemPres the pressure at this element's center
   * @param[in] elemGravCoef the depth at this element's center
   * @param[in] phaseDens the phase densities in the element
   * @param[in] dPhaseDens the derivatives of the phase densities in the element wrt pressure and component fractions
   * @param[in] phaseMob the phase mobilities in the element
   * @param[in] dPhaseMob the derivatives of the phase mobilities in the element
   * @param[in] dCompFrac_dCompDens the derivatives of the component fraction in the element wrt component density
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[out] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[out] dOneSidedVolFlux_dPres the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[out] dOneSidedVolFlux_dFacePres the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[out] dOneSidedVolFlux_dCompDens the derivatives of the vol fluxes wrt to this element's component density
   */
  template< integer NF, integer NC, integer NP >
  GEOS_HOST_DEVICE
  static void
    applyGradient( arrayView1d< real64 const > const & facePres,
                   arrayView1d< real64 const > const & faceGravCoef,
                   arraySlice1d< localIndex const > const & elemToFaces,
                   real64 const & elemPres,
                   real64 const & elemGravCoef,
                   arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & elemPhaseMassDens,
                   arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dElemPhaseMassDens,
                   arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & elemPhaseMob,
                   arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dElemPhaseMob,
                   arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dElemCompFrac_dCompDens,
                   arraySlice2d< real64 const > const & transMatrix,
                   real64 ( &oneSidedVolFlux )[NF],
                   real64 ( &dOneSidedVolFlux_dPres )[NF],
                   real64 ( &dOneSidedVolFlux_dFacePres )[NF][NF],
                   real64 ( &dOneSidedVolFlux_dCompDens )[NF][NC] );

  /**
   * @brief In a given element, compute the flux divergence, i.e, sum the fluxes at this element's faces
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] rankOffset offset of this rank
   * @param[in] faceDofNumber the dof numbers of the face pressures in the domain
   * @param[in] elemRegionList map from face to element region index
   * @param[in] elemSubRegionList map from face to element subRegion index
   * @param[in] elemList map from face to element index
   * @param[in] regionFilter set of target regions of the solver
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] phaseDens phaseDens the phase densities in the domain
   * @param[in] dPhaseDens the derivatives of the phase densities wrt pressure and component fractions
   * @param[in] phaseMob the phase mobilities in the domain
   * @param[in] dPhaseMob the derivatives of the phase mobilities in the element
   * @param[in] dCompFrac_dCompDens the derivatives of the component fractions wrt component density
   * @param[in] phaseCompFrac the phase component fractions in the domain
   * @param[in] dPhaseCompFrac the derivatives of the phase component fractions wrt pressure and component fractions
   * @param[in] elemDofNumber the dof numbers of the element in the domain
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[in] dOneSidedVolFlux_dPres the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[in] dOneSidedVolFlux_dFacePres the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[in] dOneSidedVolFlux_dCompDens the derivatives of the vol fluxes wrt to this element's component density
   * @param[in] dt the time step
   * @param[inout] localMatrix the Jacobian matrix
   * @param[inout] localRhs the residual
   */
  template< integer NF, integer NC, integer NP >
  GEOS_HOST_DEVICE
  static void
  assembleFluxDivergence( localIndex const (&localIds)[3],
                          globalIndex const rankOffset,
                          arrayView2d< localIndex const > const & elemRegionList,
                          arrayView2d< localIndex const > const & elemSubRegionList,
                          arrayView2d< localIndex const > const & elemList,
                          SortedArrayView< localIndex const > const & regionFilter,
                          arrayView1d< globalIndex const > const & faceDofNumber,
                          arrayView1d< real64 const > const & mimFaceGravCoef,
                          arraySlice1d< localIndex const > const & elemToFaces,
                          real64 const & elemGravCoef,
                          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
                          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens,
                          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
                          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
                          ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
                          ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
                          arraySlice2d< real64 const > const & transMatrixGrav,
                          real64 const (&oneSidedVolFlux)[NF],
                          real64 const (&dOneSidedVolFlux_dPres)[NF],
                          real64 const (&dOneSidedVolFlux_dFacePres)[NF][NF],
                          real64 const (&dOneSidedVolFlux_dCompDens)[NF][NC],
                          real64 const & dt,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs );

  /**
   * @brief In a given element, compute the viscous flux divergence, i.e, sum the viscous fluxes at this element's faces
   * @param[in] ifaceLoc the local index of the face
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[in] dOneSidedVolFlux_dPres the derivatives of the vol fluxes wrt to this element's pressure
   * @param[in] dOneSidedVolFlux_dFacePres the derivatives of the vol fluxes wrt to this element's face centered pressure
   * @param[in] dOneSidedVolFlux_dCompDens the derivatives of the vol fluxes wrt to this element's component density
   * @param[in] upwPhaseViscCoef the upwinded viscous transport coefficient at this face (ifaceLoc)
   * @param[in] dUpwPhaseViscCoef_dPres the derivative of the upwinded viscous transport coefficient wrt pressure
   * @param[in] dUpwPhaseViscCoef_dCompDens the derivative of the upwinded viscous transport coefficient wrt component density
   * @param[in] upwViscDofNumber degree of freedom number of the upwind element
   * @param[in] faceDofNumber degree of freedom number of the face
   * @param[in] dt the time step
   * @param[inout] divMassFluxes the divergence of the fluxes in the element
   * @param[inout] dDivMassFluxes_dElemVars the derivatives of the flux divergence wrt the element centered vars (pres and comp dens)
   * @param[inout] dDivMassFluxes_dFaceVars the derivatives of the flux divergence wrt the face centered vars
   * @param[inout] dofColIndicesElemVars degrees of freedom of the cells involved in the flux divergence
   * @param[inout] dofColIndicesFaceVars degrees of freedom of the faces involved in the flux divergence
   */
  template< integer NF, integer NC, integer NP >
  GEOS_HOST_DEVICE
  static void
    assembleViscousFlux( localIndex const ifaceLoc,
                         real64 const (&oneSidedVolFlux)[NF],
                         real64 const (&dOneSidedVolFlux_dPres)[NF],
                         real64 const (&dOneSidedVolFlux_dFacePres)[NF][NF],
                         real64 const (&dOneSidedVolFlux_dCompDens)[NF][NC],
                         real64 const (&upwPhaseViscCoef)[NP][NC],
                         real64 const (&dUpwPhaseViscCoef_dPres)[NP][NC],
                         real64 const (&dUpwPhaseViscCoef_dCompDens)[NP][NC][NC],
                         globalIndex const elemDofNumber,
                         globalIndex const neighborDofNumber,
                         globalIndex const upwViscDofNumber,
                         globalIndex const faceDofNumber,
                         real64 const & dt,
                         real64 ( &divMassFluxes )[NC],
                         real64 ( &dDivMassFluxes_dElemVars )[NC][( NC + 1 ) * ( NF + 1 )],
                         real64 ( &dDivMassFluxes_dFaceVars )[NC][NF],
                         globalIndex ( &dofColIndicesElemVars )[( NC + 1 ) * ( NF + 1 )],
                         globalIndex ( &dofColIndicesFaceVars )[NF] );

  /**
   * @brief In a given element, compute the buoyancy flux divergence, i.e, sum the buoyancy fluxes at this element's faces
   * @param[in] ifaceLoc the local index of the face
   * @param[in] phaseGravTerm the gravCoef multiplied by the difference in phase densities
   * @param[in] dPhaseGravTerm_dPres the derivatives of the gravCoef multiplied by the difference in phase densities wrt pressure
   * @param[in] dPhaseGravTerm_dCompDens the derivatives of the gravCoef multiplied by the difference in phase densities wrt component
   * density
   * @param[in] upwPhaseGravCoef the upwinded buoyancy transport coefficient at this face (ifaceLoc)
   * @param[in] dUpwPhaseGravCoef_dPres the derivative of the upwinded buoyancy transport coefficient wrt pressure
   * @param[in] dUpwPhaseGravCoef_dCompDens the derivative of the upwinded buoyancy transport coefficient wrt component density
   * @param[in] neighborDofNumber the degree of freedom number of the neighbor element
   * @param[in] dt the time step
   * @param[inout] divMassFluxes the divergence of the fluxes in the element
   * @param[inout] dDivMassFluxes_dElemVars the derivatives of the flux divergence wrt the element centered vars (pres and comp dens)
   * @param[inout] dofColIndicesElemVars degrees of freedom of the cells involved in the flux divergence
   */
  template< integer NF, integer NC, integer NP >
  GEOS_HOST_DEVICE
  static void
    assembleBuoyancyFlux( localIndex const ifaceLoc,
                          real64 const (&phaseGravTerm)[NP][NP - 1],
                          real64 const (&dPhaseGravTerm_dPres)[NP][NP - 1][2],
                          real64 const (&dPhaseGravTerm_dCompDens)[NP][NP - 1][2][NC],
                          real64 const (&upwPhaseGravCoef)[NP][NP - 1][NC],
                          real64 const (&dUpwPhaseGravCoef_dPres)[NP][NP - 1][NC][2],
                          real64 const (&dUpwPhaseGravCoef_dCompDens)[NP][NP - 1][NC][2][NC],
                          real64 const & dt,
                          real64 ( &divMassFluxes )[NC],
                          real64 ( &dDivMassFluxes_dElemVars )[NC][( NC + 1 ) * ( NF + 1 )] );

  /**
   * @brief In a given element, assemble the constraints at this element's faces
   * @param[in] faceDofNumber the dof numbers of the face pressures in the subRegion
   * @param[in] faceGhostRank the ghost ranks of the face pressures in the subRegion
   * @param[in] elemToFaces the map from one-sided face to face to access face Dof numbers
   * @param[in] elemDofNumber the dof number of this element's cell centered pressure
   * @param[in] rankOffset the offset of this rank
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[in] dOneSidedVolFlux_dPres the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[in] dOneSidedVolFlux_dFacePres the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[in] dOneSidedVolFlux_dCompDens the derivatives of the vol fluxes wrt to this element's component density
   * @param[inout] matrix the jacobian matrix
   * @param[inout] rhs the residual
   */
  template< integer NF, integer NC, integer NP >
  GEOS_HOST_DEVICE
  static void
  assembleFaceConstraints( arrayView1d< globalIndex const > const & faceDofNumber,
                           arrayView1d< integer const > const & faceGhostRank,
                           arraySlice1d< localIndex const > const & elemToFaces,
                           globalIndex const elemDofNumber,
                           globalIndex const rankOffset,
                           real64 const (&oneSidedVolFlux)[NF],
                           real64 const (&dOneSidedVolFlux_dPres)[NF],
                           real64 const (&dOneSidedVolFlux_dFacePres)[NF][NF],
                           real64 const (&dOneSidedVolFlux_dCompDens)[NF][NC],
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs );

};

}

}

#endif //GEOSX_ASSEMBLERKERNELHELPER_HPP
