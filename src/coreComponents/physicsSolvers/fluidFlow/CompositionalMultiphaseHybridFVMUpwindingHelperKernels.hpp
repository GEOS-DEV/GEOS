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
 * @file CompositionalMultiphaseHybridFVMHelperKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMUPWINDINGHELPERKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMUPWINDINGHELPERKERNELS_HPP

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mesh/MeshLevel.hpp"
#include "finiteVolume/HybridFVMInnerProduct.hpp"

namespace geosx
{

namespace CompositionalMultiphaseHybridFVMKernels
{


/******************************** UpwindingHelper ********************************/

struct UpwindingHelper
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
   * @brief At a given one-sided face, compute the upwind viscous transport coefficient
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] neighborIds region, subRegion, and element indices of the neigbhbor element
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
   * @param[in] oneSidedVolFlux the total volumetric flux at this face
   * @param[out] upwPhaseViscCoef the upwind viscous transport coef at this face
   * @param[out] dUpwPhaseViscCoef_dPres the derivative of the upwind viscous transport coef wrt pressure at this face
   * @param[out] dUpwPhaseViscCoef_dCompDens the derivatives of the upwind viscous transport coef wrt component density at this face
   * @param[out] upwViscDofNumber the dof number of the upwind cell at this face
   */
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
    UpwindViscousCoefficient( localIndex const (&localIds)[ 3 ],
                              localIndex const (&neighborIds)[ 3 ],
                              ElementViewConst< arrayView3d< real64 const > > const & phaseDens,
                              ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres,
                              ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dCompFrac,
                              ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
                              ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
                              ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dCompDens,
                              ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                              ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac,
                              ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
                              ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dCompFrac,
                              ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
                              real64 const & oneSidedVolFlux,
                              real64 ( &upwPhaseViscCoef )[ NP ][ NC ],
                              real64 ( &dUpwPhaseViscCoef_dPres )[ NP ][ NC ],
                              real64 ( &dUpwPhaseViscCoef_dCompDens )[ NP ][ NC ][ NC ],
                              globalIndex & upwViscDofNumber );

  /**
   * @brief At a given one-sided face, compute the upwind viscous transport coefficient
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] neighborIds region, subRegion, and element indices of the neigbhbor element
   * @param[in] transGravCoef
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
   * @param[inout] phaseGravTerm the gravCoef multiplied by the difference in phase densities
   * @param[inout] dPhaseGravTerm_dPres the derivatives of the gravCoef multiplied by the difference in phase densities wrt pressure
   * @param[inout] dPhaseGravTerm_dCompDens the derivatives of the gravCoef multiplied by the difference in phase densities wrt component
   * density
   * @param[inout] upwPhaseGravCoef the upwinded buoyancy transport coefficient at this face (ifaceLoc)
   * @param[inout] dUpwPhaseGravCoef_dPres the derivative of the upwinded buoyancy transport coefficient wrt pressure
   * @param[inout] dUpwPhaseGravCoef_dCompDens the derivative of the upwinded buoyancy transport coefficient wrt component density
   */
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
    UpwindBuoyancyCoefficient( localIndex const (&localIds)[ 3 ],
                               localIndex const (&neighborIds)[ 3 ],
                               real64 const & transGravCoef,
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
                               real64 ( &phaseGravTerm )[ NP ][ NP-1 ],
                               real64 ( &dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ],
                               real64 ( &dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ],
                               real64 ( &upwPhaseGravCoef )[ NP ][ NP-1 ][ NC ],
                               real64 ( &dUpwPhaseGravCoef_dPres )[ NP ][ NP-1 ][ NC ][ 2 ],
                               real64 ( &dUpwPhaseGravCoef_dCompDens )[ NP ][ NP-1 ][ NC ][ 2 ][ NC ] );

  /**
   * @brief At a given one-sided face, compute the gravCoef multiplied by the difference in phase densities
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] neighborIds region, subRegion, and element indices of the neigbhbor element
   * @param[in] transGravCoef
   * @param[in] phaseDens the phase densities in the domain (non-local)
   * @param[in] dPhaseDens_dPres the derivatives of the phase densities in the domain wrt pressure (non-local)
   * @param[in] dPhaseDens_dCompFrac the derivatives of the phase densities in the domain wrt component fraction (non-local)
   * @param[in] dCompFrac_dCompDens the derivatives of the component fraction in the domain wrt component density (non-local)
   * @param[inout] phaseGravTerm the gravCoef multiplied by the difference in phase densities
   * @param[inout] dPhaseGravTerm_dPres the derivatives of the gravCoef multiplied by the difference in phase densities wrt pressure
   * @param[inout] dPhaseGravTerm_dCompDens the derivatives of the gravCoef multiplied by the difference in phase densities wrt component
   * density
   */
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
    ComputePhaseGravTerm( localIndex const (&localIds)[ 3 ],
                          localIndex const (&neighborIds)[ 3 ],
                          real64 const & transGravCoef,
                          ElementViewConst< arrayView3d< real64 const > > const & phaseDens,
                          ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres,
                          ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dCompFrac,
                          ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                          real64 ( &phaseGravTerm )[ NP ][ NP-1 ],
                          real64 ( &dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ],
                          real64 ( &dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ] );

  /**
   * @brief At a given one-sided face, compute the upwinded total mobility
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] neighborIds region, subRegion, and element indices of the neigbhbor element
   * @param[in] phaseMob the phase mobilities in the domain (non-local)
   * @param[in] dPhaseMob_dPres the derivatives of the phase mobilities in the domain wrt pressure (non-local)
   * @param[in] dPhaseMob_dCompDens the derivatives of the phase mobilities in the domain wrt component fraction (non-local)
   * @param[in] phaseGravTerm the gravCoef multiplied by the difference in phase densities
   * @param[inout] totalMob the upwinded total mobility
   * @param[inout] dTotalMob_dPres the derivative of the upwinded total mobility wrt pressure
   * @param[inout] dTotalMob_dCompDens the derivative of the upwinded total mobility wrt component density
   */
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
    ComputeUpwindedTotalMobility( localIndex const (&localIds)[ 3 ],
                                  localIndex const (&neighborIds)[ 3 ],
                                  ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
                                  ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
                                  ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dCompDens,
                                  real64 const (&phaseGravTerm)[ NP ][ NP-1 ],
                                  real64 & totalMob,
                                  real64 ( &dTotalMob_dPres )[ 2 ],
                                  real64 ( &dTotalMob_dCompDens )[ 2 ][ NC ] );

  /**
   * @brief Set the element indices used to evaluate the mobility ratios of the buoyancy term in hybrid upwinding
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] neighborIds region, subRegion, and element indices of the neigbhbor element
   * @param[in] gravTerm the gravCoef multiplied by the difference in phase densities
   * @param[in] totalMob the upwinded total mobility
   * @param[in] eru region index of the upwind element
   * @param[in] esru subRegion index of the upwind element
   * @param[in] eiu element index of the upwind element
   * @param[in] posu position (local or neighbor) of the upwind element
   * @param[in] erd region index of the downwind element
   * @param[in] esrd subRegion index of the downwind element
   * @param[in] eid element index of the downwind element
   * @param[in] posd position (local or neighbor) of the downwind element
   */
  GEOSX_HOST_DEVICE
  static void
  SetIndicesForMobilityRatioUpwinding( localIndex const (&localIds)[ 3 ],
                                       localIndex const (&neighborIds)[ 3 ],
                                       real64 const & gravTerm,
                                       localIndex & eru, localIndex & esru, localIndex & eiu, localIndex & posu,
                                       localIndex & erd, localIndex & esrd, localIndex & eid, localIndex & posd );

  /**
   * @brief Set the element indices used to evaluate the total mobility of the buoyancy term in hybrid upwinding
   * @param[in] localIds triplet of indices for the local element
   * @param[in] neighborIds triplet of indices for the neighbor element
   * @param[in] gravTerm gravity term used to upwind
   * @param[out] totalMobIds for each phase, triplet of indices of the upwind element
   * @param[out] totalMobPos for each phase, flag specifying with the upwind element is local or neighbor
   */
  template< localIndex NP >
  GEOSX_HOST_DEVICE
  static void
    SetIndicesForTotalMobilityUpwinding( localIndex const (&localIds)[ 3 ],
                                         localIndex const (&neighborIds)[ 3 ],
                                         real64 const (&gravTerm)[ NP ][ NP-1 ],
                                         localIndex ( &totalMobIds )[ NP ][ 3 ],
                                         localIndex ( &totalMobPos )[ NP ] );

};

} // namespace CompositionalMultiphaseHybridFVMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMUPWINDINGHELPERKERNELS_HPP
