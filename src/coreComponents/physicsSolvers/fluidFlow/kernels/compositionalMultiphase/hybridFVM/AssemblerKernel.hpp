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
 * @file AssemblerKernel.hpp
 */

#ifndef GEOSX_ASSEMBLERKERNEL_HPP
#define GEOSX_ASSEMBLERKERNEL_HPP

#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "mesh/ElementRegionManager.hpp"


namespace geos
{

namespace compositionalMultiphaseHybridFVMKernels
{

using namespace constitutive;

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
   * @param[in] faceGravCoef the depth at the mesh faces
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] elemPres the pressure at this element's center
   * @param[in] elemGravCoef the depth at this element's center
   * @param[in] phaseDens the phase densities in the domain (non-local)
   * @param[in] dPhaseDens the derivatives of the phase densities in the domain wrt pressure and component fractions (non-local)
   * @param[in] phaseMob the phase mobilities in the domain (non-local)
   * @param[in] dPhaseMob the derivatives of the phase mobilities in the domain (non-local)
   * @param[in] dCompFrac_dCompDens the derivatives of the component fraction in the domain wrt component density (non-local)
   * @param[in] phaseCompFrac the phase component fractions in the domain (non-local)
   * @param[in] dPhaseCompFrac the derivatives of the phase component fractions in the domain wrt pressure and component fractions
   *(non-local)
   * @param[in] elemDofNumber the dof number of the cell centered pressures (non-local)
   * @param[in] elemGhostRank the ghost rank of the element in which we assemble the fluxes
   * @param[in] rankOffset the offset of this rank
   * @param[in] lengthTolerance tolerance used in the transmissibility matrix computation
   * @param[in] dt time step size
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[inout] matrix the system matrix
   * @param[inout] rhs the system right-hand side vector
   */
  template< integer NF, integer NC, integer NP >
  GEOS_HOST_DEVICE
  static void
  compute( localIndex const er,
           localIndex const esr,
           localIndex const ei,
           SortedArrayView< localIndex const > const & regionFilter,
           arrayView2d< localIndex const > const & elemRegionList,
           arrayView2d< localIndex const > const & elemSubRegionList,
           arrayView2d< localIndex const > const & elemList,
           arrayView1d< globalIndex const > const & faceDofNumber,
           arrayView1d< integer const > const & faceGhostRank,
           arrayView1d< real64 const > const & facePres,
           arrayView1d< real64 const > const & faceGravCoef,
           arrayView1d< real64 const > const & mimFaceGravCoef,
           arraySlice1d< localIndex const > const & elemToFaces,
           real64 const & elemPres,
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
           integer const elemGhostRank,
           globalIndex const rankOffset,
           real64 const & dt,
           arraySlice2d< real64 const > const & transMatrix,
           arraySlice2d< real64 const > const & transMatrixGrav,
           CRSMatrixView< real64, globalIndex const > const & localMatrix,
           arrayView1d< real64 > const & localRhs );
};

}

}

#endif //GEOSX_ASSEMBLERKERNEL_HPP
