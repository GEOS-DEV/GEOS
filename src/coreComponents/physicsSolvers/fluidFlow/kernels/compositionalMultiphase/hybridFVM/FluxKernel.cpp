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
 * @file CompositionalMultiphaseHybridFVMKernels.cpp
 */

#include "FluxKernel.hpp"

#include "AssemblerKernel.hpp"
#include "finiteVolume/mimeticInnerProducts/BdVLMInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/TPFAInnerProduct.hpp"
#include "physicsSolvers/fluidFlow/fields/FlowSolverBaseFields.hpp"


namespace geos
{

namespace compositionalMultiphaseHybridFVMKernels
{

/******************************** FluxKernel ********************************/

template< integer NF, integer NC, integer NP, typename IP_TYPE >
void
FluxKernel::
  launch( localIndex er, localIndex esr,
          CellElementSubRegion const & subRegion,
          constitutive::PermeabilityBase const & permeabilityModel,
          SortedArrayView< localIndex const > const & regionFilter,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
          arrayView2d< localIndex const > const & elemRegionList,
          arrayView2d< localIndex const > const & elemSubRegionList,
          arrayView2d< localIndex const > const & elemList,
          ArrayOfArraysView< localIndex const > const & faceToNodes,
          arrayView1d< globalIndex const > const & faceDofNumber,
          arrayView1d< integer const > const & faceGhostRank,
          arrayView1d< real64 const > const & facePres,
          arrayView1d< real64 const > const & faceGravCoef,
          arrayView1d< real64 const > const & mimFaceGravCoef,
          arrayView1d< real64 const > const & transMultiplier,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
          ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
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
    subRegion.getReference< array1d< real64 > >( fields::flow::pressure::key() );

  // get the element data needed for transmissibility computation
  arrayView2d< real64 const > const & elemCenter =
    subRegion.getReference< array2d< real64 > >( CellElementSubRegion::viewKeyStruct::elementCenterString() );
  arrayView1d< real64 const > const & elemVolume =
    subRegion.getReference< array1d< real64 > >( CellElementSubRegion::viewKeyStruct::elementVolumeString() );

  // TODO add this dependency to the compute function
  //arrayView3d< real64 const > const elemdPermdPres = permeabilityModel.dPerm_dPressure();

  arrayView3d< real64 const > const & elemPerm = permeabilityModel.permeability();

  // get the cell-centered depth
  arrayView1d< real64 const > const & elemGravCoef =
    subRegion.getReference< array1d< real64 > >( fields::flow::gravityCoefficient::key() );

  // assemble the residual and Jacobian element by element
  // in this loop we assemble both equation types: mass conservation in the elements and constraints at the faces
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_DEVICE ( localIndex const ei )
  {
    // transmissibility matrix
    stackArray2d< real64, NF *NF > transMatrix( NF, NF );
    stackArray2d< real64, NF *NF > transMatrixGrav( NF, NF );

    real64 const perm[ 3 ] = { elemPerm[ei][0][0], elemPerm[ei][0][1], elemPerm[ei][0][2] };

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
    compositionalMultiphaseHybridFVMKernels::AssemblerKernel::compute< NF, NC, NP >( er, esr, ei,
                                                                                     regionFilter,
                                                                                     elemRegionList,
                                                                                     elemSubRegionList,
                                                                                     elemList,
                                                                                     faceDofNumber,
                                                                                     faceGhostRank,
                                                                                     facePres,
                                                                                     faceGravCoef,
                                                                                     mimFaceGravCoef,
                                                                                     elemToFaces[ei],
                                                                                     elemPres[ei],
                                                                                     elemGravCoef[ei],
                                                                                     phaseDens,
                                                                                     dPhaseDens,
                                                                                     phaseMassDens,
                                                                                     dPhaseMassDens,
                                                                                     phaseMob,
                                                                                     dPhaseMob,
                                                                                     dCompFrac_dCompDens,
                                                                                     phaseCompFrac,
                                                                                     dPhaseCompFrac,
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

#define INST_FluxKernel( NF, NC, NP, IP_TYPE ) \
  template \
  void \
  FluxKernel:: \
    launch< NF, NC, NP, IP_TYPE >( localIndex er, localIndex esr, \
                                   CellElementSubRegion const & subRegion, \
                                   constitutive::PermeabilityBase const & permeabilityModel, \
                                   SortedArrayView< localIndex const > const & regionFilter, \
                                   arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition, \
                                   arrayView2d< localIndex const > const & elemRegionList, \
                                   arrayView2d< localIndex const > const & elemSubRegionList, \
                                   arrayView2d< localIndex const > const & elemList, \
                                   ArrayOfArraysView< localIndex const > const & faceToNodes, \
                                   arrayView1d< globalIndex const > const & faceDofNumber, \
                                   arrayView1d< integer const > const & faceGhostRank, \
                                   arrayView1d< real64 const > const & facePres, \
                                   arrayView1d< real64 const > const & faceGravCoef, \
                                   arrayView1d< real64 const > const & mimFaceGravCoef, \
                                   arrayView1d< real64 const > const & transMultiplier, \
                                   ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob, \
                                   ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob, \
                                   ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens, \
                                   ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens, \
                                   ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens, \
                                   ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens, \
                                   ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens, \
                                   ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac, \
                                   ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac, \
                                   ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber, \
                                   globalIndex const rankOffset, \
                                   real64 const lengthTolerance, \
                                   real64 const dt, \
                                   CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                                   arrayView1d< real64 > const & localRhs )

INST_FluxKernel( 4, 1, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 2, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 3, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 4, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 5, 2, mimeticInnerProduct::TPFAInnerProduct const );

INST_FluxKernel( 4, 1, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 2, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 3, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 4, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 5, 3, mimeticInnerProduct::TPFAInnerProduct const );

INST_FluxKernel( 5, 1, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 2, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 3, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 4, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 5, 2, mimeticInnerProduct::TPFAInnerProduct const );

INST_FluxKernel( 5, 1, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 2, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 3, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 4, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 5, 3, mimeticInnerProduct::TPFAInnerProduct const );

INST_FluxKernel( 6, 1, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 2, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 3, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 4, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 5, 2, mimeticInnerProduct::TPFAInnerProduct const );

INST_FluxKernel( 6, 1, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 2, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 3, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 4, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 5, 3, mimeticInnerProduct::TPFAInnerProduct const );


INST_FluxKernel( 4, 1, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 2, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 3, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 4, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 5, 2, mimeticInnerProduct::BdVLMInnerProduct const );

INST_FluxKernel( 4, 1, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 2, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 3, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 4, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 5, 3, mimeticInnerProduct::BdVLMInnerProduct const );

INST_FluxKernel( 5, 1, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 2, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 3, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 4, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 5, 2, mimeticInnerProduct::BdVLMInnerProduct const );

INST_FluxKernel( 5, 1, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 2, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 3, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 4, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 5, 3, mimeticInnerProduct::BdVLMInnerProduct const );

INST_FluxKernel( 6, 1, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 2, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 3, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 4, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 5, 2, mimeticInnerProduct::BdVLMInnerProduct const );

INST_FluxKernel( 6, 1, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 2, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 3, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 4, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 5, 3, mimeticInnerProduct::BdVLMInnerProduct const );


#undef INST_FluxKernel

} // namespace compositionalMultiphaseHybridFVMKernels

} // namespace geos
