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
 * @file AssemblerKernel.cpp
 */

#include "AssemblerKernel.hpp"

#include "AssemblerKernelHelper.hpp"

namespace geos
{

namespace compositionalMultiphaseHybridFVMKernels
{

/******************************** AssemblerKernel ********************************/

template< integer NF, integer NC, integer NP >
GEOS_HOST_DEVICE
inline
void
AssemblerKernel::
  compute( localIndex const er, localIndex const esr, localIndex const ei,
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
                                                      faceGravCoef,
                                                      elemToFaces,
                                                      elemPres,
                                                      elemGravCoef,
                                                      phaseMassDens[er][esr][ei][0],
                                                      dPhaseMassDens[er][esr][ei][0],
                                                      phaseMob[er][esr][ei],
                                                      dPhaseMob[er][esr][ei],
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
                                                                 dPhaseDens,
                                                                 phaseMassDens,
                                                                 dPhaseMassDens,
                                                                 phaseMob,
                                                                 dPhaseMob,
                                                                 dCompFrac_dCompDens,
                                                                 phaseCompFrac,
                                                                 dPhaseCompFrac,
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

#define INST_AssemblerKernel( NF, NC, NP ) \
  template \
  GEOS_HOST_DEVICE \
  void \
  AssemblerKernel:: \
    compute< NF, NC, NP >( localIndex const er, localIndex const esr, localIndex const ei, \
                           SortedArrayView< localIndex const > const & regionFilter, \
                           arrayView2d< localIndex const > const & elemRegionList, \
                           arrayView2d< localIndex const > const & elemSubRegionList, \
                           arrayView2d< localIndex const > const & elemList, \
                           arrayView1d< globalIndex const > const & faceDofNumber, \
                           arrayView1d< integer const > const & faceGhostRank, \
                           arrayView1d< real64 const > const & facePres, \
                           arrayView1d< real64 const > const & faceGravCoef, \
                           arrayView1d< real64 const > const & mimFaceGravCoef, \
                           arraySlice1d< localIndex const > const & elemToFaces, \
                           real64 const & elemPres, \
                           real64 const & elemGravCoef, \
                           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens, \
                           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens, \
                           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens, \
                           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens, \
                           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob, \
                           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob, \
                           ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens, \
                           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac, \
                           ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac, \
                           ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber, \
                           integer const elemGhostRank, \
                           globalIndex const rankOffset, \
                           real64 const & dt, \
                           arraySlice2d< real64 const > const & transMatrix, \
                           arraySlice2d< real64 const > const & transMatrixGrav, \
                           CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                           arrayView1d< real64 > const & localRhs )

INST_AssemblerKernel( 4, 1, 2 );
INST_AssemblerKernel( 4, 2, 2 );
INST_AssemblerKernel( 4, 3, 2 );
INST_AssemblerKernel( 4, 4, 2 );
INST_AssemblerKernel( 4, 5, 2 );

INST_AssemblerKernel( 4, 1, 3 );
INST_AssemblerKernel( 4, 2, 3 );
INST_AssemblerKernel( 4, 3, 3 );
INST_AssemblerKernel( 4, 4, 3 );
INST_AssemblerKernel( 4, 5, 3 );

INST_AssemblerKernel( 5, 1, 2 );
INST_AssemblerKernel( 5, 2, 2 );
INST_AssemblerKernel( 5, 3, 2 );
INST_AssemblerKernel( 5, 4, 2 );
INST_AssemblerKernel( 5, 5, 2 );

INST_AssemblerKernel( 5, 1, 3 );
INST_AssemblerKernel( 5, 2, 3 );
INST_AssemblerKernel( 5, 3, 3 );
INST_AssemblerKernel( 5, 4, 3 );
INST_AssemblerKernel( 5, 5, 3 );

INST_AssemblerKernel( 6, 1, 2 );
INST_AssemblerKernel( 6, 2, 2 );
INST_AssemblerKernel( 6, 3, 2 );
INST_AssemblerKernel( 6, 4, 2 );
INST_AssemblerKernel( 6, 5, 2 );

INST_AssemblerKernel( 6, 1, 3 );
INST_AssemblerKernel( 6, 2, 3 );
INST_AssemblerKernel( 6, 3, 3 );
INST_AssemblerKernel( 6, 4, 3 );
INST_AssemblerKernel( 6, 5, 3 );

#undef INST_AssemblerKernel

}

}
