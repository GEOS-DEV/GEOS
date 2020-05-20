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
 * @file ProppantTransportKernels.hpp
 */

#ifndef GEOSX_PROPPANTTRANSPORTKERNELS_HPP
#define GEOSX_PROPPANTTRANSPORTKERNELS_HPP

#include "common/DataTypes.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

#include "physicsSolvers/fluidFlow/ProppantTransport.hpp"

namespace geosx
{

namespace ProppantTransportKernels
{

/******************************** AccumulationKernel ********************************/

struct AccumulationKernel
{

  inline static void
  Compute( localIndex const NC,
           real64 const & proppantConcOld,
           real64 const & proppantConcNew,
           arraySlice1d< real64 const > const & componentDensOld,
           arraySlice1d< real64 const > const & componentDensNew,
           arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM( dCompDens_dPres ),
           arraySlice2d< real64 const > const & dCompDens_dCompConc,
           real64 const & volume,
           real64 const & packPoreVolume,
           arraySlice1d< real64 > const & localAccum,
           arraySlice2d< real64 > const & localAccumJacobian )
  {

    // proppant mass conservation
    localAccum[0] = (proppantConcNew - proppantConcOld) * volume;

    for( localIndex c1 = 0; c1 < NC; ++c1 )
      for( localIndex c2 = 0; c2 < NC; ++c2 )
        localAccumJacobian[c1][c2] = 0.0;


    localAccumJacobian[0][0] = volume;

    // component mass conservation
    for( localIndex c1 = 0; c1 < NC; ++c1 )
    {

      localAccum[c1+1] = ( componentDensNew[c1] * (1.0 - proppantConcNew) - componentDensOld[c1] * (1.0 - proppantConcOld) ) * volume +
                         (componentDensNew[c1] - componentDensOld[c1]) * packPoreVolume;

      for( localIndex c2 = 0; c2 < NC; ++c2 )
        localAccumJacobian[c1+1][c2+1] = dCompDens_dCompConc[c1][c2] * (1.0 - proppantConcNew) * volume + dCompDens_dCompConc[c1][c2] * packPoreVolume;

      localAccumJacobian[c1+1][0] = -componentDensNew[c1] * volume;

    }

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
  using ElementViewConst = typename ElementRegionManager::ElementViewAccessor< VIEWTYPE >::ViewTypeConst;

  template< typename VIEWTYPE >
  using ElementView = typename ElementRegionManager::ElementViewAccessor< VIEWTYPE >::ViewType;

  /**
   * @brief The type for element-based constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::MaterialViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using MaterialView = typename ElementRegionManager::MaterialViewAccessor< VIEWTYPE >::ViewTypeConst;

  /**
   * @brief launches the kernel to assemble the flux contributions to the linear system.
   * @tparam STENCIL_TYPE The type of the stencil that is being used.
   */
  template< typename STENCIL_TYPE >
  static void
  Launch( STENCIL_TYPE const & stencil,
          localIndex const numDofPerCell,
          real64 const dt,
          ElementViewConst< arrayView1d< R1Tensor const > > const & transTMultiplier,
          integer const updateProppantPacking,
          R1Tensor const & unitGravityVector,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & dPres,
          ElementViewConst< arrayView1d< real64 const > > const & proppantConc,
          ElementViewConst< arrayView1d< real64 const > > const & dProppantConc,
          ElementViewConst< arrayView3d< real64 const > > const & componentDens,
          ElementViewConst< arrayView3d< real64 const > > const & dComponentDens_dPres,
          ElementViewConst< arrayView4d< real64 const > > const & dComponentDens_dComponentConc,
          ElementViewConst< arrayView1d< real64 const > > const & gravDepth,
          ElementViewConst< arrayView2d< real64 const > > const & dens,
          ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
          ElementViewConst< arrayView2d< real64 const > > const & dDens_dProppantConc,
          ElementViewConst< arrayView3d< real64 const > > const & dDens_dComponentConc,
          ElementViewConst< arrayView2d< real64 const > > const & visc,
          ElementViewConst< arrayView2d< real64 const > > const & dVisc_dPres,
          ElementViewConst< arrayView2d< real64 const > > const & dVisc_dProppantConc,
          ElementViewConst< arrayView3d< real64 const > > const & dVisc_dComponentConc,
          ElementViewConst< arrayView2d< real64 const > > const & fluidDensity,
          ElementViewConst< arrayView2d< real64 const > > const & dFluidDens_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & dFluidDens_dComponentConc,
          ElementViewConst< arrayView1d< real64 const > > const & settlingFactor,
          ElementViewConst< arrayView1d< real64 const > > const & dSettlingFactor_dPres,
          ElementViewConst< arrayView1d< real64 const > > const & dSettlingFactor_dProppantConc,
          ElementViewConst< arrayView2d< real64 const > > const & dSettlingFactor_dComponentConc,
          ElementViewConst< arrayView1d< real64 const > > const & collisionFactor,
          ElementViewConst< arrayView1d< real64 const > > const & dCollisionFactor_dProppantConc,
          ElementViewConst< arrayView1d< integer const > > const & isProppantMobile,
          ElementViewConst< arrayView1d< real64 const > > const & proppantPackVf,
          ElementViewConst< arrayView1d< real64 const > > const & aperture,
          ElementViewConst< arrayView1d< real64 const > > const & proppantLiftFlux,
          ParallelMatrix * const jacobian,
          ParallelVector * const residual );


  template< typename STENCIL_TYPE >
  static void
  LaunchCellBasedFluxCalculation( STENCIL_TYPE const & stencil,
                                  ElementViewConst< arrayView1d< R1Tensor const > > const & transTMultiplier,
                                  R1Tensor const unitGravityVector,
                                  ElementViewConst< arrayView1d< real64 const > > const & pres,
                                  ElementViewConst< arrayView1d< real64 const > > const & gravDepth,
                                  ElementViewConst< arrayView2d< real64 const > > const & dens,
                                  ElementViewConst< arrayView2d< real64 const > > const & visc,
                                  ElementViewConst< arrayView1d< real64 const > > const & aperture,
                                  ElementViewConst< arrayView1d< real64 const > > const & proppantPackVf,
                                  ElementView< arrayView1d< R1Tensor > > const & cellBasedFlux );

  /**
   * @brief Compute flux and its derivatives for a given multi-element connector.
   *
   * This is a specialized version that flux in a single region, and uses
   * element pairing instead of a proper junction.
   */
  inline static void
  ComputeJunction( localIndex const numElems,
                   localIndex const numDofPerCell,
                   arraySlice1d< localIndex const > const & stencilElementIndices,
                   arraySlice1d< real64 const > const & stencilWeights,
                   arraySlice1d< R1Tensor const > const & stencilCellCenterToEdgeCenters,
                   arrayView1d< real64 const > const & pres,
                   arrayView1d< real64 const > const & dPres,
                   arrayView1d< real64 const > const & proppantConc,
                   arrayView1d< real64 const > const & dProppantConc,
                   arrayView3d< real64 const > const & componentDens,
                   arrayView3d< real64 const > const & dComponentDens_dPres,
                   arrayView4d< real64 const > const & dComponentDens_dComponentConc,
                   arrayView1d< real64 const > const & gravDepth,
                   arrayView2d< real64 const > const & dens,
                   arrayView2d< real64 const > const & dDens_dPres,
                   arrayView2d< real64 const > const & dDens_dProppantConc,
                   arrayView3d< real64 const > const & dDens_dComponentConc,
                   arrayView2d< real64 const > const & visc,
                   arrayView2d< real64 const > const & dVisc_dPres,
                   arrayView2d< real64 const > const & dVisc_dProppantConc,
                   arrayView3d< real64 const > const & dVisc_dComponentConc,
                   arrayView2d< real64 const > const & fluidDensity,
                   arrayView2d< real64 const > const & dFluidDens_dPres,
                   arrayView3d< real64 const > const & dFluidDens_dComponentConc,
                   arrayView1d< real64 const > const & settlingFactor,
                   arrayView1d< real64 const > const & dSettlingFactor_dPres,
                   arrayView1d< real64 const > const & dSettlingFactor_dProppantConc,
                   arrayView2d< real64 const > const & dSettlingFactor_dComponentConc,
                   arrayView1d< real64 const > const & collisionFactor,
                   arrayView1d< real64 const > const & dCollisionFactor_dProppantConc,
                   arrayView1d< integer const > const & isProppantMobile,
                   arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( proppantPackVf ),
                   arrayView1d< real64 const > const & aperture,
                   arrayView1d< real64 const > const & proppantLiftFlux,
                   R1Tensor const & unitGravityVector,
                   arrayView1d< R1Tensor const > const & transTMultiplier,
                   real64 const dt,
                   arraySlice1d< real64 > const & localFlux,
                   arraySlice2d< real64 > const & localFluxJacobian )
  {

    // We assume numElems == stencilSize;

    static constexpr real64 TINY = 1e-10;

    // working array
    constexpr localIndex maxNumFluxElems = FaceElementStencil::NUM_POINT_IN_FLUX;
    constexpr localIndex maxNumComponents = ProppantTransport::MAX_NUM_COMPONENTS;

    localIndex const NC = numDofPerCell - 1;

    stackArray1d< real64, maxNumFluxElems > weight( numElems );

    // mixture density and fluid density in each face
    stackArray1d< real64, maxNumFluxElems > mixDens( numElems );
    stackArray1d< real64, maxNumFluxElems > dMixDens_dP( numElems );
    stackArray1d< real64, maxNumFluxElems > dMixDens_dProppantC( numElems );
    stackArray2d< real64, maxNumFluxElems * maxNumComponents > dMixDens_dComponentC( numElems, NC );

    stackArray1d< real64, maxNumFluxElems > fluidDens( numElems );
    stackArray1d< real64, maxNumFluxElems > dFluidDens_dP( numElems );
    stackArray2d< real64, maxNumFluxElems * maxNumComponents > dFluidDens_dComponentC( numElems, NC );


    // realted to slip velocity calculation
    stackArray1d< real64, maxNumFluxElems > settlingFac( numElems );
    stackArray1d< real64, maxNumFluxElems > dSettlingFac_dP( numElems );
    stackArray1d< real64, maxNumFluxElems > dSettlingFac_dProppantC( numElems );
    stackArray2d< real64, maxNumFluxElems * maxNumComponents > dSettlingFac_dComponentC( numElems, NC );

    stackArray1d< real64, maxNumFluxElems > collisionFac( numElems );
    stackArray1d< real64, maxNumFluxElems > dCollisionFac_dProppantC( numElems );

    stackArray1d< real64, maxNumFluxElems > transT( numElems );
    stackArray1d< real64, maxNumFluxElems > coefs( numElems );

    stackArray1d< integer, maxNumFluxElems > isProppantMob( numElems );

    real64 edgeDensity, edgeViscosity;
    stackArray1d< real64, maxNumFluxElems > dEdgeDens_dP( numElems );
    stackArray1d< real64, maxNumFluxElems > dEdgeDens_dProppantC( numElems );
    stackArray2d< real64, maxNumFluxElems * maxNumComponents > dEdgeDens_dComponentC( numElems, NC );

    stackArray1d< real64, maxNumFluxElems > dEdgeVisc_dP( numElems );
    stackArray1d< real64, maxNumFluxElems > dEdgeVisc_dProppantC( numElems );
    stackArray2d< real64, maxNumFluxElems * maxNumComponents > dEdgeVisc_dComponentC( numElems, NC );

    stackArray1d< real64, maxNumFluxElems > dPe_dP( numElems );
    stackArray1d< real64, maxNumFluxElems > dPe_dProppantC( numElems );
    stackArray2d< real64, maxNumFluxElems * maxNumComponents > dPe_dComponentC( numElems, NC );

    stackArray1d< R1Tensor, maxNumFluxElems > Up( numElems );

    stackArray1d< real64, maxNumFluxElems > dProppantCe_dP( numElems );
    stackArray1d< real64, maxNumFluxElems > dProppantCe_dProppantC( numElems );
    stackArray2d< real64, maxNumFluxElems * maxNumComponents > dProppantCe_dComponentC( numElems, NC );

    stackArray2d< real64, maxNumFluxElems * maxNumComponents > dComponentCe_dP( numElems, NC );
    stackArray2d< real64, maxNumFluxElems * maxNumComponents > dComponentCe_dProppantC( numElems, NC );
    stackArray3d< real64, maxNumFluxElems * maxNumComponents * maxNumComponents > dComponentCe_dComponentC( numElems, NC, NC );

    stackArray1d< real64, maxNumFluxElems > edgeToFaceFlux( numElems );
    stackArray2d< real64, maxNumFluxElems * maxNumFluxElems > dEdgeToFaceFlux_dP( numElems, numElems );
    stackArray2d< real64, maxNumFluxElems * maxNumFluxElems > dEdgeToFaceFlux_dProppantC( numElems, numElems );
    stackArray3d< real64, maxNumFluxElems * maxNumFluxElems * maxNumComponents > dEdgeToFaceFlux_dComponentC( numElems, numElems, NC );

    stackArray1d< real64, maxNumFluxElems > edgeToFaceProppantFlux( numElems );
    stackArray2d< real64, maxNumFluxElems * maxNumFluxElems > dEdgeToFaceProppantFlux_dP( numElems, numElems );
    stackArray2d< real64, maxNumFluxElems * maxNumFluxElems > dEdgeToFaceProppantFlux_dProppantC( numElems, numElems );
    stackArray3d< real64, maxNumFluxElems * maxNumFluxElems * maxNumComponents > dEdgeToFaceProppantFlux_dComponentC( numElems, numElems, NC );


    stackArray1d< real64, maxNumFluxElems > edgeToFaceFluidFlux( numElems );
    stackArray2d< real64, maxNumFluxElems * maxNumFluxElems > dEdgeToFaceFluidFlux_dP( numElems, numElems );
    stackArray2d< real64, maxNumFluxElems * maxNumFluxElems > dEdgeToFaceFluidFlux_dProppantC( numElems, numElems );
    stackArray3d< real64, maxNumFluxElems * maxNumFluxElems * maxNumComponents > dEdgeToFaceFluidFlux_dComponentC( numElems, numElems, NC );

    stackArray1d< real64, maxNumFluxElems > P( numElems );
    stackArray1d< real64, maxNumFluxElems > proppantC( numElems );
    stackArray2d< real64, maxNumFluxElems * maxNumComponents > componentC( numElems, NC );
    stackArray2d< real64, maxNumFluxElems * maxNumComponents > dComponentC_dP( numElems, NC );
    stackArray2d< real64, maxNumFluxElems * maxNumComponents > dComponentC_dProppantC( numElems, NC );
    stackArray3d< real64, maxNumFluxElems * maxNumComponents * maxNumComponents > dComponentC_dComponentC( numElems, NC, NC );

    // clear working arrays

    edgeDensity = 0.0;
    edgeViscosity = 0.0;

    dEdgeDens_dP = 0.0;
    dEdgeDens_dProppantC = 0.0;
    dEdgeDens_dComponentC = 0.0;

    dEdgeVisc_dP = 0.0;
    dEdgeVisc_dProppantC = 0.0;
    dEdgeVisc_dComponentC = 0.0;

    dEdgeToFaceProppantFlux_dP = 0.0;
    dEdgeToFaceProppantFlux_dProppantC = 0.0;
    dEdgeToFaceProppantFlux_dComponentC = 0.0;

    real64 aperTerm;
    real64 sumOfWeights = 0;

    real64 stencilEdgeToFaceDownDistance;
    real64 edgeLength;

    for( localIndex i=0; i<numElems; ++i )
    {

      localIndex const ei  = stencilElementIndices[i];


      aperTerm = aperture[ei] * aperture[ei] * aperture[ei];

      sumOfWeights += stencilWeights[i];
      weight[i] = stencilWeights[i];

      transT[i] = aperTerm * stencilWeights[i];

      edgeLength = 12.0 * stencilWeights[i] * stencilCellCenterToEdgeCenters[i].L2_Norm();

      stencilEdgeToFaceDownDistance = -Dot( stencilCellCenterToEdgeCenters[i], unitGravityVector ) * edgeLength / stencilCellCenterToEdgeCenters[i].L2_Norm();

      coefs[i] =  stencilEdgeToFaceDownDistance * aperture[ei];

      if( fabs( stencilEdgeToFaceDownDistance ) > TINY )
      {

        // vertical flow component

        transT[i] *= transTMultiplier[ei][1];

      }
      else
      {

        // horizontal flow component

        transT[i] *= transTMultiplier[ei][0];

      }

    }

    for( localIndex i=0; i<numElems; ++i )
      weight[i] /= sumOfWeights;

    //get averaged edgeDensity and edgeViscosity

    for( localIndex i = 0; i < numElems; ++i )
    {

      localIndex const ei  = stencilElementIndices[i];

      edgeDensity += weight[i] * dens[ei][0];
      dEdgeDens_dP[i] = weight[i] * dDens_dPres[ei][0];
      dEdgeDens_dProppantC[i] = weight[i] * dDens_dProppantConc[ei][0];

      edgeViscosity += weight[i] * visc[ei][0];
      dEdgeVisc_dP[i] = weight[i] * dVisc_dPres[ei][0];
      dEdgeVisc_dProppantC[i] = weight[i] * dVisc_dProppantConc[ei][0];

      P[i] = pres[ei] + dPres[ei];
      proppantC[i] = proppantConc[ei] + dProppantConc[ei];

      mixDens[i] = dens[ei][0];
      dMixDens_dP[i] = dDens_dPres[ei][0];
      dMixDens_dProppantC[i] = dDens_dProppantConc[ei][0];

      fluidDens[i] = fluidDensity[ei][0];
      dFluidDens_dP[i] = dFluidDens_dPres[ei][0];

      settlingFac[i] = settlingFactor[ei];
      dSettlingFac_dP[i] = dSettlingFactor_dPres[ei];
      dSettlingFac_dProppantC[i] = dSettlingFactor_dProppantConc[ei];

      collisionFac[i] = collisionFactor[ei];
      dCollisionFac_dProppantC[i] = dCollisionFactor_dProppantConc[ei];

      isProppantMob[i] = isProppantMobile[ei];

      for( localIndex c = 0; c < NC; ++c )
      {

        dMixDens_dComponentC[i][c] = dDens_dComponentConc[ei][0][c];
        dFluidDens_dComponentC[i][c] = dFluidDens_dComponentConc[ei][0][c];

        dSettlingFac_dComponentC[i][c] = dSettlingFactor_dComponentConc[ei][c];

        componentC[i][c] = componentDens[ei][0][c];
        dComponentC_dP[i][c] = dComponentDens_dPres[ei][0][c];

        dComponentC_dProppantC[i][c] = 0.0;

        for( localIndex c2 = 0; c2 < NC; ++c2 )
          dComponentC_dComponentC[i][c][c2] = dComponentDens_dComponentConc[ei][0][c][c2];


        dEdgeDens_dComponentC[i][c] = weight[i] * dDens_dComponentConc[ei][0][c];
        dEdgeVisc_dComponentC[i][c] = weight[i] * dVisc_dComponentConc[ei][0][c];

      }

    }

    real64 transTSum = 0.0;
    real64 Pe = 0.0;
    dPe_dP = 0.0;
    dPe_dComponentC = 0.0;
    dPe_dProppantC = 0.0;

    for( localIndex i = 0; i < numElems; ++i )
    {

      localIndex const ei  = stencilElementIndices[i];

      real64 const gravD    = gravDepth[ei];
      real64 const gravTerm = edgeDensity * gravD;

      Pe += transT[i] * (P[i] - gravTerm);

      transTSum += transT[i];

      dPe_dP[i] += transT[i];

      for( localIndex j = 0; j < numElems; ++j )
      {
        dPe_dP[j] += -transT[i] * gravD * dEdgeDens_dP[j];
        dPe_dProppantC[j] += -transT[i] * gravD * dEdgeDens_dProppantC[j];

        for( localIndex c = 0; c < NC; ++c )
        {
          dPe_dComponentC[j][c] += -transT[i] * gravD * dEdgeDens_dComponentC[j][c];
        }

      }

    }

    for( localIndex i = 0; i < numElems; ++i )
    {

      dPe_dP[i] /= transTSum;
      dPe_dProppantC[i] /= transTSum;

      for( localIndex c = 0; c < NC; ++c )
        dPe_dComponentC[i][c] /= transTSum;

    }

    Pe /= transTSum;

    for( localIndex i = 0; i < numElems; ++i )
    {

      localIndex const ei  = stencilElementIndices[i];

      real64 const gravD    = gravDepth[ei];
      real64 const gravTerm = edgeDensity * gravD;

      real64 const fluxTerm = Pe - (P[i] - gravTerm);

      edgeToFaceFlux[i] = transT[i] * fluxTerm / edgeViscosity;

      dEdgeToFaceFlux_dP[i][i] += -transT[i] / edgeViscosity;

      for( localIndex j = 0; j < numElems; ++j )
      {

        dEdgeToFaceFlux_dP[i][j] += -transT[i] * fluxTerm * dEdgeVisc_dP[j] / (edgeViscosity * edgeViscosity) + transT[i] *
                                    (dPe_dP[j] + dEdgeDens_dP[j] * gravD) / edgeViscosity;

        dEdgeToFaceFlux_dProppantC[i][j] += -transT[i] * fluxTerm * dEdgeVisc_dProppantC[j] / (edgeViscosity * edgeViscosity) + transT[i] *
                                            (dPe_dProppantC[j] + dEdgeDens_dProppantC[j] * gravD) / edgeViscosity;

        for( localIndex c = 0; c < NC; ++c )
        {

          dEdgeToFaceFlux_dComponentC[i][j][c] += -transT[i] * fluxTerm * dEdgeVisc_dComponentC[j][c] / (edgeViscosity * edgeViscosity) + transT[i] *
                                                  (dPe_dComponentC[j][c] + dEdgeDens_dComponentC[j][c] * gravD) / edgeViscosity;

        }

      }


      for( localIndex j = 0; j < numElems; ++j )
      {


        dEdgeToFaceProppantFlux_dProppantC[i][j] = 0.0;

        for( localIndex c = 0; c < NC; ++c )
        {

          dEdgeToFaceProppantFlux_dComponentC[i][j][c] = 0.0;

        }

      }


      if( fabs( coefs[i] ) > TINY )
      {

        // vertical

        edgeToFaceProppantFlux[i] = (1.0 - proppantC[i]) * settlingFac[i] * coefs[i] * fluidDens[i] / mixDens[i];

        dEdgeToFaceProppantFlux_dProppantC[i][i] = (-settlingFac[i] + (1 - proppantC[i]) * dSettlingFac_dProppantC[i]) * coefs[i] * fluidDens[i] / mixDens[i];


        edgeToFaceProppantFlux[i] += edgeToFaceFlux[i];

        for( localIndex j = 0; j < numElems; ++j )
        {

          dEdgeToFaceProppantFlux_dProppantC[i][j] += dEdgeToFaceFlux_dProppantC[i][j];

          for( localIndex c = 0; c < NC; ++c )
          {

            dEdgeToFaceProppantFlux_dComponentC[i][j][c] += dEdgeToFaceFlux_dComponentC[i][j][c];

          }
        }
      }
      else
      {

        // horizontal

        edgeToFaceProppantFlux[i] = (1.0 + fluidDens[i] / mixDens[i] * (1.0 - proppantC[i]) * collisionFac[i]) * edgeToFaceFlux[i];

        dEdgeToFaceProppantFlux_dProppantC[i][i] = -fluidDens[i] / mixDens[i] * (collisionFac[i] - (1.0 - proppantC[i]) * dCollisionFac_dProppantC[i]) *
                                                   edgeToFaceFlux[i];

        for( localIndex j = 0; j < numElems; ++j )
        {

          dEdgeToFaceProppantFlux_dProppantC[i][j] += (1.0 + fluidDens[i] / mixDens[i] * (1.0 - proppantC[i]) * collisionFac[i]) *
                                                      dEdgeToFaceFlux_dProppantC[i][j];

          for( localIndex c = 0; c < NC; ++c )
          {

            dEdgeToFaceProppantFlux_dComponentC[i][j][c] += (1.0 + fluidDens[i] / mixDens[i] * (1.0 - proppantC[i]) * collisionFac[i]) *
                                                            dEdgeToFaceFlux_dComponentC[i][j][c];

          }

        }

      }

      // fluid flux for component transport

      for( localIndex j = 0; j < numElems; ++j )
      {

        dEdgeToFaceFluidFlux_dP[i][j] = 0.0;
        dEdgeToFaceFluidFlux_dProppantC[i][j] = 0.0;

        for( localIndex c = 0; c < NC; ++c )
        {

          dEdgeToFaceFluidFlux_dComponentC[i][j][c] = 0.0;

        }

      }

      // note that all the fluid properties are from previous time step

      real64 fluidFluxCoef = 1.0;

      if( isProppantMob[i] == 0 || numElems == 1 )
        fluidFluxCoef = 0.0;

      edgeToFaceFluidFlux[i] = mixDens[i] / fluidDens[i] * edgeToFaceFlux[i] - fluidFluxCoef * (mixDens[i] -  fluidDens[i] * (1.0 - proppantC[i])) /
                               fluidDens[i] * edgeToFaceProppantFlux[i];

      dEdgeToFaceFluidFlux_dProppantC[i][i] = -fluidFluxCoef * edgeToFaceProppantFlux[i];

      for( localIndex j = 0; j < numElems; ++j )
      {

        dEdgeToFaceFluidFlux_dProppantC[i][j] += mixDens[i] / fluidDens[i] * dEdgeToFaceFlux_dProppantC[i][j] - fluidFluxCoef *
                                                 (mixDens[i] -  fluidDens[i] * (1.0 - proppantC[i])) / fluidDens[i] * dEdgeToFaceProppantFlux_dProppantC[i][j];

        for( localIndex c = 0; c < NC; ++c )
        {

          dEdgeToFaceFluidFlux_dComponentC[i][j][c] += mixDens[i] / fluidDens[i] * dEdgeToFaceFlux_dComponentC[i][j][c] - fluidFluxCoef *
                                                       (mixDens[i] -  fluidDens[i] * (1.0 - proppantC[i])) / fluidDens[i] *
                                                       dEdgeToFaceProppantFlux_dComponentC[i][j][c];

        }
      }


      for( localIndex j = 0; j < numElems; ++j )
      {

        dEdgeToFaceFluidFlux_dProppantC[i][j] = 0.0;

        for( localIndex c = 0; c < NC; ++c )
        {

          dEdgeToFaceFluidFlux_dComponentC[i][j][c] = 0.0;

        }
      }
    }

    // get proppantCe

    real64 proppantCe = 0.0;
    dProppantCe_dP = 0.0;
    dProppantCe_dProppantC = 0.0;
    dProppantCe_dComponentC = 0.0;

    real64 downStreamFlux = 0.0;
    stackArray1d< real64, maxNumFluxElems > dDownStreamFlux_dP( numElems );
    stackArray1d< real64, maxNumFluxElems > dDownStreamFlux_dProppantC( numElems );
    stackArray2d< real64, maxNumFluxElems * maxNumComponents > dDownStreamFlux_dComponentC( numElems, NC );

    dDownStreamFlux_dP = 0.0;
    dDownStreamFlux_dProppantC = 0.0;
    dDownStreamFlux_dComponentC = 0.0;

    for( localIndex i = 0; i < numElems; ++i )
    {

      if( isProppantMob[i] == 0 )
        continue;

      if( edgeToFaceProppantFlux[i] >= 0.0 )
      {
        // downstream
        downStreamFlux += edgeToFaceProppantFlux[i];

        for( localIndex j = 0; j < numElems; ++j )
        {

          dDownStreamFlux_dP[j] += dEdgeToFaceProppantFlux_dP[i][j];
          dDownStreamFlux_dProppantC[j] += dEdgeToFaceProppantFlux_dProppantC[i][j];

          for( localIndex c = 0; c < NC; ++c )
          {

            dDownStreamFlux_dComponentC[j][c] += dEdgeToFaceProppantFlux_dComponentC[i][j][c];

          }

        }

      }
      else
      {

        // upstream
        proppantCe += -edgeToFaceProppantFlux[i] * proppantC[i];

        dProppantCe_dProppantC[i] += -edgeToFaceProppantFlux[i];

        for( localIndex j = 0; j < numElems; ++j )
        {

          dProppantCe_dP[j] += -dEdgeToFaceProppantFlux_dP[i][j] * proppantC[i];
          dProppantCe_dProppantC[j] += -dEdgeToFaceProppantFlux_dProppantC[i][j] * proppantC[i];

          for( localIndex c = 0; c < NC; ++c )
          {

            dProppantCe_dComponentC[j][c] += -dEdgeToFaceProppantFlux_dComponentC[i][j][c] * proppantC[i];

          }

        }

      }

    }


    if( downStreamFlux > 0.0 )
    {

      for( localIndex i = 0; i < numElems; ++i )
      {

        if( isProppantMob[i] == 0 )
          continue;

        dProppantCe_dP[i] =  dProppantCe_dP[i] / downStreamFlux - proppantCe * dDownStreamFlux_dP[i] / (downStreamFlux * downStreamFlux);
        dProppantCe_dProppantC[i] =  dProppantCe_dProppantC[i] / downStreamFlux - proppantCe * dDownStreamFlux_dProppantC[i] /
                                    (downStreamFlux * downStreamFlux);;

        for( localIndex c = 0; c < NC; ++c )
        {

          dProppantCe_dComponentC[i][c] =  dProppantCe_dComponentC[i][c] / downStreamFlux - proppantCe * dDownStreamFlux_dComponentC[i][c] /
                                          (downStreamFlux * downStreamFlux);

        }
      }

      proppantCe = proppantCe / downStreamFlux;

    }
    else
    {

      proppantCe = 0.0;
      for( localIndex i = 0; i < numElems; ++i )
      {

        if( isProppantMob[i] == 0 )
          continue;

        dProppantCe_dP[i] =  0.0;
        dProppantCe_dProppantC[i] =  weight[i];

        proppantCe += proppantC[i] * weight[i];

      }

    }


    // get componentCe

    array1d< real64 > componentCe( NC );
    componentCe = 0.0;
    dComponentCe_dP = 0.0;
    dComponentCe_dProppantC = 0.0;
    dComponentCe_dComponentC = 0.0;

    downStreamFlux = 0.0;
    dDownStreamFlux_dP = 0.0;
    dDownStreamFlux_dProppantC = 0.0;
    dDownStreamFlux_dComponentC = 0.0;

    for( localIndex i = 0; i < numElems; ++i )
    {

      if( edgeToFaceFluidFlux[i] >= 0.0 )
      {
        // downstream

        downStreamFlux += edgeToFaceFluidFlux[i];


        for( localIndex j = 0; j < numElems; ++j )
        {

          dDownStreamFlux_dP[j] += dEdgeToFaceFluidFlux_dP[i][j];
          dDownStreamFlux_dProppantC[j] += dEdgeToFaceFluidFlux_dProppantC[i][j];

          for( localIndex c = 0; c < NC; ++c )
          {

            dDownStreamFlux_dComponentC[j][c] += dEdgeToFaceFluidFlux_dComponentC[i][j][c];

          }

        }

      }
      else
      {

        // upstream
        for( localIndex c1 = 0; c1 < NC; ++c1 )
        {

          componentCe[c1] += -edgeToFaceFluidFlux[i] * componentC[i][c1];

          dComponentCe_dP[i][c1] += -edgeToFaceFluidFlux[i] * dComponentC_dP[i][c1];

          for( localIndex c2 = 0; c2 < NC; ++c2 )
            dComponentCe_dComponentC[i][c1][c2] += -edgeToFaceFluidFlux[i] * dComponentC_dComponentC[i][c1][c2];

          for( localIndex j = 0; j < numElems; ++j )
          {

            dComponentCe_dP[j][c1] += -dEdgeToFaceFluidFlux_dP[i][j] * componentC[i][c1];
            dComponentCe_dProppantC[j][c1] += -dEdgeToFaceFluidFlux_dProppantC[i][j] * componentC[i][c1];

            for( localIndex c2 = 0; c2 < NC; ++c2 )
            {

              dComponentCe_dComponentC[j][c1][c2] += -dEdgeToFaceFluidFlux_dComponentC[i][j][c2] * componentC[i][c1];

            }
          }
        }

      }

    }


    if( downStreamFlux > 0.0 )
    {

      for( localIndex c1 = 0; c1 < NC; ++c1 )
      {

        for( localIndex i = 0; i < numElems; ++i )
        {

          dComponentCe_dP[i][c1] =  dComponentCe_dP[i][c1] / downStreamFlux - componentCe[c1] * dDownStreamFlux_dP[i] / (downStreamFlux * downStreamFlux);
          dComponentCe_dProppantC[i][c1] =  dComponentCe_dProppantC[i][c1] / downStreamFlux - componentCe[c1] * dDownStreamFlux_dProppantC[i] /
                                           (downStreamFlux * downStreamFlux);

          for( localIndex c2 = 0; c2 < NC; ++c2 )
          {

            dComponentCe_dComponentC[i][c1][c2] =  dComponentCe_dComponentC[i][c1][c2] / downStreamFlux - componentCe[c1] * dDownStreamFlux_dComponentC[i][c2] /
                                                  (downStreamFlux * downStreamFlux);

          }

        }

        componentCe[c1] = componentCe[c1] / downStreamFlux;

      }

    }
    else
    {

      for( localIndex c = 0; c < NC; ++c )
      {

        componentCe[c] = 0.0;

        for( localIndex i = 0; i < numElems; ++i )
        {

          componentCe[c] += componentC[i][c] * weight[i];

          dComponentCe_dP[i][c] = dComponentC_dP[i][c] * weight[i];
          dComponentCe_dProppantC[i][c] = dComponentC_dProppantC[i][c] * weight[i];

          for( localIndex c2 = 0; c2 < NC; ++c2 )
            dComponentCe_dComponentC[i][c][c2] =  dComponentC_dComponentC[i][c][c2] * weight[i];

        }

      }

    }


    for( localIndex i = 0; i < numElems; ++i )
    {

      localIndex const ei  = stencilElementIndices[i];

      localIndex idx1 = i * numDofPerCell; // proppant

      if( isProppantMob[i] == 1 )
      {

        if( edgeToFaceProppantFlux[i] >= 0.0 )
        {


          localFlux[idx1] = -proppantCe * edgeToFaceProppantFlux[i] * dt- proppantLiftFlux[ei] * dt;

        }
        else
        {

          localFlux[idx1] = -proppantC[i] * edgeToFaceProppantFlux[i] * dt- proppantLiftFlux[ei] * dt;

        }

        for( localIndex j = 0; j < numElems; ++j )
        {

          localIndex idx2 = j * numDofPerCell;


          if( edgeToFaceProppantFlux[i] >= 0.0 )
          {

            localFluxJacobian[idx1][idx2] = -(dProppantCe_dProppantC[j] * edgeToFaceProppantFlux[i] + proppantCe * dEdgeToFaceProppantFlux_dProppantC[i][j]) *
                                            dt;

          }
          else
          {

            localFluxJacobian[idx1][idx2] = -proppantC[i] * dEdgeToFaceProppantFlux_dProppantC[i][j] * dt;

            if( i == j )
              localFluxJacobian[idx1][idx2] += -edgeToFaceProppantFlux[i] * dt;
          }

        }

      }
      else
      {

        localFlux[idx1] = 0.0;

        for( localIndex j = 0; j < numElems; ++j )
        {

          for( localIndex c = 0; c < numDofPerCell; ++c )
          {

            localIndex idx2 = j * numDofPerCell + c;
            localFluxJacobian[idx1][idx2] = 0.0;

          }

        }

      }


      // component

      if( numElems > 1 )
      {

        for( localIndex c1 = 0; c1 < NC; ++c1 )
        {

          idx1 = i * numDofPerCell + 1 + c1;

          if( edgeToFaceFluidFlux[i] >= 0.0 )
          {

            localFlux[idx1] = -componentCe[c1] * edgeToFaceFluidFlux[i] * dt;

          }
          else
          {

            localFlux[idx1] = -componentC[i][c1] * edgeToFaceFluidFlux[i] * dt;

          }

          for( localIndex j = 0; j < numElems; ++j )
          {

            localIndex idx2 = j * numDofPerCell;

            if( edgeToFaceFluidFlux[i] >= 0.0 )
            {

              for( localIndex c2 = 0; c2 < NC; ++c2 )
              {

                localFluxJacobian[idx1][idx2 + 1 +
                                        c2] =
                  -(dComponentCe_dComponentC[j][c1][c2] * edgeToFaceFluidFlux[i] + 0*componentCe[c1] * dEdgeToFaceFluidFlux_dComponentC[i][j][c2]) * dt;

              }

            }
            else
            {


              if( i == j )
              {

                for( localIndex c2 = 0; c2 < NC; ++c2 )
                {

                  localFluxJacobian[idx1][idx2 + 1 + c2] += -dComponentC_dComponentC[i][c1][c2] * edgeToFaceFluidFlux[i] * dt;

                }

              }

            }

          }

        }

      }

    }

  }

  inline static void
  ComputeCellBasedFlux( localIndex const numElems,
                        arraySlice1d< localIndex const > const & stencilElementIndices,
                        arraySlice1d< real64 const > const & stencilWeights,
                        arraySlice1d< R1Tensor const > const & stencilCellCenterToEdgeCenters,
                        arrayView1d< R1Tensor const > const & transMultiplier,
                        R1Tensor const unitGravityVector,
                        arrayView1d< real64 const > const & pres,
                        arrayView1d< real64 const > const & gravDepth,
                        arrayView2d< real64 const > const & dens,
                        arrayView2d< real64 const > const & visc,
                        arrayView1d< real64 const > const & aperture,
                        arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( proppantPackVf ),
                        arrayView1d< R1Tensor > const & cellBasedFlux )
  {

    static constexpr real64 TINY = 1e-10;

    constexpr localIndex maxNumFluxElems = FaceElementStencil::NUM_POINT_IN_FLUX;

    stackArray1d< real64, maxNumFluxElems > weight( numElems );

    stackArray1d< real64, maxNumFluxElems > transT( numElems );

    real64 edgeToFaceFlux;

    // clear working arrays

    real64 edgeDensity = 0.0;
    real64 edgeViscosity = 0.0;

    real64 aperTerm;
    real64 sumOfWeights = 0;

    real64 stencilEdgeToFaceDownDistance;
    real64 edgeLength;

    for( localIndex i=0; i<numElems; ++i )
    {

      localIndex const ei  = stencilElementIndices[i];

      aperTerm = aperture[ei] * aperture[ei] * aperture[ei];

      sumOfWeights += stencilWeights[i];
      weight[i] = stencilWeights[i];

      transT[i] = aperTerm * stencilWeights[i];

      edgeLength = 12.0 * stencilWeights[i] * stencilCellCenterToEdgeCenters[i].L2_Norm();

      stencilEdgeToFaceDownDistance = -Dot( stencilCellCenterToEdgeCenters[i], unitGravityVector ) * edgeLength / stencilCellCenterToEdgeCenters[i].L2_Norm();

      if( fabs( stencilEdgeToFaceDownDistance ) > TINY )
      {

        // vertical flow component

        transT[i] *= transMultiplier[ei][1];

      }
      else
      {

        // horizontal flow component

        transT[i] *= transMultiplier[ei][0];

      }

    }

    for( localIndex i=0; i<numElems; ++i )
      weight[i] /= sumOfWeights;

    //get averaged edgeDensity and edgeViscosity

    for( localIndex i = 0; i < numElems; ++i )
    {

      localIndex const ei  = stencilElementIndices[i];

      edgeDensity += weight[i] * dens[ei][0];

      edgeViscosity += weight[i] * visc[ei][0];

    }

    real64 transTSum = 0.0;
    real64 Pe = 0.0;

    for( localIndex i = 0; i < numElems; ++i )
    {

      localIndex const ei  = stencilElementIndices[i];

      real64 const gravD    = gravDepth[ei];
      real64 const gravTerm = edgeDensity * gravD;

      Pe += transT[i] * (pres[ei] - gravTerm);

      transTSum += transT[i];

    }

    Pe /= transTSum;

    for( localIndex i = 0; i < numElems; ++i )
    {

      localIndex const ei  = stencilElementIndices[i];

      real64 const gravD    = gravDepth[ei];
      real64 const gravTerm = edgeDensity * gravD;

      real64 const fluxTerm = Pe - (pres[ei] - gravTerm);

      edgeToFaceFlux = transT[i] * fluxTerm / edgeViscosity;

      for( localIndex k = 0; k < 3; ++k )
        cellBasedFlux[ei][k] += -edgeToFaceFlux * stencilCellCenterToEdgeCenters[i][k];

    }
  }

};

struct ProppantPackVolumeKernel
{

  template< typename VIEWTYPE >
  using ElementView = typename ElementRegionManager::ElementViewAccessor< VIEWTYPE >::ViewType;

  template< typename VIEWTYPE >
  using ElementViewConst = typename ElementRegionManager::ElementViewAccessor< VIEWTYPE >::ViewTypeConst;

  template< typename VIEWTYPE >
  using MaterialView = typename ElementRegionManager::MaterialViewAccessor< VIEWTYPE >::ViewType;

  template< typename VIEWTYPE >
  using MaterialViewConst = typename ElementRegionManager::MaterialViewAccessor< VIEWTYPE >::ViewTypeConst;

  template< typename STENCIL_TYPE >
  static void
  LaunchProppantPackVolumeCalculation( STENCIL_TYPE const & stencil,
                                       real64 const dt,
                                       real64 const proppantDensity,
                                       real64 const proppantDiameter,
                                       real64 const maxProppantConcentration,
                                       R1Tensor const unitGravityVector,
                                       real64 const criticalShieldsNumber,
                                       real64 const fricitonCoefficient,
                                       ElementView< arrayView1d< real64 > > const & conc,
                                       ElementView< arrayView1d< real64 const > > const & settlingFactor,
                                       ElementView< arrayView2d< real64 const > > const & density,
                                       ElementView< arrayView2d< real64 const > > const & fluidDensity,
                                       ElementView< arrayView2d< real64 const > > const & fluidViscosity,
                                       ElementView< arrayView1d< integer > > const & isProppantMobile,
                                       ElementView< arrayView1d< integer const > > const & isProppantBoundaryElement,
                                       ElementView< arrayView1d< real64 > > const & proppantPackVf,
                                       ElementView< arrayView1d< real64 > > const & proppantExcessPackV,
                                       ElementView< arrayView1d< real64 const > > const & aperture,
                                       ElementView< arrayView1d< real64 const > > const & volume,
                                       ElementView< arrayView1d< integer const > > const & elemGhostRank,
                                       ElementView< arrayView1d< R1Tensor const > > const & cellBasedFlux,
                                       ElementView< arrayView1d< real64 > > const & proppantLiftFlux );


  template< typename STENCIL_TYPE >
  static void
  LaunchProppantPackVolumeUpdate( STENCIL_TYPE const & stencil,
                                  R1Tensor const unitGravityVector,
                                  real64 const maxProppantConcentration,
                                  ElementView< arrayView1d< real64 > > const & conc,
                                  ElementView< arrayView1d< integer > > const & isProppantMobile,
                                  ElementView< arrayView1d< real64 > > const & proppantPackVf,
                                  ElementView< arrayView1d< real64 > > const & proppantExcessPackV );


  inline static void
  ComputeProppantPackVolume( localIndex const numElems,
                             real64 const dt,
                             real64 const proppantDensity,
                             real64 const proppantDiameter,
                             real64 const maxProppantConcentration,
                             R1Tensor const unitGravityVector,
                             real64 const criticalShieldsNumber,
                             real64 const frictionCoefficient,
                             arraySlice1d< localIndex const > const & stencilElementIndices,
                             arraySlice1d< real64 const > const & stencilWeights,
                             arraySlice1d< R1Tensor const > const & stencilCellCenterToEdgeCenters,
                             arrayView1d< real64 const > const & settlingFactor,
                             arrayView2d< real64 const > const & density,
                             arrayView2d< real64 const > const & fluidDensity,
                             //                             arrayView2d<real64 const> const & fluidViscosity,
                             arrayView2d< real64 const > const &,
                             arrayView1d< real64 const > const & volume,
                             arrayView1d< real64 const > const & aperture,
                             arrayView1d< integer const > const & elemGhostRank,
                             arrayView1d< integer const > const & isProppantBoundaryElement,
                             arrayView1d< real64 > const & conc,
                             arrayView1d< integer > const & isProppantMobile,
                             arrayView1d< real64 > const & proppantPackVf,
                             arrayView1d< real64 > const & proppantExcessPackV,
                             arrayView1d< R1Tensor const > const & cellBasedFlux,
                             arrayView1d< real64 > const & proppantLiftFlux )
  {

    integer faceIndex = -1;

    static constexpr real64 TINY = 1e-10;

    real64 edgeLength = 12.0 * stencilWeights[0] * stencilCellCenterToEdgeCenters[0].L2_Norm();

    if( numElems == 1 )
    {

      localIndex const ei  = stencilElementIndices[0];

      real64 stencilEdgeToFaceDownDistance =
        -Dot( stencilCellCenterToEdgeCenters[0], unitGravityVector ) * edgeLength / stencilCellCenterToEdgeCenters[0].L2_Norm();

      if( stencilEdgeToFaceDownDistance < -TINY && isProppantMobile[ei] == 1 )
        // bottom face element
        faceIndex = 0;

    }
    else if( numElems == 2 )
    {

      localIndex const ei0  = stencilElementIndices[0];
      localIndex const ei1  = stencilElementIndices[1];

      real64 stencilEdgeToFaceDownDistance0 =
        -Dot( stencilCellCenterToEdgeCenters[0], unitGravityVector ) * edgeLength / stencilCellCenterToEdgeCenters[0].L2_Norm();

      real64 stencilEdgeToFaceDownDistance1 =
        -Dot( stencilCellCenterToEdgeCenters[1], unitGravityVector ) * edgeLength / stencilCellCenterToEdgeCenters[1].L2_Norm();

      //0: top  1: bottom

      if( stencilEdgeToFaceDownDistance0 < -TINY && stencilEdgeToFaceDownDistance1 > TINY && isProppantMobile[ei0] == 1 && isProppantMobile[ei1] == 0 )
      {

        faceIndex = 0;

      }

      //0: bottom  1: top

      if( stencilEdgeToFaceDownDistance0 > TINY && stencilEdgeToFaceDownDistance1 < -TINY && isProppantMobile[ei1] == 1 && isProppantMobile[ei0] == 0 )
      {

        faceIndex = 1;

      }
    }


    if( faceIndex >= 0 )
    {

      localIndex const ei = stencilElementIndices[faceIndex];

      if( elemGhostRank[ei] < 0 && isProppantBoundaryElement[ei] == 0 )
      {

        real64 L = stencilCellCenterToEdgeCenters[faceIndex].L2_Norm() * 2.0;

        R1Tensor velocity = cellBasedFlux[ei];

        real64 downVelocity = Dot( cellBasedFlux[ei], unitGravityVector );

        for( localIndex idx = 0; idx < 3; ++idx )
          velocity[idx] -= downVelocity * unitGravityVector[idx];

        real64 velocityMag = velocity.L2_Norm() / volume[ei];

        real64 dH = fluidDensity[ei][0] / density[ei][0] * (1.0 - conc[ei]) * settlingFactor[ei] * conc[ei] / maxProppantConcentration * dt;

        real64 tau = 1.0/8.0 * frictionCoefficient * fluidDensity[ei][0] * velocityMag * velocityMag;

        real64 ShieldsNumber = tau / (proppantDensity - fluidDensity[ei][0]) / 9.81 / proppantDiameter;

        proppantLiftFlux[ei] = 0.0;

        if( ShieldsNumber > criticalShieldsNumber )
        {

          proppantLiftFlux[ei] = aperture[ei] *
                                 (proppantDiameter * sqrt( 9.81 * proppantDiameter * (proppantDensity - fluidDensity[ei][0]) / fluidDensity[ei][0] )) *
                                 (9.64 * pow( ShieldsNumber, 0.166 )) * pow( ShieldsNumber - criticalShieldsNumber, 1.5 );

        }

        if( proppantLiftFlux[ei] < 0.0 )
          proppantLiftFlux[ei] = 0.0;

        dH -=  proppantLiftFlux[ei] / edgeLength / aperture[ei] / maxProppantConcentration * dt;

        real64 Vf = proppantPackVf[ei];

        Vf += dH / L;

        if( Vf >= 1.0 )
        {

          isProppantMobile[ei] = 0;

          proppantExcessPackV[ei] = (Vf - 1.0) * L;

          proppantPackVf[ei] = 1.0;

          conc[ei] = maxProppantConcentration;

        }
        else
        {

          proppantPackVf[ei] = Vf;

          if( proppantPackVf[ei] < 0.0 )
            proppantPackVf[ei] = 0.0;

        }

      }

    }
  }


  inline static void
  UpdateProppantPackVolume( localIndex const numElems,
                            arraySlice1d< localIndex const > const & stencilElementIndices,
                            arraySlice1d< real64 const > const & stencilWeights,
                            arraySlice1d< R1Tensor const > const & stencilCellCenterToEdgeCenters,
                            R1Tensor const unitGravityVector,
                            real64 const maxProppantConcentration,
                            arrayView1d< real64 > const & conc,
                            arrayView1d< integer > const & isProppantMobile,
                            arrayView1d< real64 > const & proppantPackVf,
                            arrayView1d< real64 > const & proppantExcessPackV )
  {

    integer faceIndex = -1;
    real64 excessV = 0.0;

    static constexpr real64 TINY = 1e-10;

    real64 edgeLength = 12.0 * stencilWeights[0] * stencilCellCenterToEdgeCenters[0].L2_Norm();

    if( numElems == 2 )
    {

      localIndex const ei0  = stencilElementIndices[0];
      localIndex const ei1  = stencilElementIndices[1];

      real64 stencilEdgeToFaceDownDistance0 =
        -Dot( stencilCellCenterToEdgeCenters[0], unitGravityVector ) * edgeLength / stencilCellCenterToEdgeCenters[0].L2_Norm();

      real64 stencilEdgeToFaceDownDistance1 =
        -Dot( stencilCellCenterToEdgeCenters[1], unitGravityVector ) * edgeLength / stencilCellCenterToEdgeCenters[1].L2_Norm();

      //0: top  1: bottom

      if( stencilEdgeToFaceDownDistance0 < -TINY && stencilEdgeToFaceDownDistance1 > TINY && isProppantMobile[ei0] == 1 && isProppantMobile[ei1] == 0 &&
          proppantExcessPackV[ei1] > 0.0 )
      {

        faceIndex = 0;
        excessV = proppantExcessPackV[ei1];

      }

      //0: bottom  1: top

      if( stencilEdgeToFaceDownDistance0 > TINY && stencilEdgeToFaceDownDistance1 < -TINY && isProppantMobile[ei1] == 1 && isProppantMobile[ei0] == 0 &&
          proppantExcessPackV[ei0] > 0.0 )
      {

        faceIndex = 1;
        excessV = proppantExcessPackV[ei0];

      }
    }


    if( faceIndex >= 0 )
    {

      localIndex const ei = stencilElementIndices[faceIndex];

      real64 L = stencilCellCenterToEdgeCenters[faceIndex].L2_Norm() * 2.0;

      real64 Vf = proppantPackVf[ei];

      Vf += excessV / L;

      if( Vf >= 1.0 )
      {

        isProppantMobile[ei] = 0;

        proppantPackVf[ei] = 1.0;

        conc[ei] = maxProppantConcentration;

      }
      else
      {

        proppantPackVf[ei] = Vf;


      }

    }

  }

};


} // namespace ProppantTransportKernels



} // namespace geosx

#endif //GEOSX_PROPPANTTRANSPORTKERNELS_HPP
