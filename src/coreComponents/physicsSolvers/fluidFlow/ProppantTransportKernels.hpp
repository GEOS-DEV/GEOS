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
 * @file ProppantTransportKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_PROPPANTTRANSPORTKERNELS_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_PROPPANTTRANSPORTKERNELS_HPP_

#include "common/DataTypes.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "physicsSolvers/fluidFlow/ProppantTransport.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace ProppantTransportKernels
{

/******************************** FluidUpdateKernel ********************************/

struct FluidUpdateKernel
{
  template< typename FLUID_WRAPPER >
  static void Launch( FLUID_WRAPPER const & fluidWrapper,
                      arrayView1d< real64 const > const & pres,
                      arrayView1d< real64 const > const & dPres,
                      arrayView2d< real64 const > const & componentConcentration,
                      arrayView2d< real64 const > const & dComponentConcentration )
  {
    forAll< parallelDevicePolicy<> >( fluidWrapper.numElems(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const NC = fluidWrapper.numFluidComponents();
      stackArray1d< real64, constitutive::SlurryFluidBase::MAX_NUM_COMPONENTS > compConc( NC );

      for( localIndex c = 0; c < NC; ++c )
      {
        compConc[c] = componentConcentration[a][c] + dComponentConcentration[a][c];
      }

      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.UpdateFluidProperty( a, q,
                                          pres[a] + dPres[a],
                                          compConc,
                                          0.0 );
      }
    } );
  }
};

/******************************** ComponentDensityUpdateKernel ********************************/

struct ComponentDensityUpdateKernel
{
  template< typename FLUID_WRAPPER >
  static void Launch( FLUID_WRAPPER const & fluidWrapper,
                      arrayView1d< real64 const > const & pres,
                      arrayView1d< real64 const > const & dPres,
                      arrayView2d< real64 const > const & componentConcentration,
                      arrayView2d< real64 const > const & dComponentConcentration )
  {
    forAll< parallelDevicePolicy<> >( fluidWrapper.numElems(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const NC = fluidWrapper.numFluidComponents();
      stackArray1d< real64, constitutive::SlurryFluidBase::MAX_NUM_COMPONENTS > compConc( NC );

      for( localIndex c = 0; c < NC; ++c )
      {
        compConc[c] = componentConcentration[a][c] + dComponentConcentration[a][c];
      }

      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.UpdateComponentDensity( a, q,
                                             pres[a] + dPres[a],
                                             compConc );
      }
    } );
  }
};

/******************************** ProppantUpdateKernel ********************************/

struct ProppantUpdateKernel
{
  template< typename PROPPANT_WRAPPER >
  static void Launch( PROPPANT_WRAPPER const & proppantWrapper,
                      arrayView1d< real64 const > const & proppantConc,
                      arrayView1d< real64 const > const & dProppantConc,
                      arrayView2d< real64 const > const & fluidDens,
                      arrayView2d< real64 const > const & dFluidDensDPres,
                      arrayView3d< real64 const > const & dFluidDensDCompConc,
                      arrayView2d< real64 const > const & fluidVisc,
                      arrayView2d< real64 const > const & dFluidViscDPres,
                      arrayView3d< real64 const > const & dFluidViscDCompConc )
  {
    forAll< parallelDevicePolicy<> >( proppantWrapper.numElems(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      proppantWrapper.Update( a,
                              proppantConc[a] + dProppantConc[a],
                              fluidDens[a][0],
                              dFluidDensDPres[a][0],
                              dFluidDensDCompConc[a][0],
                              fluidVisc[a][0],
                              dFluidViscDPres[a][0],
                              dFluidViscDCompConc[a][0] );
    } );
  }
};

/******************************** AccumulationKernel ********************************/

struct AccumulationKernel
{
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  Compute( localIndex const nc,
           real64 const proppantConcOld,
           real64 const proppantConcNew,
           arraySlice1d< real64 const > const & componentDensOld,
           arraySlice1d< real64 const > const & componentDensNew,
           arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM( dCompDens_dPres ),
           arraySlice2d< real64 const > const & dCompDensDCompConc,
           real64 const volume,
           real64 const packPoreVolume,
           real64 const proppantLiftVolume,
           arraySlice1d< real64 > const & localAccum,
           arraySlice2d< real64 > const & localAccumJacobian );

  static void
  Launch( localIndex const size,
          localIndex const nc,
          localIndex const ndof,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & elemGhostRank,
          arrayView1d< real64 const > const & proppantConc,
          arrayView1d< real64 const > const & dProppantConc,
          arrayView2d< real64 const > const & componentDensOld,
          arrayView3d< real64 const > const & componentDens,
          arrayView3d< real64 const > const & dCompDensDPres,
          arrayView4d< real64 const > const & dCompDensDCompConc,
          arrayView1d< real64 const > const & volume,
          arrayView1d< real64 const > const & proppantPackVf,
          arrayView1d< real64 const > const & proppantLiftFlux,
          real64 const dt,
          real64 const maxProppantConcentration,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );
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
  using ElementViewConst = typename ElementRegionManager::ElementViewConst< VIEWTYPE >;

  template< typename VIEWTYPE >
  using ElementView = typename ElementRegionManager::ElementView< VIEWTYPE >;

  /**
   * @brief launches the kernel to assemble the flux contributions to the linear system.
   * @tparam STENCIL_TYPE The type of the stencil that is being used.
   */
  template< typename STENCIL_TYPE >
  static void
  Launch( STENCIL_TYPE const & stencil,
          localIndex const numDofPerCell,
          real64 const dt,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< R1Tensor const > > const & transTMultiplier,
          integer const updateProppantPacking,
          R1Tensor const & unitGravityVector,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & dPres,
          ElementViewConst< arrayView1d< real64 const > > const & proppantConc,
          ElementViewConst< arrayView1d< real64 const > > const & dProppantConc,
          ElementViewConst< arrayView3d< real64 const > > const & componentDens,
          ElementViewConst< arrayView3d< real64 const > > const & dComponentDensDPres,
          ElementViewConst< arrayView4d< real64 const > > const & dComponentDensDComponentConc,
          ElementViewConst< arrayView1d< real64 const > > const & gravDepth,
          ElementViewConst< arrayView2d< real64 const > > const & dens,
          ElementViewConst< arrayView2d< real64 const > > const & dDensDPres,
          ElementViewConst< arrayView2d< real64 const > > const & dDensDProppantConc,
          ElementViewConst< arrayView3d< real64 const > > const & dDensDComponentConc,
          ElementViewConst< arrayView2d< real64 const > > const & visc,
          ElementViewConst< arrayView2d< real64 const > > const & dViscDPres,
          ElementViewConst< arrayView2d< real64 const > > const & dViscDProppantConc,
          ElementViewConst< arrayView3d< real64 const > > const & dViscDComponentConc,
          ElementViewConst< arrayView2d< real64 const > > const & fluidDensity,
          ElementViewConst< arrayView2d< real64 const > > const & dFluidDensDPres,
          ElementViewConst< arrayView3d< real64 const > > const & dFluidDensDComponentConc,
          ElementViewConst< arrayView1d< real64 const > > const & settlingFactor,
          ElementViewConst< arrayView1d< real64 const > > const & dSettlingFactorDPres,
          ElementViewConst< arrayView1d< real64 const > > const & dSettlingFactorDProppantConc,
          ElementViewConst< arrayView2d< real64 const > > const & dSettlingFactorDComponentConc,
          ElementViewConst< arrayView1d< real64 const > > const & collisionFactor,
          ElementViewConst< arrayView1d< real64 const > > const & dCollisionFactorDProppantConc,
          ElementViewConst< arrayView1d< integer const > > const & isProppantMobile,
          ElementViewConst< arrayView1d< real64 const > > const & proppantPackVf,
          ElementViewConst< arrayView1d< real64 const > > const & aperture,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );


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
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
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
                   arrayView3d< real64 const > const & dComponentDensDPres,
                   arrayView4d< real64 const > const & dComponentDensDComponentConc,
                   arrayView1d< real64 const > const & gravDepth,
                   arrayView2d< real64 const > const & dens,
                   arrayView2d< real64 const > const & dDensDPres,
                   arrayView2d< real64 const > const & dDensDProppantConc,
                   arrayView3d< real64 const > const & dDensDComponentConc,
                   arrayView2d< real64 const > const & visc,
                   arrayView2d< real64 const > const & dViscDPres,
                   arrayView2d< real64 const > const & dViscDProppantConc,
                   arrayView3d< real64 const > const & dViscDComponentConc,
                   arrayView2d< real64 const > const & fluidDensity,
                   arrayView2d< real64 const > const & dFluidDensDPres,
                   arrayView3d< real64 const > const & dFluidDensDComponentConc,
                   arrayView1d< real64 const > const & settlingFactor,
                   arrayView1d< real64 const > const & dSettlingFactorDPres,
                   arrayView1d< real64 const > const & dSettlingFactorDProppantConc,
                   arrayView2d< real64 const > const & dSettlingFactorDComponentConc,
                   arrayView1d< real64 const > const & collisionFactor,
                   arrayView1d< real64 const > const & dCollisionFactorDProppantConc,
                   arrayView1d< integer const > const & isProppantMobile,
                   arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( proppantPackVf ),
                   arrayView1d< real64 const > const & aperture,
                   R1Tensor const & unitGravityVector,
                   arrayView1d< R1Tensor const > const & transTMultiplier,
                   real64 const dt,
                   arraySlice1d< real64 > const & localFlux,
                   arraySlice2d< real64 > const & localFluxJacobian );

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
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
                        arrayView1d< R1Tensor > const & cellBasedFlux );
};

struct ProppantPackVolumeKernel
{

  template< typename VIEWTYPE >
  using ElementView = ElementRegionManager::ElementView< VIEWTYPE >;

  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

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
                                       ElementView< arrayView1d< real64 const > > const & settlingFactor,
                                       ElementView< arrayView2d< real64 const > > const & density,
                                       ElementView< arrayView2d< real64 const > > const & fluidDensity,
                                       ElementView< arrayView2d< real64 const > > const & fluidViscosity,
                                       ElementView< arrayView1d< integer const > > const & isProppantMobile,
                                       ElementView< arrayView1d< integer const > > const & isProppantBoundaryElement,
                                       ElementView< arrayView1d< real64 const > > const & aperture,
                                       ElementView< arrayView1d< real64 const > > const & volume,
                                       ElementView< arrayView1d< integer const > > const & elemGhostRank,
                                       ElementView< arrayView1d< R1Tensor const > > const & cellBasedFlux,
                                       ElementView< arrayView1d< real64 > > const & conc,
                                       ElementView< arrayView1d< real64 > > const & proppantPackVf,
                                       ElementView< arrayView1d< real64 > > const & proppantExcessPackV,
                                       ElementView< arrayView1d< real64 > > const & proppantLiftFlux );

  template< typename STENCIL_TYPE >
  static void
  LaunchProppantPackVolumeUpdate( STENCIL_TYPE const & stencil,
                                  R1Tensor const unitGravityVector,
                                  real64 const maxProppantConcentration,
                                  ElementView< arrayView1d< integer const > > const & isProppantMobile,
                                  ElementView< arrayView1d< real64 const > > const & proppantExcessPackV,
                                  ElementView< arrayView1d< real64 > > const & conc,
                                  ElementView< arrayView1d< real64 > > const & proppantPackVf );

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
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
                             arrayView2d< real64 const > const &,
                             arrayView1d< real64 const > const & volume,
                             arrayView1d< real64 const > const & aperture,
                             arrayView1d< integer const > const & elemGhostRank,
                             arrayView1d< integer const > const & isProppantBoundaryElement,
                             arrayView1d< integer const > const & isProppantMobile,
                             arrayView1d< R1Tensor const > const & cellBasedFlux,
                             arrayView1d< real64 > const & conc,
                             arrayView1d< real64 > const & proppantPackVf,
                             arrayView1d< real64 > const & proppantExcessPackV,
                             arrayView1d< real64 > const & proppantLiftFlux );

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  UpdateProppantPackVolume( localIndex const numElems,
                            arraySlice1d< localIndex const > const & stencilElementIndices,
                            arraySlice1d< real64 const > const & stencilWeights,
                            arraySlice1d< R1Tensor const > const & stencilCellCenterToEdgeCenters,
                            R1Tensor const unitGravityVector,
                            real64 const maxProppantConcentration,
                            arrayView1d< integer const > const & isProppantMobile,
                            arrayView1d< real64 const > const & proppantExcessPackV,
                            arrayView1d< real64 > const & conc,
                            arrayView1d< real64 > const & proppantPackVf );
};

} // namespace ProppantTransportKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_PROPPANTTRANSPORTKERNELS_HPP_
