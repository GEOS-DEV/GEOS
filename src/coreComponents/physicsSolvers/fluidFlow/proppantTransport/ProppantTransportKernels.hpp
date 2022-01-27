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
 * @file ProppantTransportKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_PROPPANTTRANSPORT_PROPPANTTRANSPORTKERNELS_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_PROPPANTTRANSPORT_PROPPANTTRANSPORTKERNELS_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "constitutive/fluid/ParticleFluidExtrinsicData.hpp"
#include "constitutive/fluid/SlurryFluidBase.hpp"
#include "constitutive/fluid/SlurryFluidExtrinsicData.hpp"
#include "constitutive/fluid/CellBasedFluxSlurryFluidAccessors.hpp"
#include "constitutive/permeability/PermeabilityExtrinsicData.hpp"
#include "constitutive/permeability/PermeabilityAccessors.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/CellBasedFluxFlowAccessors.hpp"
#include "physicsSolvers/fluidFlow/proppantTransport/ProppantTransportExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"

namespace geosx
{

namespace ProppantTransportKernels
{

/******************************** FluidUpdateKernel ********************************/

struct FluidUpdateKernel
{
  template< typename FLUID_WRAPPER >
  static void launch( FLUID_WRAPPER const & fluidWrapper,
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
        fluidWrapper.updateFluidProperty( a, q,
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
  static void launch( FLUID_WRAPPER const & fluidWrapper,
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
        fluidWrapper.updateComponentDensity( a, q,
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
  static void launch( PROPPANT_WRAPPER const & proppantWrapper,
                      arrayView1d< real64 const > const & proppantConc,
                      arrayView1d< real64 const > const & dProppantConc,
                      arrayView2d< real64 const > const & fluidDens,
                      arrayView2d< real64 const > const & dFluidDens_dPres,
                      arrayView3d< real64 const > const & dFluidDens_dCompConc,
                      arrayView2d< real64 const > const & fluidVisc,
                      arrayView2d< real64 const > const & dFluidVisc_dPres,
                      arrayView3d< real64 const > const & dFluidVisc_dCompConc )
  {
    forAll< parallelDevicePolicy<> >( proppantWrapper.numElems(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      proppantWrapper.update( a,
                              proppantConc[a] + dProppantConc[a],
                              fluidDens[a][0],
                              dFluidDens_dPres[a][0],
                              dFluidDens_dCompConc[a][0],
                              fluidVisc[a][0],
                              dFluidVisc_dPres[a][0],
                              dFluidVisc_dCompConc[a][0] );
    } );
  }
};

/******************************** AccumulationKernel ********************************/

struct AccumulationKernel
{
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  compute( localIndex const NC,
           real64 const proppantConcOld,
           real64 const proppantConcNew,
           arraySlice1d< real64 const > const & componentDensOld,
           arraySlice1d< real64 const > const & componentDensNew,
           arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM( dCompDens_dPres ),
           arraySlice2d< real64 const > const & dCompDens_dCompConc,
           real64 const volume,
           real64 const packPoreVolume,
           real64 const proppantLiftVolume,
           arraySlice1d< real64 > const & localAccum,
           arraySlice2d< real64 > const & localAccumJacobian );

  static void
  launch( localIndex const size,
          localIndex const NC,
          localIndex const NDOF,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & elemGhostRank,
          arrayView1d< real64 const > const & proppantConc,
          arrayView1d< real64 const > const & dProppantConc,
          arrayView2d< real64 const > const & componentDensOld,
          arrayView3d< real64 const > const & componentDens,
          arrayView3d< real64 const > const & dCompDens_dPres,
          arrayView4d< real64 const > const & dCompDens_dCompConc,
          arrayView1d< real64 const > const & volume,
          arrayView1d< real64 const > const & proppantPackVolFrac,
          arrayView1d< real64 const > const & proppantLiftFlux,
          real64 const dt,
          real64 const maxProppantConcentration,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );
};

/******************************** FluxKernel ********************************/

struct FluxKernel
{
  using FlowAccessors =
    StencilAccessors< extrinsicMeshData::ghostRank,
                      extrinsicMeshData::elementAperture,
                      extrinsicMeshData::flow::pressure,
                      extrinsicMeshData::flow::deltaPressure,
                      extrinsicMeshData::flow::gravityCoefficient,
                      extrinsicMeshData::proppant::proppantConcentration,
                      extrinsicMeshData::proppant::deltaProppantConcentration,
                      extrinsicMeshData::proppant::isProppantMobile >;

  using CellBasedFluxFlowAccessors = CellBasedFluxFlowAccessorsImpl;

  using ParticleFluidAccessors =
    StencilAccessors< extrinsicMeshData::particlefluid::settlingFactor,
                      extrinsicMeshData::particlefluid::dSettlingFactor_dPressure,
                      extrinsicMeshData::particlefluid::dSettlingFactor_dProppantConcentration,
                      extrinsicMeshData::particlefluid::dSettlingFactor_dComponentConcentration,
                      extrinsicMeshData::particlefluid::collisionFactor,
                      extrinsicMeshData::particlefluid::dCollisionFactor_dProppantConcentration >;

  using SlurryFluidAccessors =
    StencilAccessors< extrinsicMeshData::slurryfluid::density,
                      extrinsicMeshData::slurryfluid::dDensity_dPressure,
                      extrinsicMeshData::slurryfluid::dDensity_dProppantConcentration,
                      extrinsicMeshData::slurryfluid::dDensity_dComponentConcentration,
                      extrinsicMeshData::slurryfluid::viscosity,
                      extrinsicMeshData::slurryfluid::dViscosity_dPressure,
                      extrinsicMeshData::slurryfluid::dViscosity_dProppantConcentration,
                      extrinsicMeshData::slurryfluid::dViscosity_dComponentConcentration,
                      extrinsicMeshData::slurryfluid::componentDensity,
                      extrinsicMeshData::slurryfluid::dComponentDensity_dPressure,
                      extrinsicMeshData::slurryfluid::dComponentDensity_dComponentConcentration,
                      extrinsicMeshData::slurryfluid::fluidDensity,
                      extrinsicMeshData::slurryfluid::dFluidDensity_dPressure,
                      extrinsicMeshData::slurryfluid::dFluidDensity_dComponentConcentration >;

  using CellBasedFluxSlurryFluidAccessors = CellBasedFluxSlurryFluidAccessorsImpl;

  using PermeabilityAccessors = PermeabilityAccessorsImpl;

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

  static void
  launch( SurfaceElementStencilWrapper const & stencilWrapper,
          localIndex const numDofPerCell,
          real64 const dt,
          globalIndex const rankOffset,
          integer const updateProppantPacking,
          R1Tensor const & unitGravityVector,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
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
          ElementViewConst< arrayView3d< real64 const > > const & permeability,
          ElementViewConst< arrayView3d< real64 const > > const & permeabilityMultiplier,
          ElementViewConst< arrayView1d< real64 const > > const & aperture,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );

  static void
  launchCellBasedFluxCalculation( SurfaceElementStencilWrapper const & stencilWrapper,
                                  R1Tensor const & unitGravityVector,
                                  ElementViewConst< arrayView1d< real64 const > > const & pres,
                                  ElementViewConst< arrayView1d< real64 const > > const & gravDepth,
                                  ElementViewConst< arrayView2d< real64 const > > const & dens,
                                  ElementViewConst< arrayView2d< real64 const > > const & visc,
                                  ElementViewConst< arrayView3d< real64 const > > const & permeability,
                                  ElementViewConst< arrayView3d< real64 const > > const & permeabilityMultiplier,
                                  ElementViewConst< arrayView1d< real64 const > > const & aperture,
                                  ElementView< arrayView2d< real64 > > const & cellBasedFlux );

  /**
   * @brief Compute flux and its derivatives for a given multi-element connector.
   *
   * This is a specialized version that flux in a single region, and uses
   * element pairing instead of a proper junction.
   */
  template< localIndex MAX_NUM_FLUX_ELEMS >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  computeJunction( localIndex const numElems,
                   localIndex const numDofPerCell,
                   arraySlice1d< localIndex const > const & stencilElementIndices,
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
                   real64 const (&transmissibility)[MAX_NUM_FLUX_ELEMS],
                   real64 const (&apertureWeight)[MAX_NUM_FLUX_ELEMS],
                   real64 const (&geometricWeight)[MAX_NUM_FLUX_ELEMS],
                   real64 const dt,
                   arraySlice1d< real64 > const & localFlux,
                   arraySlice2d< real64 > const & localFluxJacobian );


  template< localIndex MAX_NUM_FLUX_ELEMS >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  computeCellBasedFlux( localIndex const numElems,
                        arraySlice1d< localIndex const > const & stencilElementIndices,
                        arrayView1d< real64 const > const & pres,
                        arrayView1d< real64 const > const & gravDepth,
                        arrayView2d< real64 const > const & dens,
                        arrayView2d< real64 const > const & visc,
                        arraySlice1d< R1Tensor const > const & cellCenterToEdgeCenters,
                        real64 const (&transmissibility)[MAX_NUM_FLUX_ELEMS],
                        real64 const (&apertureWeight)[MAX_NUM_FLUX_ELEMS],
                        real64 const (&geometricWeight)[MAX_NUM_FLUX_ELEMS],
                        arrayView2d< real64 > const & cellBasedFlux );
};

struct ProppantPackVolumeKernel
{

  using FlowAccessors =
    StencilAccessors< extrinsicMeshData::ghostRank,
                      extrinsicMeshData::elementAperture,
                      extrinsicMeshData::elementVolume,
                      extrinsicMeshData::proppant::cellBasedFlux,
                      extrinsicMeshData::proppant::isProppantMobile,
                      extrinsicMeshData::proppant::isProppantBoundary >;

  using ParticleFluidAccessors =
    StencilAccessors< extrinsicMeshData::particlefluid::settlingFactor >;

  using SlurryFluidAccessors =
    StencilAccessors< extrinsicMeshData::slurryfluid::density,
                      extrinsicMeshData::slurryfluid::fluidDensity,
                      extrinsicMeshData::slurryfluid::fluidViscosity >;


  template< typename VIEWTYPE >
  using ElementView = ElementRegionManager::ElementView< VIEWTYPE >;

  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  static void
  launchProppantPackVolumeCalculation( SurfaceElementStencil const & stencil,
                                       real64 const dt,
                                       real64 const proppantDensity,
                                       real64 const proppantDiameter,
                                       real64 const maxProppantConcentration,
                                       R1Tensor const & unitGravityVector,
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
                                       ElementView< arrayView2d< real64 const > > const & cellBasedFlux,
                                       ElementView< arrayView1d< real64 > > const & conc,
                                       ElementView< arrayView1d< real64 > > const & proppantPackVolFrac,
                                       ElementView< arrayView1d< real64 > > const & proppantExcessPackVolume,
                                       ElementView< arrayView1d< real64 > > const & proppantLiftFlux );

  static void
  launchProppantPackVolumeUpdate( SurfaceElementStencil const & stencil,
                                  R1Tensor const & unitGravityVector,
                                  real64 const maxProppantConcentration,
                                  ElementView< arrayView1d< integer const > > const & isProppantMobile,
                                  ElementView< arrayView1d< real64 const > > const & proppantExcessPackVolume,
                                  ElementView< arrayView1d< real64 > > const & conc,
                                  ElementView< arrayView1d< real64 > > const & proppantPackVolFrac );

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  computeProppantPackVolume( localIndex const numElems,
                             real64 const dt,
                             real64 const proppantDensity,
                             real64 const proppantDiameter,
                             real64 const maxProppantConcentration,
                             R1Tensor const & unitGravityVector,
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
                             arrayView2d< real64 const > const & cellBasedFlux,
                             arrayView1d< real64 > const & conc,
                             arrayView1d< real64 > const & proppantPackVolFrac,
                             arrayView1d< real64 > const & proppantExcessPackVolume,
                             arrayView1d< real64 > const & proppantLiftFlux );

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  updateProppantPackVolume( localIndex const numElems,
                            arraySlice1d< localIndex const > const & stencilElementIndices,
                            arraySlice1d< real64 const > const & stencilWeights,
                            arraySlice1d< R1Tensor const > const & stencilCellCenterToEdgeCenters,
                            R1Tensor const & unitGravityVector,
                            real64 const maxProppantConcentration,
                            arrayView1d< integer const > const & isProppantMobile,
                            arrayView1d< real64 const > const & proppantExcessPackVolume,
                            arrayView1d< real64 > const & conc,
                            arrayView1d< real64 > const & proppantPackVolFrac );
};

} // namespace ProppantTransportKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_PROPPANTTRANSPORT_PROPPANTTRANSPORTKERNELS_HPP_
