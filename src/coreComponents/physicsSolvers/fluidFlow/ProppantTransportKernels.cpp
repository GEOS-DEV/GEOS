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
 * @file ProppantTransportKernels.cpp
 */

#include "ProppantTransportKernels.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"

namespace geosx
{

namespace ProppantTransportKernels
{

inline void addLocalContributionsToGlobalSystem( localIndex const numFluxElems,
                                                 localIndex const stencilSize,
                                                 globalIndex const * const eqnRowIndices,
                                                 globalIndex const * const dofColIndices,
                                                 real64 const * const localFluxJacobian,
                                                 real64 const * const localFlux,
                                                 ParallelMatrix * const jacobian,
                                                 ParallelVector * const residual )
{

  // Add to global residual/jacobian
  jacobian->add( eqnRowIndices,
                 dofColIndices,
                 localFluxJacobian,
                 numFluxElems,
                 stencilSize );

  residual->add( eqnRowIndices,
                 localFlux,
                 numFluxElems );

}


template<>
void FluxKernel::
  Launch< CellElementStencilTPFA >( CellElementStencilTPFA const & GEOSX_UNUSED_PARAM( stencil ),
                                    localIndex const GEOSX_UNUSED_PARAM( numDofPerCell ),
                                    real64 const GEOSX_UNUSED_PARAM( dt ),
                                    FluxKernel::ElementViewConst< arrayView1d< R1Tensor const > > const & GEOSX_UNUSED_PARAM( transTMultiplier ),
                                    integer const GEOSX_UNUSED_PARAM( updateProppantPacking ),
                                    R1Tensor const & GEOSX_UNUSED_PARAM( unitGravityVector ),
                                    FluxKernel::ElementViewConst< arrayView1d< globalIndex const > > const & GEOSX_UNUSED_PARAM( dofNumber ),
                                    FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( pres ),
                                    FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dPres ),
                                    FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( proppantConc ),
                                    FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dProppantConc ),
                                    FluxKernel::ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( componentDens ),
                                    FluxKernel::ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dComponentDens_dPres ),
                                    FluxKernel::ElementViewConst< arrayView4d< real64 const > > const & GEOSX_UNUSED_PARAM( dComponentDens_dComponentConc ),
                                    FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( gravDepth ),
                                    FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dens ),
                                    FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dDens_dPres ),
                                    FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dDens_dProppantConc ),
                                    FluxKernel::ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dDens_dComponentConc ),
                                    FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( visc ),
                                    FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dVisc_dPres ),
                                    FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dVisc_dProppantConc ),
                                    FluxKernel::ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dVisc_dComponentConc ),
                                    FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( fluidDensity ),
                                    FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dFluidDens_dPres ),
                                    FluxKernel::ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dFluidDens_dComponentConc ),
                                    FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( settlingFactor ),
                                    FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dSettlingFactor_dPres ),
                                    FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dSettlingFactor_dProppantConc ),
                                    FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dSettlingFactor_dComponentConc ),
                                    FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( collisionFactor ),
                                    FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dCollisionFactor_dProppantConc ),
                                    FluxKernel::ElementViewConst< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM( isProppantMobile ),
                                    FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( proppantPackVf ),
                                    FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( aperture ),
                                    FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( proppantLiftFlux ),
                                    ParallelMatrix * const GEOSX_UNUSED_PARAM( jacobian ),
                                    ParallelVector * const GEOSX_UNUSED_PARAM( residual ) )
{}

template<>
void FluxKernel::
  Launch< FaceElementStencil >( FaceElementStencil const & stencil,
                                localIndex const numDofPerCell,
                                real64 const dt,
                                FluxKernel::ElementViewConst< arrayView1d< R1Tensor const > > const & transTMultiplier,
                                integer const updateProppantPacking,
                                R1Tensor const & unitGravityVector,
                                FluxKernel::ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
                                FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & pres,
                                FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & dPres,
                                FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & proppantConc,
                                FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & dProppantConc,
                                FluxKernel::ElementViewConst< arrayView3d< real64 const > > const & componentDens,
                                FluxKernel::ElementViewConst< arrayView3d< real64 const > > const & dComponentDens_dPres,
                                FluxKernel::ElementViewConst< arrayView4d< real64 const > > const & dComponentDens_dComponentConc,
                                FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & gravDepth,
                                FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & dens,
                                FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                                FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & dDens_dProppantConc,
                                FluxKernel::ElementViewConst< arrayView3d< real64 const > > const & dDens_dComponentConc,
                                FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & visc,
                                FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & dVisc_dPres,
                                FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & dVisc_dProppantConc,
                                FluxKernel::ElementViewConst< arrayView3d< real64 const > > const & dVisc_dComponentConc,
                                FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & fluidDensity,
                                FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & dFluidDens_dPres,
                                FluxKernel::ElementViewConst< arrayView3d< real64 const > > const & dFluidDens_dComponentConc,
                                FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & settlingFactor,
                                FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & dSettlingFactor_dPres,
                                FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & dSettlingFactor_dProppantConc,
                                FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & dSettlingFactor_dComponentConc,
                                FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & collisionFactor,
                                FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & dCollisionFactor_dProppantConc,
                                FluxKernel::ElementViewConst< arrayView1d< integer const > > const & isProppantMobile,
                                FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & proppantPackVf,
                                FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & aperture,
                                FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & proppantLiftFlux,
                                ParallelMatrix * const jacobian,
                                ParallelVector * const residual )
{
  constexpr localIndex maxNumFluxElems = FaceElementStencil::NUM_POINT_IN_FLUX;
  constexpr localIndex maxStencilSize = FaceElementStencil::MAX_STENCIL_SIZE;

  typename FaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename FaceElementStencil::WeightContainerViewConstType const & weights = stencil.getWeights();

  ArrayOfArraysView< R1Tensor const > const & cellCenterToEdgeCenters = stencil.getCellCenterToEdgeCenters();

  ArrayOfArraysView< integer const > const & isGhostConnectors = stencil.getIsGhostConnectors();

  constexpr localIndex DOF1 = maxNumFluxElems * constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;
  constexpr localIndex DOF2 = maxStencilSize * constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;

  forAll< serialPolicy >( stencil.size(), [=] ( localIndex iconn )
  {

    localIndex const numFluxElems = stencil.stencilSize( iconn );

    if( isGhostConnectors[iconn][0] < 0 && !(updateProppantPacking == 0 && numFluxElems <= 1) )
    {

      localIndex const stencilSize  = numFluxElems;

      localIndex const DOF = numFluxElems * numDofPerCell;

      // working arrays
      stackArray1d< globalIndex, DOF1 > eqnRowIndices( DOF );
      stackArray1d< globalIndex, DOF2 > dofColIndices( DOF );

      stackArray1d< real64, DOF1 > localFlux( DOF );
      stackArray2d< real64, DOF1 *DOF2 > localFluxJacobian( DOF, DOF );

      localIndex const er = seri[iconn][0];
      localIndex const esr = sesri[iconn][0];

      FluxKernel::ComputeJunction( numFluxElems,
                                   numDofPerCell,
                                   sei[iconn],
                                   weights[iconn],
                                   cellCenterToEdgeCenters[iconn],
                                   pres[er][esr],
                                   dPres[er][esr],
                                   proppantConc[er][esr],
                                   dProppantConc[er][esr],
                                   componentDens[er][esr],
                                   dComponentDens_dPres[er][esr],
                                   dComponentDens_dComponentConc[er][esr],
                                   gravDepth[er][esr],
                                   dens[er][esr],
                                   dDens_dPres[er][esr],
                                   dDens_dProppantConc[er][esr],
                                   dDens_dComponentConc[er][esr],
                                   visc[er][esr],
                                   dVisc_dPres[er][esr],
                                   dVisc_dProppantConc[er][esr],
                                   dVisc_dComponentConc[er][esr],
                                   fluidDensity[er][esr],
                                   dFluidDens_dPres[er][esr],
                                   dFluidDens_dComponentConc[er][esr],
                                   settlingFactor[er][esr],
                                   dSettlingFactor_dPres[er][esr],
                                   dSettlingFactor_dProppantConc[er][esr],
                                   dSettlingFactor_dComponentConc[er][esr],
                                   collisionFactor[er][esr],
                                   dCollisionFactor_dProppantConc[er][esr],
                                   isProppantMobile[er][esr],
                                   proppantPackVf[er][esr],
                                   aperture[er][esr],
                                   proppantLiftFlux[er][esr],
                                   unitGravityVector,
                                   transTMultiplier[er][esr],
                                   dt,
                                   localFlux,
                                   localFluxJacobian );


      for( localIndex i = 0; i < numFluxElems; ++i )
      {

        for( localIndex j = 0; j < numDofPerCell; ++j )
        {

          eqnRowIndices[i * numDofPerCell + j] = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] + j;

        }

      }

      for( localIndex i = 0; i < stencilSize; ++i )
      {

        for( localIndex j = 0; j < numDofPerCell; ++j )
        {

          dofColIndices[i * numDofPerCell + j] = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] + j;

        }


      }

      addLocalContributionsToGlobalSystem( numFluxElems * numDofPerCell,
                                           stencilSize * numDofPerCell,
                                           eqnRowIndices.data(),
                                           dofColIndices.data(),
                                           localFluxJacobian.data(),
                                           localFlux.data(),
                                           jacobian,
                                           residual );

    }
  } );
}

template<>
void FluxKernel::
  LaunchCellBasedFluxCalculation< CellElementStencilTPFA >( CellElementStencilTPFA const & GEOSX_UNUSED_PARAM( stencil ),
                                                            FluxKernel::ElementViewConst< arrayView1d< R1Tensor const > > const & GEOSX_UNUSED_PARAM(
                                                              transTMultiplier ),
                                                            R1Tensor const GEOSX_UNUSED_PARAM( unitGravityVector ),
                                                            FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( pres ),

                                                            FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( gravDepth ),
                                                            FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dens ),
                                                            FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( visc ),
                                                            FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( aperture ),
                                                            FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM(
                                                              proppantPackVf ),
                                                            FluxKernel::ElementView< arrayView1d< R1Tensor > > const & GEOSX_UNUSED_PARAM( cellBasedFlux ) )
{}

template<>
void FluxKernel::
  LaunchCellBasedFluxCalculation< FaceElementStencil >( FaceElementStencil const & stencil,
                                                        FluxKernel::ElementViewConst< arrayView1d< R1Tensor const > > const & transTMultiplier,
                                                        R1Tensor const unitGravityVector,
                                                        FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & pres,
                                                        FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & gravDepth,
                                                        FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & dens,
                                                        FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & visc,
                                                        FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & aperture,
                                                        FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & proppantPackVf,
                                                        FluxKernel::ElementView< arrayView1d< R1Tensor > > const & cellBasedFlux )
{

  typename FaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename FaceElementStencil::WeightContainerViewConstType const & weights = stencil.getWeights();

  ArrayOfArraysView< R1Tensor const > const & cellCenterToEdgeCenters = stencil.getCellCenterToEdgeCenters();


  forAll< serialPolicy >( stencil.size(), [=] ( localIndex iconn )
  {

    localIndex const numFluxElems = stencil.stencilSize( iconn );

    localIndex const er = seri[iconn][0];
    localIndex const esr = sesri[iconn][0];

    FluxKernel::ComputeCellBasedFlux( numFluxElems,
                                      sei[iconn],
                                      weights[iconn],
                                      cellCenterToEdgeCenters[iconn],
                                      transTMultiplier[er][esr],
                                      unitGravityVector,
                                      pres[er][esr],
                                      gravDepth[er][esr],
                                      dens[er][esr],
                                      visc[er][esr],
                                      aperture[er][esr],
                                      proppantPackVf[er][esr],
                                      cellBasedFlux[er][esr] );

  } );
}


template<>
void ProppantPackVolumeKernel::
  LaunchProppantPackVolumeCalculation< CellElementStencilTPFA >( CellElementStencilTPFA const & GEOSX_UNUSED_PARAM( stencil ),
                                                                 real64 const GEOSX_UNUSED_PARAM( dt ),
                                                                 real64 const GEOSX_UNUSED_PARAM( proppantDensity ),
                                                                 real64 const GEOSX_UNUSED_PARAM( proppantDiameter ),
                                                                 real64 const GEOSX_UNUSED_PARAM( maxProppantConcentration ),
                                                                 R1Tensor const GEOSX_UNUSED_PARAM( unitGravityVector ),
                                                                 real64 const GEOSX_UNUSED_PARAM( criticalShieldsNumber ),
                                                                 real64 const GEOSX_UNUSED_PARAM( frictionCoefficient ),
                                                                 ProppantPackVolumeKernel::ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM(
                                                                   conc ),
                                                                 ProppantPackVolumeKernel::ElementView< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM(
                                                                   settlingFactor ),
                                                                 ProppantPackVolumeKernel::ElementView< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM(
                                                                   density ),
                                                                 ProppantPackVolumeKernel::ElementView< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM(
                                                                   fluidDensity ),
                                                                 ProppantPackVolumeKernel::ElementView< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM(
                                                                   fluidViscosity ),
                                                                 ProppantPackVolumeKernel::ElementView< arrayView1d< integer > > const & GEOSX_UNUSED_PARAM(
                                                                   isProppantMobile ),
                                                                 ProppantPackVolumeKernel::ElementViewConst< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM(
                                                                   isProppantBoundaryElement ),
                                                                 ProppantPackVolumeKernel::ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM(
                                                                   proppantPackVf ),
                                                                 ProppantPackVolumeKernel::ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM(
                                                                   proppantExcessPackV ),
                                                                 ProppantPackVolumeKernel::ElementView< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM(
                                                                   aperture ),
                                                                 ProppantPackVolumeKernel::ElementView< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM(
                                                                   volume ),
                                                                 ProppantPackVolumeKernel::ElementView< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM(
                                                                   elemGhostRank ),
                                                                 ProppantPackVolumeKernel::ElementView< arrayView1d< R1Tensor const > > const & GEOSX_UNUSED_PARAM(
                                                                   cellBasedFlux ),
                                                                 ProppantPackVolumeKernel::ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM(
                                                                   proppantLiftFlux ) )
{}


template<>
void ProppantPackVolumeKernel::
  LaunchProppantPackVolumeCalculation< FaceElementStencil >( FaceElementStencil const & stencil,
                                                             real64 const dt,
                                                             real64 const proppantDensity,
                                                             real64 const proppantDiameter,
                                                             real64 const maxProppantConcentration,
                                                             R1Tensor const unitGravityVector,
                                                             real64 const criticalShieldsNumber,
                                                             real64 const frictionCoefficient,
                                                             ProppantPackVolumeKernel::ElementView< arrayView1d< real64 > > const & conc,
                                                             ProppantPackVolumeKernel::ElementView< arrayView1d< real64 const > > const & settlingFactor,
                                                             ProppantPackVolumeKernel::ElementView< arrayView2d< real64 const > > const & density,
                                                             ProppantPackVolumeKernel::ElementView< arrayView2d< real64 const > > const & fluidDensity,
                                                             ProppantPackVolumeKernel::ElementView< arrayView2d< real64 const > > const & fluidViscosity,
                                                             ProppantPackVolumeKernel::ElementView< arrayView1d< integer > > const & isProppantMobile,
                                                             ProppantPackVolumeKernel::ElementView< arrayView1d< integer const > > const & isProppantBoundaryElement,
                                                             ProppantPackVolumeKernel::ElementView< arrayView1d< real64 > > const & proppantPackVf,
                                                             ProppantPackVolumeKernel::ElementView< arrayView1d< real64 > > const & proppantExcessPackV,
                                                             ProppantPackVolumeKernel::ElementView< arrayView1d< real64 const > > const & aperture,
                                                             ProppantPackVolumeKernel::ElementView< arrayView1d< real64 const > > const & volume,
                                                             ProppantPackVolumeKernel::ElementView< arrayView1d< integer const > > const & elemGhostRank,
                                                             ProppantPackVolumeKernel::ElementView< arrayView1d< R1Tensor const > > const & cellBasedFlux,
                                                             ProppantPackVolumeKernel::ElementView< arrayView1d< real64 > > const & proppantLiftFlux )
{

  typename FaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename FaceElementStencil::WeightContainerViewConstType const & weights = stencil.getWeights();

  ArrayOfArraysView< R1Tensor const > const & cellCenterToEdgeCenters = stencil.getCellCenterToEdgeCenters();

  forAll< serialPolicy >( stencil.size(), [=] ( localIndex iconn )
  {

    localIndex const numFluxElems = stencil.stencilSize( iconn );

    localIndex const er = seri[iconn][0];
    localIndex const esr = sesri[iconn][0];

    ProppantPackVolumeKernel::ComputeProppantPackVolume( numFluxElems,
                                                         dt,
                                                         proppantDensity,
                                                         proppantDiameter,
                                                         maxProppantConcentration,
                                                         unitGravityVector,
                                                         criticalShieldsNumber,
                                                         frictionCoefficient,
                                                         sei[iconn],
                                                         weights[iconn],
                                                         cellCenterToEdgeCenters[iconn],
                                                         settlingFactor[er][esr],
                                                         density[er][esr],
                                                         fluidDensity[er][esr],
                                                         fluidViscosity[er][esr],
                                                         volume[er][esr],
                                                         aperture[er][esr],
                                                         elemGhostRank[er][esr],
                                                         isProppantBoundaryElement[er][esr],
                                                         conc[er][esr],
                                                         isProppantMobile[er][esr],
                                                         proppantPackVf[er][esr],
                                                         proppantExcessPackV[er][esr],
                                                         cellBasedFlux[er][esr],
                                                         proppantLiftFlux[er][esr] );

  } );

}


template<>
void ProppantPackVolumeKernel::
  LaunchProppantPackVolumeUpdate< CellElementStencilTPFA >( CellElementStencilTPFA const & GEOSX_UNUSED_PARAM( stencil ),
                                                            R1Tensor const GEOSX_UNUSED_PARAM( unitGravityVector ),
                                                            real64 const GEOSX_UNUSED_PARAM( maxProppantConcentration ),
                                                            ProppantPackVolumeKernel::ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM( conc ),
                                                            ProppantPackVolumeKernel::ElementView< arrayView1d< integer > > const & GEOSX_UNUSED_PARAM(
                                                              isProppantMobile ),
                                                            ProppantPackVolumeKernel::ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM(
                                                              proppantPackVf ),
                                                            ProppantPackVolumeKernel::ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM(
                                                              proppantExcessPackV ) )
{}

template<>
void ProppantPackVolumeKernel::
  LaunchProppantPackVolumeUpdate< FaceElementStencil >( FaceElementStencil const & stencil,
                                                        R1Tensor const unitGravityVector,
                                                        real64 const maxProppantConcentration,
                                                        ProppantPackVolumeKernel::ElementView< arrayView1d< real64 > > const & conc,
                                                        ProppantPackVolumeKernel::ElementView< arrayView1d< integer > > const & isProppantMobile,
                                                        ProppantPackVolumeKernel::ElementView< arrayView1d< real64 > > const & proppantPackVf,
                                                        ProppantPackVolumeKernel::ElementView< arrayView1d< real64 > > const & proppantExcessPackV )
{

  typename FaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename FaceElementStencil::WeightContainerViewConstType const & weights = stencil.getWeights();

  ArrayOfArraysView< R1Tensor const > const & cellCenterToEdgeCenters = stencil.getCellCenterToEdgeCenters();

  forAll< serialPolicy >( stencil.size(), [=] ( localIndex iconn )
  {

    localIndex const numFluxElems = stencil.stencilSize( iconn );

    localIndex const er = seri[iconn][0];
    localIndex const esr = sesri[iconn][0];

    ProppantPackVolumeKernel::UpdateProppantPackVolume( numFluxElems,
                                                        sei[iconn],
                                                        weights[iconn],
                                                        cellCenterToEdgeCenters[iconn],
                                                        unitGravityVector,
                                                        maxProppantConcentration,
                                                        conc[er][esr],
                                                        isProppantMobile[er][esr],
                                                        proppantPackVf[er][esr],
                                                        proppantExcessPackV[er][esr] );

  } );

}

} // namespace ProppantTransportKernels

} // namespace geosx
