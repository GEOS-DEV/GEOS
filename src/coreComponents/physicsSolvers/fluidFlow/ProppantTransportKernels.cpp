/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
                 numFluxElems);

}


template<>
void FluxKernel::
Launch<CellElementStencilTPFA>( CellElementStencilTPFA const & GEOSX_UNUSED_ARG(stencil),
				localIndex const GEOSX_UNUSED_ARG(numDofPerCell),
				real64 const GEOSX_UNUSED_ARG(dt),
				localIndex const GEOSX_UNUSED_ARG(fluidIndex),
				localIndex const GEOSX_UNUSED_ARG(proppantIndex),
				real64 const GEOSX_UNUSED_ARG(proppantPackPermeability),
                                R1Tensor const & GEOSX_UNUSED_ARG(unitGravityVector),                                
				FluxKernel::ElementViewConst< arrayView1d<globalIndex const > > const & GEOSX_UNUSED_ARG( dofNumber ),
				FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & GEOSX_UNUSED_ARG( pres ),
				FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & GEOSX_UNUSED_ARG( dPres ),
				FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & GEOSX_UNUSED_ARG( proppantConc ),
				FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & GEOSX_UNUSED_ARG( dProppantConc ),
				FluxKernel::MaterialView< arrayView3d<real64 const> > const & GEOSX_UNUSED_ARG( componentDens ),
				FluxKernel::MaterialView< arrayView3d<real64 const> > const & GEOSX_UNUSED_ARG( dComponentDens_dPres ),
				FluxKernel::MaterialView< arrayView4d<real64 const> > const & GEOSX_UNUSED_ARG( dComponentDens_dComponentConc ),
                                FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & GEOSX_UNUSED_ARG( gravDepth ),
				FluxKernel::MaterialView< arrayView2d<real64 const> > const & GEOSX_UNUSED_ARG( dens ),
				FluxKernel::MaterialView< arrayView2d<real64 const> > const & GEOSX_UNUSED_ARG( dDens_dPres ),
				FluxKernel::MaterialView< arrayView2d<real64 const> > const & GEOSX_UNUSED_ARG( dDens_dProppantConc ),
				FluxKernel::MaterialView< arrayView3d<real64 const> > const & GEOSX_UNUSED_ARG( dDens_dComponentConc ),
				FluxKernel::MaterialView< arrayView2d<real64 const> > const & GEOSX_UNUSED_ARG( visc ),
				FluxKernel::MaterialView< arrayView2d<real64 const> > const & GEOSX_UNUSED_ARG( dVisc_dPres ),
				FluxKernel::MaterialView< arrayView2d<real64 const> > const & GEOSX_UNUSED_ARG( dVisc_dProppantConc ),
				FluxKernel::MaterialView< arrayView3d<real64 const> > const & GEOSX_UNUSED_ARG( dVisc_dComponentConc ),
				FluxKernel::MaterialView< arrayView2d<real64 const> > const & GEOSX_UNUSED_ARG( fluidDensity ),
				FluxKernel::MaterialView< arrayView2d<real64 const> > const & GEOSX_UNUSED_ARG( dFluidDens_dPres ),
				FluxKernel::MaterialView< arrayView3d<real64 const> > const & GEOSX_UNUSED_ARG( dFluidDens_dComponentConc ),
				FluxKernel::MaterialView< arrayView1d<real64 const> > const & GEOSX_UNUSED_ARG( settlingFactor ),
				FluxKernel::MaterialView< arrayView1d<real64 const> > const & GEOSX_UNUSED_ARG( dSettlingFactor_dPres ),
				FluxKernel::MaterialView< arrayView1d<real64 const> > const & GEOSX_UNUSED_ARG( dSettlingFactor_dProppantConc ),
				FluxKernel::MaterialView< arrayView2d<real64 const> > const & GEOSX_UNUSED_ARG( dSettlingFactor_dComponentConc ),
				FluxKernel::MaterialView< arrayView1d<real64 const> > const & GEOSX_UNUSED_ARG( collisionFactor ),
				FluxKernel::MaterialView< arrayView1d<real64 const> > const & GEOSX_UNUSED_ARG( dCollisionFactor_dProppantConc ),
				FluxKernel::ElementViewConst< arrayView1d<integer const> > const & GEOSX_UNUSED_ARG( isProppantMobile ),
				FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & GEOSX_UNUSED_ARG( proppantPackVf ),
				FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & GEOSX_UNUSED_ARG( aperture ),
				FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & GEOSX_UNUSED_ARG( proppantLiftFlux ),
				FluxKernel::ElementViewConst< arrayView1d<integer const> > const & GEOSX_UNUSED_ARG( isInterfaceElement ),
				ParallelMatrix * const GEOSX_UNUSED_ARG(jacobian),
				ParallelVector * const GEOSX_UNUSED_ARG(residual) )
{

}

template<>
void FluxKernel::
Launch<FaceElementStencil>( FaceElementStencil const & stencil,
			    localIndex const numDofPerCell,
			    real64 const dt,
			    localIndex const fluidIndex,
			    localIndex const proppantIndex,
                            real64 const proppantPackPermeability,
                            R1Tensor const & unitGravityVector,
			    FluxKernel::ElementViewConst< arrayView1d<globalIndex const > > const & dofNumber,
			    FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & pres,
			    FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & dPres,
			    FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & proppantConc,
			    FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & dProppantConc,
			    FluxKernel::MaterialView< arrayView3d<real64 const> > const & componentDens,
			    FluxKernel::MaterialView< arrayView3d<real64 const> > const & dComponentDens_dPres,
			    FluxKernel::MaterialView< arrayView4d<real64 const> > const & dComponentDens_dComponentConc,                        
			    FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & gravDepth,
			    FluxKernel::MaterialView< arrayView2d<real64 const> > const & dens,
			    FluxKernel::MaterialView< arrayView2d<real64 const> > const & dDens_dPres,
			    FluxKernel::MaterialView< arrayView2d<real64 const> > const & dDens_dProppantConc,
			    FluxKernel::MaterialView< arrayView3d<real64 const> > const & dDens_dComponentConc,                            
			    FluxKernel::MaterialView< arrayView2d<real64 const> > const & visc,
			    FluxKernel::MaterialView< arrayView2d<real64 const> > const & dVisc_dPres,
			    FluxKernel::MaterialView< arrayView2d<real64 const> > const & dVisc_dProppantConc,
			    FluxKernel::MaterialView< arrayView3d<real64 const> > const & dVisc_dComponentConc,                            
			    FluxKernel::MaterialView< arrayView2d<real64 const> > const & fluidDensity,
			    FluxKernel::MaterialView< arrayView2d<real64 const> > const & dFluidDens_dPres,
			    FluxKernel::MaterialView< arrayView3d<real64 const> > const & dFluidDens_dComponentConc,	                              
			    FluxKernel::MaterialView< arrayView1d<real64 const> > const & settlingFactor,
			    FluxKernel::MaterialView< arrayView1d<real64 const> > const & dSettlingFactor_dPres,
			    FluxKernel::MaterialView< arrayView1d<real64 const> > const & dSettlingFactor_dProppantConc,
			    FluxKernel::MaterialView< arrayView2d<real64 const> > const & dSettlingFactor_dComponentConc,                            
			    FluxKernel::MaterialView< arrayView1d<real64 const> > const & collisionFactor,
			    FluxKernel::MaterialView< arrayView1d<real64 const> > const & dCollisionFactor_dProppantConc,
			    FluxKernel::ElementViewConst< arrayView1d<integer const> > const & isProppantMobile,	  
			    FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & proppantPackVf,
			    FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & aperture,
			    FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & proppantLiftFlux,
			    FluxKernel::ElementViewConst< arrayView1d<integer const> > const & isInterfaceElement,
                            ParallelMatrix * const jacobian,
                            ParallelVector * const residual )
{
  constexpr localIndex maxNumFluxElems = FaceElementStencil::NUM_POINT_IN_FLUX;
  constexpr localIndex maxStencilSize = FaceElementStencil::MAX_STENCIL_SIZE;

  typename FaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename FaceElementStencil::WeightContainerViewConstType const & weights = stencil.getWeights();

  ArrayOfArraysView<R1Tensor const> const & cellCenterToEdgeCenters = stencil.getCellCenterToEdgeCenters();

  ArrayOfArraysView<integer const> const & isGhostConnectors = stencil.getIsGhostConnectors();  
  
  constexpr  localIndex DOF1 = maxNumFluxElems * constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;
  constexpr  localIndex DOF2 = maxStencilSize * constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;

  forall_in_range<serialPolicy>( 0, stencil.size(), GEOSX_LAMBDA ( localIndex iconn )
  {

    if(isGhostConnectors[iconn][0] < 0) 
      {

        localIndex const numFluxElems = stencil.stencilSize(iconn);
        localIndex const stencilSize  = numFluxElems;

        localIndex const DOF = numFluxElems * numDofPerCell;
    
        // working arrays
        stackArray1d<globalIndex, DOF1> eqnRowIndices(DOF);
        stackArray1d<globalIndex, DOF2> dofColIndices(DOF);

        stackArray1d<real64, DOF1> localFlux(DOF);
        stackArray2d<real64, DOF1*DOF2> localFluxJacobian(DOF, DOF);

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
                                     componentDens[er][esr][fluidIndex],
                                     dComponentDens_dPres[er][esr][fluidIndex],
                                     dComponentDens_dComponentConc[er][esr][fluidIndex],
                                     gravDepth[er][esr],
                                     dens[er][esr][fluidIndex],
                                     dDens_dPres[er][esr][fluidIndex],
                                     dDens_dProppantConc[er][esr][fluidIndex],
                                     dDens_dComponentConc[er][esr][fluidIndex],
                                     visc[er][esr][fluidIndex],
                                     dVisc_dPres[er][esr][fluidIndex],
                                     dVisc_dProppantConc[er][esr][fluidIndex],
                                     dVisc_dComponentConc[er][esr][fluidIndex],
                                     fluidDensity[er][esr][fluidIndex],
                                     dFluidDens_dPres[er][esr][fluidIndex],
                                     dFluidDens_dComponentConc[er][esr][fluidIndex],
                                     settlingFactor[er][esr][proppantIndex],
                                     dSettlingFactor_dPres[er][esr][proppantIndex],
                                     dSettlingFactor_dProppantConc[er][esr][proppantIndex],
                                     dSettlingFactor_dComponentConc[er][esr][proppantIndex],
                                     collisionFactor[er][esr][proppantIndex],
                                     dCollisionFactor_dProppantConc[er][esr][proppantIndex],
                                     isProppantMobile[er][esr],
                                     proppantPackVf[er][esr],
                                     aperture[er][esr],
                                     proppantLiftFlux[er][esr],
                                     isInterfaceElement[er][esr],
                                     unitGravityVector,
                                     proppantPackPermeability,
                                     dt,
                                     localFlux,
                                     localFluxJacobian);


          for (localIndex i = 0; i < numFluxElems; ++i)
            {

              for (localIndex j = 0; j < numDofPerCell; ++j)
                {

                  eqnRowIndices[i * numDofPerCell + j] = dofNumber[seri(iconn,i)][sesri(iconn,i)][sei(iconn,i)] + j;

                }
	
            }

        for (localIndex i = 0; i < stencilSize; ++i)
          {

            for (localIndex j = 0; j < numDofPerCell; ++j)
              {

                dofColIndices[i * numDofPerCell + j] = dofNumber[seri(iconn,i)][sesri(iconn,i)][sei(iconn,i)] + j;

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
LaunchCellBasedFluxCalculation<CellElementStencilTPFA>( CellElementStencilTPFA const & GEOSX_UNUSED_ARG(stencil),
                                                    localIndex const GEOSX_UNUSED_ARG(fluidIndex),
                                                    real64 const GEOSX_UNUSED_ARG(proppantPackPermeability),
                                                    R1Tensor const GEOSX_UNUSED_ARG(unitGravityVector),
                                                    FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & GEOSX_UNUSED_ARG(pres),

                                                    FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & GEOSX_UNUSED_ARG(gravDepth),
                                                    FluxKernel::MaterialView< arrayView2d<real64 const> > const & GEOSX_UNUSED_ARG(dens),
                                                    FluxKernel::MaterialView< arrayView2d<real64 const> > const & GEOSX_UNUSED_ARG(visc),
                                                    FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & GEOSX_UNUSED_ARG(aperture),
                                                    FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & GEOSX_UNUSED_ARG(proppantPackVf),
                                                    FluxKernel::ElementView< arrayView1d<R1Tensor> > const & GEOSX_UNUSED_ARG(cellBasedFlux) )
{

}
  
template<>
void FluxKernel::
LaunchCellBasedFluxCalculation<FaceElementStencil>( FaceElementStencil const & stencil,
                                                    localIndex const fluidIndex,
                                                    real64 const proppantPackPermeability,
                                                    R1Tensor const unitGravityVector,                                                    
                                                    FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & pres,

                                                    FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & gravDepth,
                                                    FluxKernel::MaterialView< arrayView2d<real64 const> > const & dens,
                                                    FluxKernel::MaterialView< arrayView2d<real64 const> > const & visc,
                                                    FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & aperture,
                                                    FluxKernel::ElementViewConst< arrayView1d<real64 const> > const & proppantPackVf,
                                                    FluxKernel::ElementView< arrayView1d<R1Tensor> > const & cellBasedFlux)
{

  typename FaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename FaceElementStencil::WeightContainerViewConstType const & weights = stencil.getWeights();

  ArrayOfArraysView<R1Tensor const> const & cellCenterToEdgeCenters = stencil.getCellCenterToEdgeCenters();

  
  forall_in_range<serialPolicy>( 0, stencil.size(), GEOSX_LAMBDA ( localIndex iconn )
  {

    localIndex const numFluxElems = stencil.stencilSize(iconn);

    localIndex const er = seri[iconn][0];
    localIndex const esr = sesri[iconn][0];

    FluxKernel::ComputeCellBasedFlux( numFluxElems,
                                      sei[iconn],
                                      weights[iconn],
                                      cellCenterToEdgeCenters[iconn],
                                      proppantPackPermeability,
                                      unitGravityVector,
                                      pres[er][esr],
                                      gravDepth[er][esr],
                                      dens[er][esr][fluidIndex],
                                      visc[er][esr][fluidIndex],
                                      aperture[er][esr],
                                      proppantPackVf[er][esr],
                                      cellBasedFlux[er][esr]);

  } );
}


template<>
void ProppantPackVolumeKernel::
LaunchProppantPackVolumeCalculation<CellElementStencilTPFA>( CellElementStencilTPFA const & GEOSX_UNUSED_ARG(stencil),
                                                         real64 const GEOSX_UNUSED_ARG(dt),
                                                         localIndex const GEOSX_UNUSED_ARG(fluidIndex),
                                                         localIndex const GEOSX_UNUSED_ARG(proppantIndex),
                                                         real64 const GEOSX_UNUSED_ARG(proppantDensity),
                                                         real64 const GEOSX_UNUSED_ARG(proppantDiameter),
                                                         real64 const GEOSX_UNUSED_ARG(maxProppantConcentration),
                                                         R1Tensor const GEOSX_UNUSED_ARG(unitGravityVector),
                                                         real64 const GEOSX_UNUSED_ARG(criticalShieldsNumber),
                                                         ProppantPackVolumeKernel::ElementView < arrayView1d<real64> > const & GEOSX_UNUSED_ARG(conc),
                                                         ProppantPackVolumeKernel::MaterialViewConst< arrayView1d<real64 const> > const & GEOSX_UNUSED_ARG(settlingFactor),
                                                         ProppantPackVolumeKernel::MaterialViewConst< arrayView2d<real64 const> > const & GEOSX_UNUSED_ARG(fluidDensity),
                                                         ProppantPackVolumeKernel::MaterialViewConst< arrayView2d<real64 const> > const & GEOSX_UNUSED_ARG(fluidViscosity),
                                                         ProppantPackVolumeKernel::ElementView < arrayView1d<integer> > const & GEOSX_UNUSED_ARG(isProppantMobile),
                                                         ProppantPackVolumeKernel::ElementViewConst < arrayView1d<integer const> > const & GEOSX_UNUSED_ARG(isProppantBoundaryElement),
                                                         ProppantPackVolumeKernel::ElementView < arrayView1d<real64> > const & GEOSX_UNUSED_ARG(proppantPackVf),
                                                         ProppantPackVolumeKernel::ElementView < arrayView1d<real64> > const & GEOSX_UNUSED_ARG(proppantExcessPackV),
                                                         ProppantPackVolumeKernel::ElementView < arrayView1d<real64 const> > const & GEOSX_UNUSED_ARG(aperture),
                                                         ProppantPackVolumeKernel::ElementView < arrayView1d<real64 const> > const & GEOSX_UNUSED_ARG(volume),
                                                         ProppantPackVolumeKernel::ElementView < arrayView1d<integer const> > const & GEOSX_UNUSED_ARG(elemGhostRank),
                                                         ProppantPackVolumeKernel::ElementView < arrayView1d<R1Tensor> > const & GEOSX_UNUSED_ARG(cellBasedFlux),
                                                         ProppantPackVolumeKernel::ElementView < arrayView1d<real64> > const & GEOSX_UNUSED_ARG(proppantLiftFlux) )
{


}


template<>
void ProppantPackVolumeKernel::
LaunchProppantPackVolumeCalculation<FaceElementStencil>( FaceElementStencil const & stencil,
                                                         real64 const dt,
                                                         localIndex const fluidIndex,
                                                         localIndex const proppantIndex,
                                                         real64 const proppantDensity,
                                                         real64 const proppantDiameter,
                                                         real64 const maxProppantConcentration,
                                                         R1Tensor const unitGravityVector,
                                                         real64 const criticalShieldsNumber,
                                                         ProppantPackVolumeKernel::ElementView < arrayView1d<real64> > const & conc,
                                                         ProppantPackVolumeKernel::MaterialViewConst< arrayView1d<real64 const> > const & settlingFactor,
                                                         ProppantPackVolumeKernel::MaterialViewConst< arrayView2d<real64 const> > const & fluidDensity,
                                                         ProppantPackVolumeKernel::MaterialViewConst< arrayView2d<real64 const> > const & fluidViscosity,
                                                         ProppantPackVolumeKernel::ElementView < arrayView1d<integer> > const & isProppantMobile,
                                                         ProppantPackVolumeKernel::ElementViewConst < arrayView1d<integer const> > const & isProppantBoundaryElement,                                                         
                                                         ProppantPackVolumeKernel::ElementView < arrayView1d<real64> > const & proppantPackVf,
                                                         ProppantPackVolumeKernel::ElementView < arrayView1d<real64> > const & proppantExcessPackV,
                                                         ProppantPackVolumeKernel::ElementView < arrayView1d<real64 const> > const & aperture,
                                                         ProppantPackVolumeKernel::ElementView < arrayView1d<real64 const> > const & volume,
                                                         ProppantPackVolumeKernel::ElementView < arrayView1d<integer const> > const & elemGhostRank,          
                                                         ProppantPackVolumeKernel::ElementView < arrayView1d<R1Tensor> > const & cellBasedFlux,
                                                         ProppantPackVolumeKernel::ElementView < arrayView1d<real64> > const & proppantLiftFlux)
{

  typename FaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename FaceElementStencil::WeightContainerViewConstType const & weights = stencil.getWeights();

  ArrayOfArraysView<R1Tensor const> const & cellCenterToEdgeCenters = stencil.getCellCenterToEdgeCenters();

  forall_in_range<serialPolicy>( 0, stencil.size(), GEOSX_LAMBDA ( localIndex iconn )
  {

    localIndex const numFluxElems = stencil.stencilSize(iconn);

    localIndex const er = seri[iconn][0];
    localIndex const esr = sesri[iconn][0];

    ProppantPackVolumeKernel::ComputeProppantPackVolume( numFluxElems,
                                                         dt,
                                                         proppantDensity,
                                                         proppantDiameter,
                                                         maxProppantConcentration,
                                                         unitGravityVector,
                                                         criticalShieldsNumber,
                                                         sei[iconn],
                                                         weights[iconn],
                                                         cellCenterToEdgeCenters[iconn],
                                                         settlingFactor[er][esr][proppantIndex],
                                                         fluidDensity[er][esr][fluidIndex],
                                                         fluidViscosity[er][esr][fluidIndex],
                                                         volume[er][esr],
                                                         aperture[er][esr],
                                                         elemGhostRank[er][esr],
                                                         isProppantBoundaryElement[er][esr],
                                                         conc[er][esr],
                                                         isProppantMobile[er][esr],                                                         
                                                         proppantPackVf[er][esr],
                                                         proppantExcessPackV[er][esr],
                                                         cellBasedFlux[er][esr],
                                                         proppantLiftFlux[er][esr]);    

  } );
  
}


template<>
void ProppantPackVolumeKernel::
LaunchProppantPackVolumeUpdate<CellElementStencilTPFA>( CellElementStencilTPFA const & GEOSX_UNUSED_ARG(stencil),
                                                    R1Tensor const GEOSX_UNUSED_ARG(unitGravityVector),
                                                    real64 const GEOSX_UNUSED_ARG(maxProppantConcentration),                                 
                                                    ProppantPackVolumeKernel::ElementView < arrayView1d<real64> > const & GEOSX_UNUSED_ARG(conc),
                                                    ProppantPackVolumeKernel::ElementView < arrayView1d<integer> > const & GEOSX_UNUSED_ARG(isProppantMobile),
                                                    ProppantPackVolumeKernel::ElementView < arrayView1d<real64> > const & GEOSX_UNUSED_ARG(proppantPackVf),
                                                    ProppantPackVolumeKernel::ElementView < arrayView1d<real64> > const & GEOSX_UNUSED_ARG(proppantExcessPackV) )
{


}  

template<>
void ProppantPackVolumeKernel::
LaunchProppantPackVolumeUpdate<FaceElementStencil>( FaceElementStencil const & stencil,
                                                    R1Tensor const unitGravityVector,
                                                    real64 const maxProppantConcentration,
                                                    ProppantPackVolumeKernel::ElementView < arrayView1d<real64> > const & conc,
                                                    ProppantPackVolumeKernel::ElementView < arrayView1d<integer> > const & isProppantMobile,
                                                    ProppantPackVolumeKernel::ElementView < arrayView1d<real64> > const & proppantPackVf,
                                                    ProppantPackVolumeKernel::ElementView < arrayView1d<real64> > const & proppantExcessPackV)
{

  typename FaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename FaceElementStencil::WeightContainerViewConstType const & weights = stencil.getWeights();

  ArrayOfArraysView<R1Tensor const> const & cellCenterToEdgeCenters = stencil.getCellCenterToEdgeCenters();

  forall_in_range<serialPolicy>( 0, stencil.size(), GEOSX_LAMBDA ( localIndex iconn )
  {

    localIndex const numFluxElems = stencil.stencilSize(iconn);

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
                                                        proppantExcessPackV[er][esr]);

  } );
  
}

template<>
void ProppantPackVolumeKernel::
LaunchInterfaceElementUpdate<CellElementStencilTPFA>( CellElementStencilTPFA const & GEOSX_UNUSED_ARG(stencil),
                                                  R1Tensor const GEOSX_UNUSED_ARG(unitGravityVector),
                                                  ProppantPackVolumeKernel::ElementView < arrayView1d<integer> > const & GEOSX_UNUSED_ARG(isProppantMobile),
                                                  ProppantPackVolumeKernel::ElementView < arrayView1d<integer> > const & GEOSX_UNUSED_ARG(isInterfaceElement) )
{


}

template<>
void ProppantPackVolumeKernel::
LaunchInterfaceElementUpdate<FaceElementStencil>( FaceElementStencil const & stencil,
                                                  R1Tensor const unitGravityVector,
                                                  ProppantPackVolumeKernel::ElementView < arrayView1d<integer> > const & isProppantMobile,
                                                  ProppantPackVolumeKernel::ElementView < arrayView1d<integer> > const & isInterfaceElement )
{

  typename FaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename FaceElementStencil::WeightContainerViewConstType const & weights = stencil.getWeights();

  ArrayOfArraysView<R1Tensor const> const & cellCenterToEdgeCenters = stencil.getCellCenterToEdgeCenters();

  forall_in_range<serialPolicy>( 0, stencil.size(), GEOSX_LAMBDA ( localIndex iconn )
  {

    localIndex const numFluxElems = stencil.stencilSize(iconn);

    localIndex const er = seri[iconn][0];
    localIndex const esr = sesri[iconn][0];

    ProppantPackVolumeKernel::UpdateInterfaceElement( numFluxElems,
                                                      sei[iconn],
                                                      weights[iconn],
                                                      cellCenterToEdgeCenters[iconn],
                                                      unitGravityVector,                                                      
                                                      isProppantMobile[er][esr],
                                                      isInterfaceElement[er][esr]);

  } );
  
}

  
} // namespace ProppantTransportKernels

} // namespace geosx
