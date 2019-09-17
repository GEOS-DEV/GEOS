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
 * @file ProppantTransportKernels.hpp
 */

#ifndef GEOSX_PROPPANTTRANSPORTKERNELS_HPP
#define GEOSX_PROPPANTTRANSPORTKERNELS_HPP

#include "common/DataTypes.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

namespace geosx
{

namespace ProppantTransportKernels
{

/******************************** AccumulationKernel ********************************/

struct AccumulationKernel
{

  inline static void
  Compute( real64 const & densNew,
           real64 const & densOld,
           real64 const & dDens_dPres,
           real64 const & dDens_dConc,
	   real64 const & concNew,
           real64 const & concOld,	   
           real64 const & volume,
           arraySlice1d<real64> const & localAccum,
           arraySlice2d<real64> const & localAccumJacobian )
  {
  
        // fluid-mixture mass conservation
        localAccum[0] = (densNew  - densOld) * volume;

        // proppant mass conservation
        localAccum[1] = (concNew - concOld) * volume;

        // Derivative of residual wrt to pressure and concentration in the cell
        localAccumJacobian[0][0] = dDens_dPres * volume;
        localAccumJacobian[0][1] = dDens_dConc * volume;

        localAccumJacobian[1][0] = 0.0;
        localAccumJacobian[1][1] = volume;

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
  using ElementView = typename ElementRegionManager::ElementViewAccessor<VIEWTYPE>::ViewTypeConst;

  /**
   * @brief The type for element-based constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::MaterialViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using MaterialView = typename ElementRegionManager::MaterialViewAccessor<VIEWTYPE>::ViewTypeConst;

  /**
   * @brief launches the kernel to assemble the flux contributions to the linear system.
   * @tparam STENCIL_TYPE The type of the stencil that is being used.
   * @param[in] stencil The stencil object.
   * @param[in] dt The timestep for the integration step.
   * @param[in] fluidIndex The index of the fluid being fluxed.
   * @param[in] gravityFlag Flag to indicate whether or not to use gravity.
   * @param[in] dofNumber The dofNumbers for each element
   * @param[in] pres The pressures in each element
   * @param[in] dPres The change in pressure for each element
   * @param[in] gravDepth The factor for gravity calculations (g*H)
   * @param[in] dens The material density in each element
   * @param[in] dDens_dPres The change in material density for each element
   * @param[in] mob The fluid mobility in each element
   * @param[in] dMob_dPres The derivative of mobility wrt pressure in each element
   * @param[out] jacobian The linear system matrix
   * @param[out] residual The linear system residual
   */
  template< typename STENCIL_TYPE >
  static void
  Launch( STENCIL_TYPE const & stencil,
	  localIndex const numDofPerCell,	  
          real64 const dt,
          localIndex const fluidIndex,
          localIndex const proppantIndex,
	  bool updateProppantMobilityFlag,
	  bool updatePermeabilityFlag,
          ElementView < arrayView1d<globalIndex const > > const & dofNumber,
          ElementView < arrayView1d<real64 const> > const & pres,
          ElementView < arrayView1d<real64 const> > const & dPres,
          ElementView < arrayView1d<real64 const> > const & Conc,
          ElementView < arrayView1d<real64 const> > const & ConcOld,	  
          ElementView < arrayView1d<real64 const> > const & dConc,
	  ElementView < arrayView1d<real64 const> > const & gravDepth,
          MaterialView< arrayView2d<real64 const> > const & dens,
          MaterialView< arrayView2d<real64 const> > const & dDens_dPres,
          MaterialView< arrayView2d<real64 const> > const & dDens_dConc,
	  MaterialView< arrayView2d<real64 const> > const & visc,
          MaterialView< arrayView2d<real64 const> > const & dVisc_dPres,
          MaterialView< arrayView2d<real64 const> > const & dVisc_dConc,
	  MaterialView< arrayView2d<real64 const> > const & fluidDensity,
          MaterialView< arrayView2d<real64 const> > const & dFluidDens_dPres,	  
          MaterialView< arrayView1d<real64 const> > const & settlingFactor,
          MaterialView< arrayView1d<real64 const> > const & dSettlingFactor_dConc,
          MaterialView< arrayView1d<real64 const> > const & collisionFactor,
          MaterialView< arrayView1d<real64 const> > const & dCollisionFactor_dConc,
          MaterialView< arrayView1d<bool const> > const & isProppantMobile,	  
          MaterialView< arrayView1d<real64 const> > const & proppantPackPermeability,
          ElementView < arrayView1d<real64 const> > const & aperture0,
          ElementView < arrayView1d<real64 const> > const & aperture,
          ParallelMatrix * const jacobian,
          ParallelVector * const residual );


  /**
     * @brief Compute flux and its derivatives for a given multi-element connector.
     *
     * This is a specialized version that flux in a single region, and uses
     * element pairing instead of a proper junction.
     */
  inline static void
  ComputeJunction( localIndex const numElems,
		   localIndex const numDofPerCell,		   
                   arraySlice1d<localIndex const> const & stencilElementIndices,
                   arraySlice1d<real64 const> const & stencilWeights,
                   arraySlice1d<real64 const> const & stencilEdgeToFaceDownDistances,
		   arrayView1d<real64 const> const & pres,
                   arrayView1d<real64 const> const & dPres,
		   arrayView1d<real64 const> const & conc,
		   arrayView1d<real64 const> const & concOld,		   
                   arrayView1d<real64 const> const & dConc,		   
		   arrayView1d<real64 const> const & gravDepth,
                   arrayView2d<real64 const> const & dens,
                   arrayView2d<real64 const> const & dDens_dPres,
                   arrayView2d<real64 const> const & dDens_dConc,
                   arrayView2d<real64 const> const & visc,
                   arrayView2d<real64 const> const & dVisc_dPres,
                   arrayView2d<real64 const> const & dVisc_dConc,		   
                   arrayView2d<real64 const> const & fluidDensity,
                   arrayView2d<real64 const> const & dFluidDens_dPres,
                   arrayView1d<real64 const> const & settlingFactor,
                   arrayView1d<real64 const> const & dSettlingFactor_dConc,
                   arrayView1d<real64 const> const & collisionFactor,
                   arrayView1d<real64 const> const & dCollisionFactor_dConc,
                   arrayView1d<bool const> const & isProppantMobile,
                   arrayView1d<real64 const> const & proppantPackPermeability,		   
                   arrayView1d<real64 const> const &,
                   arrayView1d<real64 const> const & aperture,
		   bool updateProppantMobilityFlag,
		   bool updatePermeabilityFlag,
                   real64 const dt,
                   arraySlice1d<real64> const & localFlux,
                   arraySlice2d<real64> const & localFluxJacobian)
  {

    // We assume numElems == stencilSize;

    // working array
    constexpr localIndex maxNumFluxElems = FaceElementStencil::NUM_POINT_IN_FLUX;
    
    stackArray1d<real64, maxNumFluxElems> weight(numElems);

    // mixture density and fluid density in each face
    stackArray1d<real64, maxNumFluxElems> mixDens(numElems);
    stackArray1d<real64, maxNumFluxElems> dMixDens_dP(numElems);
    stackArray1d<real64, maxNumFluxElems> dMixDens_dC(numElems);

    stackArray1d<real64, maxNumFluxElems> fluidDens(numElems);
    stackArray1d<real64, maxNumFluxElems> dFluidDens_dP(numElems);


    // realted to slip velocity calculation
    stackArray1d<real64, maxNumFluxElems> settlingFac(numElems);
    stackArray1d<real64, maxNumFluxElems> dSettlingFac_dC(numElems);

    stackArray1d<real64, maxNumFluxElems> collisionFac(numElems);
    stackArray1d<real64, maxNumFluxElems> dCollisionFac_dC(numElems);

    stackArray1d<real64, maxNumFluxElems> transT(numElems);
    stackArray1d<real64, maxNumFluxElems> coefs(numElems);
    stackArray1d<bool, maxNumFluxElems> isProppantMob(numElems);

    real64 edgeDensity, edgeViscosity;
    stackArray1d<real64, maxNumFluxElems> dEdgeDens_dP(numElems);
    stackArray1d<real64, maxNumFluxElems> dEdgeDens_dC(numElems);

    stackArray1d<real64, maxNumFluxElems> dEdgeVisc_dP(numElems);
    stackArray1d<real64, maxNumFluxElems> dEdgeVisc_dC(numElems);

    stackArray1d<real64, maxNumFluxElems> dPe_dP(numElems);
    stackArray1d<real64, maxNumFluxElems> dPe_dC(numElems);
    stackArray1d<R1Tensor, maxNumFluxElems> Up(numElems);

    stackArray1d<real64, maxNumFluxElems> dCe_dP(numElems);
    stackArray1d<real64, maxNumFluxElems> dCe_dC(numElems);

    stackArray1d<real64, maxNumFluxElems> edgeToFaceFlux(numElems);
    stackArray2d<real64, maxNumFluxElems * maxNumFluxElems> dEdgeToFaceFlux_dP(numElems,numElems);
    stackArray2d<real64, maxNumFluxElems * maxNumFluxElems> dEdgeToFaceFlux_dC(numElems,numElems);

    stackArray1d<real64, maxNumFluxElems> edgeToFaceProppantFlux(numElems);
    stackArray2d<real64, maxNumFluxElems * maxNumFluxElems> dEdgeToFaceProppantFlux_dP(numElems,numElems);
    stackArray2d<real64, maxNumFluxElems * maxNumFluxElems> dEdgeToFaceProppantFlux_dC(numElems,numElems);

    stackArray1d<real64, maxNumFluxElems> P(numElems);
    stackArray1d<real64, maxNumFluxElems> C(numElems);

    // clear working arrays

    edgeDensity = 0.0;
    edgeViscosity = 0.0;

    dEdgeDens_dP = 0.0;
    dEdgeDens_dC = 0.0;

    dEdgeVisc_dP = 0.0;
    dEdgeVisc_dC = 0.0;

    dEdgeToFaceProppantFlux_dP = 0.0;
    dEdgeToFaceProppantFlux_dC = 0.0;

    real64 aperTerm;
    real64 sumOfWeights = 0;

    for( localIndex i=0 ; i<numElems ; ++i )
    {

      localIndex const ei  = stencilElementIndices[i];      


      aperTerm = aperture[ei] * aperture[ei] * aperture[ei];
      
      sumOfWeights += stencilWeights[i];
      weight[i] = stencilWeights[i];

      transT[i] = aperTerm * stencilWeights[i];
 	
    }

    for( localIndex i=0 ; i<numElems ; ++i )
      weight[i] /= sumOfWeights;

    //get averaged edgeDensity and edgeViscosity

    for (localIndex i = 0; i < numElems; ++i) {

      localIndex const ei  = stencilElementIndices[i];

      edgeDensity += weight[i] * dens[ei][0];
      dEdgeDens_dP[i] = weight[i] * dDens_dPres[ei][0];
      dEdgeDens_dC[i] = weight[i] * dDens_dConc[ei][0];

      edgeViscosity += weight[i] * visc[ei][0];
      dEdgeVisc_dP[i] = weight[i] * dVisc_dPres[ei][0];
      dEdgeVisc_dC[i] = weight[i] * dVisc_dConc[ei][0];

      P[i] = pres[ei] + dPres[ei];
      C[i] = conc[ei] + dConc[ei];

      mixDens[i] = dens[ei][0];
      dMixDens_dP[i] = dDens_dPres[ei][0];
      dMixDens_dC[i] = dDens_dConc[ei][0];

      fluidDens[i] = fluidDensity[ei][0];
      dFluidDens_dP[i] = dFluidDens_dPres[ei][0];

      settlingFac[i] = settlingFactor[ei];
      dSettlingFac_dC[i] = dSettlingFactor_dConc[ei];

      collisionFac[i] = collisionFactor[ei];
      dCollisionFac_dC[i] = dCollisionFactor_dConc[ei];

      isProppantMob[i] = isProppantMobile[ei];

      if(proppantPackPermeability[ei] > 0.0 && updatePermeabilityFlag)
      {
        transT[i] = transT[i] * (1.0 - concOld[ei]) + proppantPackPermeability[ei]  * 12 * aperture[ei] * stencilWeights[i] * concOld[ei];

      }

      coefs[i] =  stencilEdgeToFaceDownDistances[i] * aperture[ei];
      
    }

    real64 transTSum = 0.0;
    real64 Pe = 0.0;
    dPe_dP = 0.0;
    dPe_dC = 0.0;

    for (localIndex i = 0; i < numElems; ++i)
    {

      localIndex const ei  = stencilElementIndices[i];
      
      real64 const gravD    = gravDepth[ei];
      real64 const gravTerm = edgeDensity * gravD;

      Pe += transT[i] * (P[i] - gravTerm);

      transTSum += transT[i];

      dPe_dP[i] += transT[i];

      for (localIndex j = 0; j < numElems; ++j) {
        dPe_dP[j] += -transT[i] * gravD * dEdgeDens_dP[j];
        dPe_dC[j] += -transT[i] * gravD * dEdgeDens_dC[j];
      }

    }

    for (localIndex i = 0; i < numElems; ++i)
    {

      dPe_dP[i] /= transTSum;
      dPe_dC[i] /= transTSum;

    }

    Pe /= transTSum;

    for (localIndex i = 0; i < numElems; ++i)
    {

      localIndex const ei  = stencilElementIndices[i];      

      real64 const gravD    = gravDepth[ei];
      real64 const gravTerm = edgeDensity * gravD;

      real64 const fluxTerm = Pe - (P[i] - gravTerm);
      
      edgeToFaceFlux[i] = transT[i] * fluxTerm / edgeViscosity;

      dEdgeToFaceFlux_dP[i][i] += -transT[i] / edgeViscosity;

      for (localIndex j = 0; j < numElems; ++j)
      {

        dEdgeToFaceFlux_dP[i][j] += -transT[i] * fluxTerm * dEdgeVisc_dP[j] / (edgeViscosity * edgeViscosity) + transT[i] * (dPe_dP[j] + dEdgeDens_dP[j] * gravD) / edgeViscosity;

        dEdgeToFaceFlux_dC[i][j] += -transT[i] * fluxTerm * dEdgeVisc_dC[j] / (edgeViscosity * edgeViscosity) + transT[i] * (dPe_dC[j] + dEdgeDens_dC[j] * gravD) / edgeViscosity;


      }

      //contribution from paricle slipping

      edgeToFaceProppantFlux[i] = (1.0 + fluidDens[i] / mixDens[i] * (1.0 - C[i]) * collisionFac[i]) *  edgeToFaceFlux[i] + fluidDens[i] / mixDens[i] * (1.0 - C[i]) * settlingFac[i] * coefs[i];

      dEdgeToFaceProppantFlux_dP[i][i] = (dFluidDens_dP[i] / mixDens[i] - fluidDens[i] * dMixDens_dP[i] / (mixDens[i] * mixDens[i])) * (1.0 - C[i]) * (collisionFac[i] * edgeToFaceFlux[i] + settlingFac[i] * coefs[i]);

      dEdgeToFaceProppantFlux_dC[i][i] = -fluidDens[i] / mixDens[i] * (1.0 + dMixDens_dC[i] / mixDens[i] * (1.0 - C[i])) * (collisionFac[i] *  edgeToFaceFlux[i] + settlingFac[i] * coefs[i]) + fluidDens[i] / mixDens[i] * (1.0 - C[i]) * (dCollisionFac_dC[i] *  edgeToFaceFlux[i] + dSettlingFac_dC[i] * coefs[i]);

      for (localIndex j = 0; j < numElems; ++j)
      {

        dEdgeToFaceProppantFlux_dP[i][j] += (1.0 + fluidDens[i] / mixDens[i] * (1.0 - C[i]) * collisionFac[i]) *  dEdgeToFaceFlux_dP[i][j];

        dEdgeToFaceProppantFlux_dC[i][j] += (1.0 + fluidDens[i] / mixDens[i] * (1.0 - C[i]) * collisionFac[i]) *  dEdgeToFaceFlux_dC[i][j];
      }

    }

    // get Ce

    real64 Ce = 0.0;
    dCe_dP = 0.0;
    dCe_dC = 0.0;

    real64 downStreamFlux = 0.0;
    stackArray1d<real64, maxNumFluxElems> dDownStreamFlux_dP(numElems);
    stackArray1d<real64, maxNumFluxElems> dDownStreamFlux_dC(numElems);

    dDownStreamFlux_dP = 0.0;
    dDownStreamFlux_dC = 0.0;

    for (localIndex i = 0; i < numElems; ++i)
    {

      if(!isProppantMob[i] && updateProppantMobilityFlag)
        continue;

      if(edgeToFaceProppantFlux[i] >= 0.0)
      {
        // downstream
        downStreamFlux += edgeToFaceProppantFlux[i];

        for(localIndex j = 0; j < numElems; ++j)
        {

          dDownStreamFlux_dP[j] += dEdgeToFaceProppantFlux_dP[i][j];
          dDownStreamFlux_dC[j] += dEdgeToFaceProppantFlux_dC[i][j];

        }

      }
      else
      {

        // upstream
        Ce += -edgeToFaceProppantFlux[i] * C[i];

        dCe_dC[i] += -edgeToFaceProppantFlux[i];

        for(localIndex j = 0; j < numElems; ++j)
        {

          dCe_dP[j] += -dEdgeToFaceProppantFlux_dP[i][j] * C[i];
          dCe_dC[j] += -dEdgeToFaceProppantFlux_dC[i][j] * C[i];

        }

      }

    }

  
    if(downStreamFlux > 0.0)
    {

      for (localIndex i = 0; i < numElems; ++i)
      {

        if(!isProppantMob[i] && updateProppantMobilityFlag)
          continue;

        dCe_dP[i] =  dCe_dP[i] / downStreamFlux - Ce * dDownStreamFlux_dP[i] / (downStreamFlux * downStreamFlux);
        dCe_dC[i] =  dCe_dC[i] / downStreamFlux - Ce * dDownStreamFlux_dC[i] / (downStreamFlux * downStreamFlux);;

      }

      Ce = Ce / downStreamFlux;

    }
    else
    {

      Ce = 0.0;
      for (localIndex i = 0; i < numElems; ++i)
      {

        if(!isProppantMob[i] && updateProppantMobilityFlag)
          continue;

        dCe_dP[i] =  0.0;
        dCe_dC[i] =  weight[i];
        Ce += C[i] * weight[i];

      }

    }

    for (localIndex i = 0; i < numElems; ++i)
    {

      localIndex idx1 = i * numDofPerCell;

      localFlux[idx1] = -edgeDensity * edgeToFaceFlux[i] * dt;

      for (localIndex j = 0; j < numElems; ++j)
      {

        localIndex idx2 = j * numDofPerCell;

        //dP

        localFluxJacobian[idx1][idx2]  = -(dEdgeDens_dP[j] * edgeToFaceFlux[i] + edgeDensity * dEdgeToFaceFlux_dP[i][j]) * dt;

        //dC

        localFluxJacobian[idx1][idx2 + 1]  = -(dEdgeDens_dC[j] * edgeToFaceFlux[i] + edgeDensity * dEdgeToFaceFlux_dC[i][j]) * dt;

      }


      if(isProppantMob[i] || !updateProppantMobilityFlag)
      {

        if(edgeToFaceProppantFlux[i] >= 0.0)
        {


          localFlux[idx1 + 1] = -Ce * edgeToFaceProppantFlux[i] * dt;
        }
        else
        {

          localFlux[idx1 + 1] = -C[i] * edgeToFaceProppantFlux[i] * dt;

        }

        for (localIndex j = 0; j < numElems; ++j)
        {

          localIndex idx2 = j * numDofPerCell;


          if(edgeToFaceProppantFlux[i] >= 0.0)
          {


            localFluxJacobian[idx1 + 1][idx2] = -(dCe_dP[j] * edgeToFaceProppantFlux[i] + Ce * dEdgeToFaceProppantFlux_dP[i][j]) * dt;

            localFluxJacobian[idx1 + 1][idx2 + 1] = -(dCe_dC[j] * edgeToFaceProppantFlux[i] + Ce * dEdgeToFaceProppantFlux_dC[i][j]) * dt;

          }
          else
          {


            localFluxJacobian[idx1 + 1][idx2] = -C[i] * dEdgeToFaceProppantFlux_dP[i][j] * dt;
            localFluxJacobian[idx1 + 1][idx2 + 1] = -C[i] * dEdgeToFaceProppantFlux_dC[i][j] * dt;

            if(i == j)
              localFluxJacobian[idx1 + 1][idx2 + 1] += -edgeToFaceProppantFlux[i] * dt;

          }

        }

      }
      else
      {

        localFlux[idx1 + 1] = 0.0;
        for (localIndex j = 0; j < numElems; ++j)
        {

          localIndex idx2 = j * numDofPerCell;
          localFluxJacobian[idx1 + 1][idx2] = 0.0;
          localFluxJacobian[idx1 + 1][idx2 + 1] = 0.0;
        }
      }
    }

  }

};


} // namespace ProppantTransportKernels

} // namespace geosx

#endif //GEOSX_PROPPANTTRANSPORTKERNELS_HPP
