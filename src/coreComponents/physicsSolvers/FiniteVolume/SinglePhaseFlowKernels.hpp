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
 * @file SinglePhaseFlowKernels.hpp
 */

#ifndef GEOSX_SINGLEPHASEFLOWKERNELS_HPP
#define GEOSX_SINGLEPHASEFLOWKERNELS_HPP

#include "common/DataTypes.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"


class Epetra_FECrsMatrix;
class Epetra_FEVector;


namespace geosx
{


namespace SinglePhaseFlowKernels
{

/******************************** MobilityKernel ********************************/

struct MobilityKernel
{
  static void
  Compute( real64 const & dens,
           real64 const & dDens_dPres,
           real64 const & visc,
           real64 const & dVisc_dPres,
           real64 & mob,
           real64 & dMob_dPres );

  static void
  Compute( real64 const & dens,
           real64 const & visc,
           real64 & mob );

  static void Launch( localIndex begin, localIndex end,
                      arrayView2d<real64 const> const & dens,
                      arrayView2d<real64 const> const & dDens_dPres,
                      arrayView2d<real64 const> const & visc,
                      arrayView2d<real64 const> const & dVisc_dPres,
                      arrayView1d<real64> const & mob,
                      arrayView1d<real64> const & dMob_dPres );

  static void Launch( set<localIndex> targetSet,
                      arrayView2d<real64 const> const & dens,
                      arrayView2d<real64 const> const & dDens_dPres,
                      arrayView2d<real64 const> const & visc,
                      arrayView2d<real64 const> const & dVisc_dPres,
                      arrayView1d<real64> const & mob,
                      arrayView1d<real64> const & dMob_dPres );

  static void Launch( localIndex begin, localIndex end,
                      arrayView2d<real64 const> const & dens,
                      arrayView2d<real64 const> const & visc,
                      arrayView1d<real64> const & mob );

  static void Launch( set<localIndex> targetSet,
                      arrayView2d<real64 const> const & dens,
                      arrayView2d<real64 const> const & visc,
                      arrayView1d<real64> const & mob );
};

/******************************** AccumulationKernel ********************************/

template<bool ISPORO>
struct AssembleAccumulationTermsHelper;

template<>
struct AssembleAccumulationTermsHelper<true>
{
  inline static constexpr void
  porosityUpdate( real64 & poro,
                  real64 & dPoro_dPres,
                  real64 const biotCoefficient,
                  real64 const poroOld,
                  real64 const bulkModulus,
                  real64 const totalMeanStress,
                  real64 const oldTotalMeanStress,
                  real64 const dPres,
                  real64 const poroRef,
                  real64 const pvmult,
                  real64 const dPVMult_dPres )
  {
    dPoro_dPres = (biotCoefficient - poroOld) / bulkModulus;
    poro = poroOld + dPoro_dPres * (totalMeanStress - oldTotalMeanStress + dPres);
  }
};

template<>
struct AssembleAccumulationTermsHelper<false>
{
  inline static constexpr void
  porosityUpdate( real64 & poro,
                  real64 & dPoro_dPres,
                  real64 const biotCoefficient,
                  real64 const poroOld,
                  real64 const bulkModulus,
                  real64 const totalMeanStress,
                  real64 const oldTotalMeanStress,
                  real64 const dPres,
                  real64 const poroRef,
                  real64 const pvmult,
                  real64 const dPVMult_dPres )
  {
    poro = poroRef * pvmult;
    dPoro_dPres = dPVMult_dPres * poroRef;
  }
};

struct AccumulationKernel
{

  template<bool COUPLED>
  inline static void
  Compute( real64 const & dPres,
           real64 const & densNew,
           real64 const & densOld,
           real64 const & dDens_dPres,
           real64 const & volume,
           real64 const & dVol,
           real64 const & poroRef,
           real64 const & poroOld,
           real64 const & pvMult,
           real64 const & dPVMult_dPres,
           real64 const & biotCoefficient,
           real64 const & bulkModulus,
           real64 const & totalMeanStress,
           real64 const & oldTotalMeanStress,
           real64 & poroNew,
           real64 & localAccum,
           real64 & localAccumJacobian )
  {
    real64 const volNew = volume + dVol;

    // TODO porosity update needs to be elsewhere...
    real64 dPoro_dPres;
    AssembleAccumulationTermsHelper<COUPLED>::porosityUpdate( poroNew,
                                                              dPoro_dPres,
                                                              biotCoefficient,
                                                              poroOld,
                                                              bulkModulus,
                                                              totalMeanStress,
                                                              oldTotalMeanStress,
                                                              dPres,
                                                              poroRef,
                                                              pvMult,
                                                              dPVMult_dPres );


    // Residual contribution is mass conservation in the cell
    localAccum = poroNew * densNew * volNew - poroOld * densOld * volume;
    //localAccum = poroOld * densNew * volNew - poroOld * densOld * volume;

    // Derivative of residual wrt to pressure in the cell
    localAccumJacobian = (dPoro_dPres * densNew + dDens_dPres * poroNew) * volNew;
    //localAccumJacobian = (0 * densNew + dDens_dPres * poroOld) * volNew;
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
          real64 const dt,
          localIndex const fluidIndex,
          integer const gravityFlag,
          ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex> > const & dofNumber,
          FluxKernel::ElementView < arrayView1d<real64 const> > const & pres,
          FluxKernel::ElementView < arrayView1d<real64 const> > const & dPres,
          FluxKernel::ElementView < arrayView1d<real64 const> > const & gravDepth,
          FluxKernel::MaterialView< arrayView2d<real64 const> > const & dens,
          FluxKernel::MaterialView< arrayView2d<real64 const> > const & dDens_dPres,
          FluxKernel::ElementView < arrayView1d<real64 const> > const & mob,
          FluxKernel::ElementView < arrayView1d<real64 const> > const & dMob_dPres,
          Epetra_FECrsMatrix * const jacobian,
          Epetra_FEVector * const residual );


  /**
   * @brief Compute flux and its derivatives for a given connection
   *
   * This is a general version that assumes different element regions.
   * See below for a specialized version for fluxes within a region.
   */
  inline static void
  Compute( localIndex const stencilSize,
           arraySlice1d<localIndex const> const & seri,
           arraySlice1d<localIndex const> const & sesri,
           arraySlice1d<localIndex const> const & sei,
           arraySlice1d<real64 const> const & stencilWeights,
           ElementView <arrayView1d<real64 const>> const & pres,
           ElementView <arrayView1d<real64 const>> const & dPres,
           ElementView <arrayView1d<real64 const>> const & gravDepth,
           MaterialView<arrayView2d<real64 const>> const & dens,
           MaterialView<arrayView2d<real64 const>> const & dDens_dPres,
           ElementView <arrayView1d<real64 const>> const & mob,
           ElementView <arrayView1d<real64 const>> const & dMob_dPres,
           localIndex const fluidIndex,
           integer const gravityFlag,
           real64 const dt,
           arraySlice1d<real64> const & flux,
           arraySlice2d<real64> const & fluxJacobian )
  {
    localIndex constexpr numElems = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
    localIndex constexpr maxStencil = CellElementStencilTPFA::MAX_STENCIL_SIZE;

    stackArray1d<real64, numElems>   densWeight(numElems);
    stackArray1d<real64, maxStencil> dDensMean_dP(stencilSize);
    stackArray1d<real64, maxStencil> dFlux_dP(stencilSize);

    // clear working arrays
    dDensMean_dP = 0.0;

    // density averaging weights
    densWeight = 0.5;

    // calculate quantities on primary connected cells
    real64 densMean = 0.0;
    for (localIndex ke = 0; ke < numElems; ++ke)
    {
      // density
      real64 const density = dens[seri[ke]][sesri[ke]][fluidIndex][sei[ke]][0];
      real64 const dDens_dP = dDens_dPres[seri[ke]][sesri[ke]][fluidIndex][sei[ke]][0];

      // average density
      densMean        += densWeight[ke] * density;
      dDensMean_dP[ke] = densWeight[ke] * dDens_dP;
    }

    // compute potential difference MPFA-style
    real64 potDif = 0.0;
    for (localIndex ke = 0; ke < stencilSize; ++ke)
    {
      localIndex const er  = seri[ke];
      localIndex const esr = sesri[ke];
      localIndex const ei  = sei[ke];

      real64 weight = stencilWeights[ke];

      real64 const gravD = gravDepth[er][esr][ei];
      real64 const gravTerm = gravityFlag ? densMean * gravD : 0.0;
      real64 const dGrav_dP = gravityFlag ? dDensMean_dP[ke] * gravD : 0.0;

      potDif += weight * (pres[er][esr][ei] + dPres[er][esr][ei] - gravTerm);
      dFlux_dP[ke] = weight * (1.0 - dGrav_dP);
    }

    // upwinding of fluid properties (make this an option?)
    localIndex const k_up = (potDif >= 0) ? 0 : 1;

    localIndex er_up  = seri[k_up];
    localIndex esr_up = sesri[k_up];
    localIndex ei_up  = sei[k_up];

    real64 const mobility     = mob[er_up][esr_up][ei_up];
    real64 const dMobility_dP = dMob_dPres[er_up][esr_up][ei_up];

    // compute the final flux and derivatives
    real64 const fluxVal = mobility * potDif;
    for (localIndex ke = 0; ke < stencilSize; ++ke)
    {
      dFlux_dP[ke] *= mobility;
    }

    dFlux_dP[k_up] += dMobility_dP * potDif;

    // populate local flux vector and derivatives
    flux[0] = dt * fluxVal;
    flux[1] = -flux[0];

    for (localIndex ke = 0; ke < stencilSize; ++ke)
    {
      fluxJacobian[0][ke] = dt * dFlux_dP[ke];
      fluxJacobian[1][ke] = -fluxJacobian[0][ke];
    }
  }


  /**
   * @brief Compute flux and its derivatives for a given connection
   *.
   * This is a specialized version for fluxes within the same region.
   * See above for a general version.
   */
  inline static void
  ComputeCellTPFA( localIndex const stencilSize,
                   arraySlice1d<localIndex> const & stencilElementIndices,
                   arraySlice1d<real64> const & stencilWeights,
                   arrayView1d<real64 const> const & pres,
                   arrayView1d<real64 const> const & dPres,
                   arrayView1d<real64 const> const & gravDepth,
                   arrayView2d<real64 const> const & dens,
                   arrayView2d<real64 const> const & dDens_dPres,
                   arrayView1d<real64 const> const & mob,
                   arrayView1d<real64 const> const & dMob_dPres,
                   localIndex const fluidIndex,
                   integer const gravityFlag,
                   real64 const dt,
                   arraySlice1d<real64> const & flux,
                   arraySlice2d<real64> const & fluxJacobian )
  {
    localIndex constexpr numElems = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
    localIndex constexpr maxStencil = CellElementStencilTPFA::MAX_STENCIL_SIZE;

    stackArray1d<real64, numElems> densWeight(numElems);
    stackArray1d<real64, maxStencil> dDensMean_dP(stencilSize);
    stackArray1d<real64, maxStencil> dFlux_dP(stencilSize);

    // clear working arrays
    dDensMean_dP = 0.0;

    // density averaging weights
    densWeight = 1.0 / numElems;

    // calculate quantities on primary connected cells
    real64 densMean = 0.0;
    for (localIndex i = 0; i < numElems; ++i)
    {
      // density
      real64 const density = dens[stencilElementIndices[i]][0];
      real64 const dDens_dP = dDens_dPres[stencilElementIndices[i]][0];

      // average density
      densMean += densWeight[i] * density;
      dDensMean_dP[i] = densWeight[i] * dDens_dP;
    }

    // compute potential difference MPFA-style
    real64 potDif = 0.0;
    for (localIndex ke = 0; ke < stencilSize; ++ke)
    {
      localIndex const ei = stencilElementIndices[ke];
      real64 const weight = stencilWeights[ke];

      real64 const gravD = gravDepth[ei];
      real64 const gravTerm = gravityFlag ? densMean * gravD : 0.0;
      real64 const dGrav_dP = gravityFlag ? dDensMean_dP[ke] * gravD : 0.0;

      potDif += weight * (pres[ei] + dPres[ei] - gravTerm);
      dFlux_dP[ke] = weight * (1.0 - dGrav_dP);
    }

    // upwinding of fluid properties (make this an option?)
    localIndex const k_up = (potDif >= 0) ? 0 : 1;

    localIndex ei_up  = stencilElementIndices[k_up];

    real64 const mobility     = mob[ei_up];
    real64 const dMobility_dP = dMob_dPres[ei_up];

    // compute the final flux and derivatives
    real64 const fluxVal = mobility * potDif;
    for (localIndex ke = 0; ke < stencilSize; ++ke)
    {
      dFlux_dP[ke] *= mobility;
    }

    dFlux_dP[k_up] += dMobility_dP * potDif;

    // populate local flux vector and derivatives
    flux[0] = dt * fluxVal;
    flux[1] = -flux[0];

    for (localIndex ke = 0; ke < stencilSize; ++ke)
    {
      fluxJacobian[0][ke] = dt * dFlux_dP[ke];
      fluxJacobian[1][ke] = -fluxJacobian[0][ke];
    }
  }

  /**
     * @brief Compute flux and its derivatives for a given multi-element connector.
     *
     * This is a specialized version that flux in a single region, and uses
     * element pairing instead of a proper junction.
     */
  inline static void
  ComputeJunction( localIndex const numFluxElems,
                   arraySlice1d<localIndex const> const & stencilElementIndices,
                   arraySlice1d<real64 const> const & stencilWeights,
                   arrayView1d<real64 const> const & pres,
                   arrayView1d<real64 const> const & dPres,
                   arrayView1d<real64 const> const & gravDepth,
                   arrayView2d<real64 const> const & dens,
                   arrayView2d<real64 const> const & dDens_dPres,
                   arrayView1d<real64 const> const & mob,
                   arrayView1d<real64 const> const & dMob_dPres,
                   localIndex const fluidIndex,
                   integer const gravityFlag,
                   real64 const dt,
                   arraySlice1d<real64> const & flux,
                   arraySlice2d<real64> const & fluxJacobian )
  {
    real64 sumOfWeights = 0;
    for( localIndex k=0 ; k<numFluxElems ; ++k )
    {
      sumOfWeights += stencilWeights[k];
    }

    localIndex k[2];
    for( k[0]=0 ; k[0]<numFluxElems ; ++k[0] )
    {
      for( k[1]=k[0]+1 ; k[1]<numFluxElems ; ++k[1] )
      {
        real64  dFlux_dP[2] = {0,0};

        localIndex const ei[2] = { stencilElementIndices[k[0]],
                                   stencilElementIndices[k[1]] };

        real64 const weight[2] = {   stencilWeights[k[0]] * stencilWeights[k[1]] / sumOfWeights,
                                     - stencilWeights[k[0]] * stencilWeights[k[1]] / sumOfWeights };


        // average density
        real64 const densMean = 0.5 * ( dens[ei[0]][0] + dens[ei[1]][0] );

        real64 const dDensMean_dP[2] = { 0.5 * dDens_dPres[ei[0]][0],
                                         0.5 * dDens_dPres[ei[1]][0] };

        real64 potDif = 0.0;
        for( localIndex i = 0 ; i < 2 ; ++i )
        {
          real64 const gravD = gravDepth[ei[i]];
          real64 const gravTerm = gravityFlag ? densMean * gravD : 0.0;
          real64 const dGrav_dP = gravityFlag ? dDensMean_dP[i] * gravD : 0.0;

          potDif += weight[i] * (pres[ei[i]] + dPres[ei[i]] - gravTerm);
          dFlux_dP[i] = weight[i] * (1.0 - dGrav_dP);
        }

        // upwinding of fluid properties (make this an option?)
        localIndex const k_up = (potDif >= 0) ? 0 : 1;

        localIndex ei_up  = stencilElementIndices[k[k_up]];

        real64 const mobility     = mob[ei_up];
        real64 const dMobility_dP = dMob_dPres[ei_up];

        // compute the final flux and derivatives
        real64 const fluxVal = mobility * potDif;
        dFlux_dP[0] *= mobility;
        dFlux_dP[1] *= mobility;

        dFlux_dP[k_up] += dMobility_dP * potDif;

        // populate local flux vector and derivatives
        flux[k[0]] += dt * fluxVal;
        flux[k[1]] -= dt * fluxVal;

        fluxJacobian[k[0]][k[0]] += dt * dFlux_dP[0];
        fluxJacobian[k[1]][k[0]] -= dt * dFlux_dP[0];
        fluxJacobian[k[0]][k[1]] += dt * dFlux_dP[1];
        fluxJacobian[k[1]][k[1]] -= dt * dFlux_dP[1];

      }
    }
  }

};

} // namespace SinglePhaseFlowKernels

} // namespace geosx

#endif //GEOSX_SINGLEPHASEFLOWKERNELS_HPP
