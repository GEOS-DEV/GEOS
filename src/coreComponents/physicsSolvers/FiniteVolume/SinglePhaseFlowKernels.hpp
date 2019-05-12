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

namespace geosx
{


namespace SinglePhaseFlowKernels
{

/******************************** MobilityKernel ********************************/

struct MobilityKernel
{
  static inline RAJA_HOST_DEVICE void
  Compute( real64 const & dens,
           real64 const & dDens_dPres,
           real64 const & visc,
           real64 const & dVisc_dPres,
           real64 & mob,
           real64 & dMob_dPres );

  static inline RAJA_HOST_DEVICE void
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

    // Derivative of residual wrt to pressure in the cell
    localAccumJacobian = (dPoro_dPres * densNew + dDens_dPres * poroNew) * volNew;
  }
};

/******************************** FluxKernel ********************************/

/**
 * @brief Compute flux and its derivatives for a given connection represented by a stencil object
 * @param stencil
 * @param pres
 * @param dPres
 * @param gravDepth
 * @param dens
 * @param dDens_dPres
 * @param visc
 * @param dVisc_dPres
 * @param fluidIndex
 * @param gravityFlag
 * @param dt
 * @param flux
 * @param fluxJacobian
 *
 * This is a general version that assumes different element regions.
 * See below for a specialized version for fluxes within a region.
 */
inline static void
MakeFlux( localIndex const stencilSize,
          FluxApproximationBase::CellStencil::Entry const * stencil,
          ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & pres,
          ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & dPres,
          ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & gravDepth,
          ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dens,
          ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dDens_dPres,
          ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & mob,
          ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & dMob_dPres,
          localIndex const fluidIndex,
          integer const gravityFlag,
          real64 const dt,
          arraySlice1d<real64> const & flux,
          arraySlice2d<real64> const & fluxJacobian )
{
  localIndex constexpr numElems = FluxApproximationBase::CellStencil::NUM_POINT_IN_FLUX;
  localIndex constexpr maxStencil = FluxApproximationBase::CellStencil::MAX_STENCIL_SIZE;

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
    CellDescriptor const & cell = stencil[ke].index;

    // density
    real64 const density = dens[cell.region][cell.subRegion][fluidIndex][cell.index][0];
    real64 const dDens_dP = dDens_dPres[cell.region][cell.subRegion][fluidIndex][cell.index][0];

    // average density
    densMean        += densWeight[ke] * density;
    dDensMean_dP[ke] = densWeight[ke] * dDens_dP;
  }

  // compute potential difference MPFA-style
  real64 potDif = 0.0;
  for (localIndex ke = 0; ke < stencilSize; ++ke)
  {
    CellDescriptor const & cell = stencil[ke].index;
    localIndex const er  = cell.region;
    localIndex const esr = cell.subRegion;
    localIndex const ei  = cell.index;

    real64 weight = stencil[ke].weight;

    real64 const gravD = gravDepth[er][esr][ei];
    real64 const gravTerm = gravityFlag ? densMean * gravD : 0.0;
    real64 const dGrav_dP = gravityFlag ? dDensMean_dP[ke] * gravD : 0.0;

    potDif += weight * (pres[er][esr][ei] + dPres[er][esr][ei] - gravTerm);
    dFlux_dP[ke] = weight * (1.0 - dGrav_dP);
  }

  // upwinding of fluid properties (make this an option?)
  localIndex const k_up = (potDif >= 0) ? 0 : 1;

  CellDescriptor const & cell_up = stencil[k_up].index;
  localIndex er_up  = cell_up.region;
  localIndex esr_up = cell_up.subRegion;
  localIndex ei_up  = cell_up.index;

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
 * @brief Compute flux and its derivatives for a given connection represented by a stencil object
 * @param stencil
 * @param pres
 * @param dPres
 * @param gravDepth
 * @param dens
 * @param dDens_dPres
 * @param visc
 * @param dVisc_dPres
 * @param fluidIndex
 * @param gravityFlag
 * @param dt
 * @param flux
 * @param fluxJacobian
 *.
 * This is a specialized version for fluxes within the same region.
 * See above for a general version.
 */
inline static void
MakeFlux( localIndex const stencilSize,
          FluxApproximationBase::CellStencil::Entry const * stencil,
          arrayView1d<real64 const> const & pres,
          arrayView1d<real64 const> const & dPres,
          arrayView1d<real64 const> const & gravDepth,
          arrayView2d<real64 const> const & dens,
          arrayView2d<real64 const> const & dDens_dPres,
          arrayView1d<real64 const> const & mob,
          arrayView1d<real64 const> const & dMob_dPres,
          integer const gravityFlag,
          real64 const dt,
          arraySlice1d<real64> const & flux,
          arraySlice2d<real64> const & fluxJacobian )
{
  localIndex constexpr numElems = FluxApproximationBase::CellStencil::NUM_POINT_IN_FLUX;
  localIndex constexpr maxStencil = FluxApproximationBase::CellStencil::MAX_STENCIL_SIZE;

  stackArray1d<real64, numElems> densWeight(numElems);
  stackArray1d<real64, maxStencil> dDensMean_dP(stencilSize);
  stackArray1d<real64, maxStencil> dFlux_dP(stencilSize);

  // clear working arrays
  dDensMean_dP = 0.0;

  // density averaging weights
  densWeight = 0.5;

  // calculate quantities on primary connected cells
  real64 densMean = 0.0;
  for (localIndex i = 0; i < numElems; ++i)
  {
    CellDescriptor const & cell = stencil[i].index;

    // density
    real64 const density = dens[cell.index][0];
    real64 const dDens_dP = dDens_dPres[cell.index][0];

    // average density
    densMean += densWeight[i] * density;
    dDensMean_dP[i] = densWeight[i] * dDens_dP;
  }

  // compute potential difference MPFA-style
  real64 potDif = 0.0;
  for (localIndex ke = 0; ke < stencilSize; ++ke)
  {
    CellDescriptor const & cell = stencil[ke].index;
    localIndex const ei = cell.index;

    real64 const weight = stencil[ke].weight;

    real64 const gravD = gravDepth[ei];
    real64 const gravTerm = gravityFlag ? densMean * gravD : 0.0;
    real64 const dGrav_dP = gravityFlag ? dDensMean_dP[ke] * gravD : 0.0;

    potDif += weight * (pres[ei] + dPres[ei] - gravTerm);
    dFlux_dP[ke] = weight * (1.0 - dGrav_dP);
  }

  // upwinding of fluid properties (make this an option?)
  localIndex const k_up = (potDif >= 0) ? 0 : 1;

  CellDescriptor const & cell_up = stencil[k_up].index;
  localIndex ei_up  = cell_up.index;

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

} // namespace SinglePhaseFlowKernels

} // namespace geosx

#endif //GEOSX_SINGLEPHASEFLOWKERNELS_HPP
