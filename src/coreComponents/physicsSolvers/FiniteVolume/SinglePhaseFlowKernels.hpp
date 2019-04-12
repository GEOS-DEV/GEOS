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
 * @file SinglePhaseFlow_kernels.hpp
 */

#ifndef GEOSX_SINGLEPHASEFLOWKERNELS_HPP
#define GEOSX_SINGLEPHASEFLOWKERNELS_HPP

#include "common/DataTypes.hpp"

namespace geosx
{


namespace SinglePhaseFlowKernels
{

template<bool ISPORO>
struct AssembleAccumulationTermsHelper;

template<>
struct AssembleAccumulationTermsHelper<true>
{
  inline static constexpr void porosityUpdate( real64 & poro,
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
  inline static constexpr void porosityUpdate( real64 & poro,
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

template<bool COUPLED>
inline static void MakeAccumulation( real64 const & dPres,
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
  localAccumJacobian = dPoro_dPres * densNew * volNew + dDens_dPres * poroNew * volNew;
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
 *
 * This is a general version that assumes different element regions.
 * See below for a specialized version for fluxes within a region.
 */
inline static void MakeFlux( StencilCollection<CellDescriptor, real64>::Accessor stencil,
                             ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & pres,
                             ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & dPres,
                             ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & gravDepth,
                             ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dens,
                             ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dDens_dPres,
                             ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & visc,
                             ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dVisc_dPres,
                             localIndex const fluidIndex,
                             integer const gravityFlag,
                             real64 const dt,
                             arraySlice1d<real64> const & flux,
                             arraySlice2d<real64> const & fluxJacobian )
{
  localIndex constexpr numElems = StencilCollection<CellDescriptor, real64>::NUM_POINT_IN_FLUX;
  localIndex constexpr maxStencil = StencilCollection<CellDescriptor, real64>::MAX_STENCIL_SIZE;
  localIndex const stencilSize = stencil.size();

  stackArray1d<real64, numElems> densWeight(numElems);
  stackArray1d<real64, numElems> mobility(numElems);
  stackArray1d<real64, numElems> dMobility_dP(numElems);
  stackArray1d<real64, maxStencil> dDensMean_dP(stencilSize);
  stackArray1d<real64, maxStencil> dFlux_dP(stencilSize);

  // clear working arrays
  dDensMean_dP = 0.0;

  // density averaging weights
  densWeight = 0.5;

  // calculate quantities on primary connected cells
  real64 densMean = 0.0;
  stencil.forConnected([&](CellDescriptor const & cell, localIndex i)
  {
    localIndex const er = cell.region;
    localIndex const esr = cell.subRegion;
    localIndex const ei = cell.index;

    // density
    real64 const density = dens[er][esr][fluidIndex][ei][0];
    real64 const dDens_dP = dDens_dPres[er][esr][fluidIndex][ei][0];

    // viscosity
    real64 const viscosity = visc[er][esr][fluidIndex][ei][0];
    real64 const dVisc_dP = dVisc_dPres[er][esr][fluidIndex][ei][0];

    // mobility
    mobility[i] = density / viscosity;
    dMobility_dP[i] = dDens_dP / viscosity - mobility[i] / viscosity * dVisc_dP;

    // average density
    densMean += densWeight[i] * density;
    dDensMean_dP[i] = densWeight[i] * dDens_dP;
  });

  // compute potential difference MPFA-style
  real64 potDif = 0.0;
  stencil.forAll([&](CellDescriptor const & cell, real64 w, localIndex i)
  {
    localIndex const er = cell.region;
    localIndex const esr = cell.subRegion;
    localIndex const ei = cell.index;

    real64 const gravD = gravDepth[er][esr][ei];
    real64 const gravTerm = gravityFlag ? densMean * gravD : 0.0;
    real64 const dGrav_dP = gravityFlag ? dDensMean_dP[i] * gravD : 0.0;

    potDif += w * (pres[er][esr][ei] + dPres[er][esr][ei] - gravTerm);
    dFlux_dP[i] = w * (1.0 - dGrav_dP);
  });

  // upwinding of fluid properties (make this an option?)
  localIndex const k_up = (potDif >= 0) ? 0 : 1;

  // compute the final flux and derivatives
  real64 const fluxVal = mobility[k_up] * potDif;
  for (localIndex ke = 0; ke < stencilSize; ++ke)
  {
    dFlux_dP[ke] *= mobility[k_up];
  }

  dFlux_dP[k_up] += dMobility_dP[k_up] * potDif;

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
inline static void MakeFlux( StencilCollection<localIndex, real64>::Accessor stencil,
                             arrayView1d<real64 const> const & pres,
                             arrayView1d<real64 const> const & dPres,
                             arrayView1d<real64 const> const & gravDepth,
                             arrayView2d<real64 const> const & dens,
                             arrayView2d<real64 const> const & dDens_dPres,
                             arrayView2d<real64 const> const & visc,
                             arrayView2d<real64 const> const & dVisc_dPres,
                             integer const gravityFlag,
                             real64 const dt,
                             arraySlice1d<real64> const & flux,
                             arraySlice2d<real64> const & fluxJacobian )
{
  localIndex constexpr numElems = StencilCollection<CellDescriptor, real64>::NUM_POINT_IN_FLUX;
  localIndex constexpr maxStencil = StencilCollection<CellDescriptor, real64>::MAX_STENCIL_SIZE;
  localIndex const stencilSize = stencil.size();

  stackArray1d<real64, numElems> densWeight(numElems);
  stackArray1d<real64, numElems> mobility(numElems);
  stackArray1d<real64, numElems> dMobility_dP(numElems);
  stackArray1d<real64, maxStencil> dDensMean_dP(stencilSize);
  stackArray1d<real64, maxStencil> dFlux_dP(stencilSize);

  // clear working arrays
  dDensMean_dP = 0.0;

  // density averaging weights
  densWeight = 0.5;

  // calculate quantities on primary connected cells
  real64 densMean = 0.0;
  stencil.forConnected([&](localIndex const ei, localIndex const i)
  {
    // density
    real64 const density = dens[ei][0];
    real64 const dDens_dP = dDens_dPres[ei][0];

    // viscosity
    real64 const viscosity = visc[ei][0];
    real64 const dVisc_dP = dVisc_dPres[ei][0];

    // mobility
    mobility[i] = density / viscosity;
    dMobility_dP[i] = dDens_dP / viscosity - mobility[i] / viscosity * dVisc_dP;

    // average density
    densMean += densWeight[i] * density;
    dDensMean_dP[i] = densWeight[i] * dDens_dP;
  });

  // compute potential difference MPFA-style
  real64 potDif = 0.0;
  stencil.forAll([&](localIndex const ei, real64 const w, localIndex const i)
  {
    real64 const gravD = gravDepth[ei];
    real64 const gravTerm = gravityFlag ? densMean * gravD : 0.0;
    real64 const dGrav_dP = gravityFlag ? dDensMean_dP[i] * gravD : 0.0;

    potDif += w * (pres[ei] + dPres[ei] - gravTerm);
    dFlux_dP[i] = w * (1.0 - dGrav_dP);
  });

  // upwinding of fluid properties (make this an option?)
  localIndex const k_up = (potDif >= 0) ? 0 : 1;

  // compute the final flux and derivatives
  real64 const fluxVal = mobility[k_up] * potDif;
  for (localIndex ke = 0; ke < stencilSize; ++ke)
  {
    dFlux_dP[ke] *= mobility[k_up];
  }

  dFlux_dP[k_up] += dMobility_dP[k_up] * potDif;

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
