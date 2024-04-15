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
 * @file WellElementKernellUtilities.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLELEMENTKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLELEMENTKERNELS_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/capillaryPressure/layouts.hpp"
 


namespace geos
{

namespace wellElementKernelUtilities
{

// TODO make input parameter
static constexpr real64 epsC1PPU = 5000;


using Deriv = constitutive::multifluid::DerivativeOffset;

struct PotGrad
{
  template< integer NC, integer ISTHERMAL >
  GEOS_HOST_DEVICE
  static void
  compute (real64 const & gravCoef,
           real64 const & gravCoefNext,
           real64 const & pres,
           real64 const & presNext,
           real64 const & totalMassDens,
           real64 const & totalMassDensNext,
           arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dTotalMassDens,
           arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dTotalMassDensNext,
           arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dTotalMassDens_dCompDens,
           arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dTotalMassDens_dCompDensNext,
           real64 & potDiff,
           real64 ( & dPotDiff)[2][NC+1+IS_THERMAL] )
 
  {
      // local working variables and arrays
  real64 pres[2]{};
  real64 dPres_dP[2]{};
  real64 dPres_dC[2][NC]{};
  real64 dFlux_dP[2]{};
  real64 dFlux_dC[2][NC]{};
  real64 dMult_dP[2]{};
  real64 dMult_dC[2][NC]{};
  real64 dPotDiff_dP[2]{};
  real64 dPotDiff_dC[2][NC]{};

  real64 multiplier[2]{};
  multiplier[ElemTag::CURRENT] = 1.0;
  multiplier[ElemTag::NEXT]    = -1.0;
  
  pres[ElemTag::CURRENT] = pres;
  pres[ElemTag::NEXT]    = presNext;


    // local working variables and arrays
  real64 dAvgMassDens_dCompCurrent[NC]{};
  real64 dAvgMassDens_dCompNext[NC]{};

  // compute the average density at the interface between well elements
  real64 const avgMassDens = 0.5 * ( totalMassDensNext + totalMassDens );

  real64 const gravD = gravCoefNext - gravCoef;
  pres[0] += 

  real64 const dAvgMassDens_dPresNext    = 0.5 * dTotalMassDensNext[Deriv::dP];
  real64 const dAvgMassDens_dPresCurrent = 0.5 * dTotalMassDens[Deriv::dP];
  for( integer ic = 0; ic < NC; ++ic )
  {
    dAvgMassDens_dCompNext[ic]    = 0.5 * dTotalMassDensNext[Deriv::dC+ic];
    dAvgMassDens_dCompCurrent[ic] = 0.5 * dTotalMassDens[Deriv::dC+ic];
  }

  // compute depth diff times acceleration
  real64 const gravD = gravCoefNext - gravCoef;

  potDiff  = presNext - pres  - avgMassDens * gravD;
  
  // TODO: add friction and acceleration terms
 
  /*
  localPresRel = ( presNext - pres - avgMassDens * gravD );
  dpot_dp_next = ( 1 - dAvgMassDens_dPresNext * gravD );
  dpot_dp_cur  = ( -1 - dAvgMassDens_dPresCurrent * gravD );
  for( integer ic = 0; ic < NC; ++ic )
  {
    localPresRelJacobian[TAG::NEXT *(NC+1+IS_THERMAL) + ic+1]    = -dAvgMassDens_dCompNext[ic] * gravD;
    localPresRelJacobian[TAG::CURRENT *(NC+1) + ic+1] = -dAvgMassDens_dCompCurrent[ic] * gravD;
  }
  if constexpr ( IS_THERMAL )
  {
    localPresRelJacobian[TAG::NEXT *(NC+1+IS_THERMAL)+NC+1]    =  0.5 * dTotalMassDensNext[Deriv::dT];
    localPresRelJacobian[TAG::CURRENT *(NC+1)+1] = 0.5 * dTotalMassDens[Deriv::dT];
  }
  */
  }

};



struct WellElementPhaseFlux
{
  /**
   * @brief Form the PhasePotentialUpwind from pressure gradient and gravitational head
   * @tparam numComp number of components
   * @param numPhase number of phases
   * @param ip phase index
   * @param pres pressure
   * @param gravCoef gravitational coefficient
   * @param dPhaseVolFrac derivative of phase volume fraction wrt pressure, temperature, comp density
   * @param dCompFrac_dCompDens derivative of component fraction wrt component density
   * @param phaseMassDens phase mass density
   * @param dPhaseMassDens derivative of phase mass density wrt pressure, temperature, comp fraction
   * @param potGrad potential gradient for this phase
   * @param phaseFlux phase flux
   * @param dPhaseFlux  derivatives of phase flux  
   */
  template< integer numComp, integer IS_THERMAL >
  GEOS_HOST_DEVICE
  static void
  compute( integer const ip,
           real64 const & gravCoef,
           real64 const & gravCoefNext,
           real64 const & pres,
           real64 const & presNext,
           real64 const & totalMassDens,
           real64 const & totalMassDensNext,
           arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dTotalMassDens,
           arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dTotalMassDensNext,
           arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dTotalMassDens_dCompDens,
           arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dTotalMassDens_dCompDensNext,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFrac,
           arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dPhaseVolFrac,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseMassDens,
           arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dPhaseMassDens,
           real64 & potGrad,
           real64 ( &phaseFlux ),
           real64 ( & dPhaseFlux )[numComp+1+IS_THERMAL] )
  {
    real64 dPresGrad_dP[numFluxSupportPoints]{};
    real64 dPresGrad_dC[numFluxSupportPoints][numComp]{};
    real64 dGravHead_dP[numFluxSupportPoints]{};
    real64 dGravHead_dC[numFluxSupportPoints][numComp]{};
    PotGrad::compute< numComp, IS_THERMAL >( gravCoef 
                                             ,gravCoefNext
                                             , pres
                                             ,presNext
                                             ,totalMassDens
                                             , totalMassDensNext
                                             , dTotalMassDens
                                             , dTotalMassDensNext
                                             , dTotalMassDens_dCompDens
                                             , dTotalMassDens_dCompDensNex );

    // gravity head
    real64 gravHead = gravCoef - gravCoefNext
    for( integer i = 0; i < numFluxSupportPoints; i++ )
    {
      localIndex const er  = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei  = sei[i];

      real64 const gravD = trans[i] * gravCoef[er][esr][ei];

      gravHead += gravD;
    }

    // *** upwinding ***

    // phase flux and derivatives

    // assuming TPFA in the code below

    real64 Ttrans = fabs( trans[0] );
    potGrad = potGrad / Ttrans;

    real64 const mobility_i = phaseMob[seri[0]][sesri[0]][sei[0]][ip];
    real64 const mobility_j = phaseMob[seri[1]][sesri[1]][sei[1]][ip];

    // compute phase flux, see Eqs. (66) and (69) from the reference above
    real64 smoEps = epsC1PPU;
    if( fabs( gravHead ) <= 1e-20 )
      smoEps = 1000;
    real64 const tmpSqrt = sqrt( potGrad * potGrad + smoEps * smoEps );
    real64 const smoMax = 0.5 * (-potGrad + tmpSqrt);

    phaseFlux = potGrad * mobility_i - smoMax * (mobility_j - mobility_i) );

    // derivativess

    // first part, mobility derivative

    // dP
    {
      real64 const dMob_dP = dPhaseMob[seri[0]][sesri[0]][sei[0]][ip][Deriv::dP];
      dPhaseFlux_dP[0] += Ttrans * potGrad * dMob_dP;
    }

    // dC
    {
      arraySlice1d< real64 const, compflow::USD_PHASE_DC - 2 >
      dPhaseMobSub = dPhaseMob[seri[0]][sesri[0]][sei[0]][ip];
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseFlux_dC[0][jc] += Ttrans * potGrad * dPhaseMobSub[Deriv::dC + jc];
      }
    }

    real64 const tmpInv = 1.0 / tmpSqrt;
    real64 const dSmoMax_x = 0.5 * (1.0 - potGrad * tmpInv);

    // pressure gradient and mobility difference depend on all points in the stencil
    real64 const dMobDiff_sign[numFluxSupportPoints] = {-1.0, 1.0};
    for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      // dP

      real64 const dPotGrad_dP = dPresGrad_dP[ke] - dGravHead_dP[ke];

      // first part
      dPhaseFlux_dP[ke] += dPotGrad_dP * mobility_i;

      // second part
      real64 const dSmoMax_dP = -dPotGrad_dP * dSmoMax_x;
      dPhaseFlux_dP[ke] += -dSmoMax_dP * (mobility_j - mobility_i);

      real64 const dMob_dP = dPhaseMob[seri[ke]][sesri[ke]][sei[ke]][ip][Deriv::dP];
      dPhaseFlux_dP[ke] += -Ttrans * smoMax * dMobDiff_sign[ke] * dMob_dP;

      // dC

      arraySlice1d< real64 const, compflow::USD_PHASE_DC - 2 >
      dPhaseMobSub = dPhaseMob[seri[ke]][sesri[ke]][sei[ke]][ip];

      for( integer jc = 0; jc < numComp; ++jc )
      {
        real64 const dPotGrad_dC = dPresGrad_dC[ke][jc] - dGravHead_dC[ke][jc];

        // first part
        dPhaseFlux_dC[ke][jc] += dPotGrad_dC * mobility_i;

        // second part
        real64 const dSmoMax_dC = -dPotGrad_dC * dSmoMax_x;
        dPhaseFlux_dC[ke][jc] += -dSmoMax_dC * (mobility_j - mobility_i);
        dPhaseFlux_dC[ke][jc] += -Ttrans * smoMax * dMobDiff_sign[ke] * dPhaseMobSub[Deriv::dC + jc];
      }
    }


  }
};


} // namespace wellElementKernelUtilities

} // namespace geosx


#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLELEMENTKERNELS_HPP
