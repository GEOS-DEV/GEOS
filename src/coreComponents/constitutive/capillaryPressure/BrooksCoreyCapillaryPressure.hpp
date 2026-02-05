/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file BrooksCoreyCapillaryPressure.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_BROOKSCOREYCAPILLARYPRESSURE_HPP
#define GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_BROOKSCOREYCAPILLARYPRESSURE_HPP

#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"

namespace geos
{

namespace constitutive
{

class BrooksCoreyCapillaryPressureUpdate final : public CapillaryPressureBaseUpdate
{
public:

  BrooksCoreyCapillaryPressureUpdate( arrayView1d< real64 const > const & phaseMinVolumeFraction,
                                      arrayView1d< real64 const > const & phaseCapPressureExponentInv,
                                      arrayView1d< real64 const > const & phaseEntryPressure,
                                      real64 const capPressureEpsilon,
                                      real64 const volFracScale,
                                      arrayView1d< integer const > const & phaseTypes,
                                      arrayView1d< integer const > const & phaseOrder,
                                      arrayView3d< real64, cappres::USD_CAPPRES > const & phaseCapPressure,
                                      arrayView4d< real64, cappres::USD_CAPPRES_DS > const & dPhaseCapPressure_dPhaseVolFrac )
    : CapillaryPressureBaseUpdate( phaseTypes,
                                   phaseOrder,
                                   phaseCapPressure,
                                   dPhaseCapPressure_dPhaseVolFrac ),
    m_phaseMinVolumeFraction( phaseMinVolumeFraction ),
    m_phaseCapPressureExponentInv( phaseCapPressureExponentInv ),
    m_phaseEntryPressure( phaseEntryPressure ),
    m_capPressureEpsilon( capPressureEpsilon ),
    m_volFracScale( volFracScale )
  {}

  GEOS_HOST_DEVICE
  void compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                arraySlice1d< real64, cappres::USD_CAPPRES - 2 > const & phaseCapPres,
                arraySlice2d< real64, cappres::USD_CAPPRES_DS - 2 > const & dPhaseCapPres_dPhaseVolFrac ) const;

  GEOS_HOST_DEVICE
  void computeInv( arraySlice1d< real64, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                   arraySlice1d< real64 const, cappres::USD_CAPPRES - 2 > const & phaseCapPres,
                   arraySlice2d< real64, cappres::USD_CAPPRES_DS - 2 > const & dPhaseCapPres_dPhaseVolFrac ) const;

  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const override
  {
    compute( phaseVolFraction,
             m_phaseCapPressure[k][q],
             m_dPhaseCapPressure_dPhaseVolFrac[k][q] );
  }

private:

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  evaluateBrooksCoreyFunction( real64 const scaledWettingVolFrac,
                               real64 const dScaledWettingPhaseVolFrac_dVolFrac,
                               real64 const exponentInv,
                               real64 const entryPressure,
                               real64 const eps,
                               real64 & phaseCapPressure,
                               real64 & dPhaseCapPressure_dVolFrac );


  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  evaluateBrooksCoreyFunctionInv( real64 const phaseCapPressure,
                                  int const ip,
                                  real64 const volFracScaleInv,
                                  real64 const exponentInv,
                                  real64 const entryPressure,
                                  real64 const maxCapPres_eps,
                                  real64 const phaseMinVolumeFraction,
                                  arrayView1d< integer const > const phaseOrder,
                                  real64 & phaseVolFraction,
                                  real64 & dPhaseCapPressure_dVolFrac );

  arrayView1d< real64 const > m_phaseMinVolumeFraction;
  arrayView1d< real64 const > m_phaseCapPressureExponentInv;
  arrayView1d< real64 const > m_phaseEntryPressure;

  real64 m_capPressureEpsilon;
  real64 m_volFracScale;

};



class BrooksCoreyCapillaryPressure : public CapillaryPressureBase
{
public:

  BrooksCoreyCapillaryPressure( string const & name,
                                dataRepository::Group * const parent );

  static string catalogName() { return "BrooksCoreyCapillaryPressure"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = BrooksCoreyCapillaryPressureUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  struct viewKeyStruct : CapillaryPressureBase::viewKeyStruct
  {
    static constexpr char const * phaseCapPressureExponentInvString() { return "phaseCapPressureExponentInv"; }
    static constexpr char const * phaseEntryPressureString() { return "phaseEntryPressure"; }
    static constexpr char const * capPressureEpsilonString() { return "capPressureEpsilon"; }
    static constexpr char const * volFracScaleString() { return "volFracScale"; }
  };

protected:

  virtual void postInputInitialization() override;

  array1d< real64 > m_phaseCapPressureExponentInv;
  array1d< real64 > m_phaseEntryPressure;

  real64 m_capPressureEpsilon;
  real64 m_volFracScale;
};


GEOS_HOST_DEVICE
inline void
BrooksCoreyCapillaryPressureUpdate::
  compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
           arraySlice1d< real64, cappres::USD_CAPPRES - 2 > const & phaseCapPres,
           arraySlice2d< real64, cappres::USD_CAPPRES_DS - 2 > const & dPhaseCapPres_dPhaseVolFrac ) const
{
  LvArray::forValuesInSlice( dPhaseCapPres_dPhaseVolFrac, []( real64 & val ){ val = 0.0; } );

  real64 const volFracScaleInv = 1.0 / m_volFracScale;

  // the Brooks-Corey model does not support volFracScaled = 0,
  // hence we need an epsilon value to avoid a division by zero
  // TODO: for S < epsilon, replace the original unbounded BC curve with a bounded power-law extension
  real64 const eps = m_capPressureEpsilon;


  // compute first water-oil capillary pressure as a function of water-phase vol fraction
  integer const ip_water = m_phaseOrder[CapillaryPressureBase::PhaseType::WATER];
  if( ip_water >= 0 )
  {
    real64 const volFracScaled = (phaseVolFraction[ip_water] - m_phaseMinVolumeFraction[ip_water]) * volFracScaleInv;
    real64 const exponentInv   = m_phaseCapPressureExponentInv[ip_water];
    real64 const entryPressure = m_phaseEntryPressure[ip_water];

    real64 const wettingVolFracScaled           = volFracScaled;
    real64 const dWettingVolFracScaled_dVolFrac = volFracScaleInv;

    evaluateBrooksCoreyFunction( wettingVolFracScaled,
                                 dWettingVolFracScaled_dVolFrac,
                                 exponentInv,
                                 entryPressure,
                                 eps,
                                 phaseCapPres[ip_water],
                                 dPhaseCapPres_dPhaseVolFrac[ip_water][ip_water] );

  }


  // compute first gas-oil capillary pressure as a function of gas-phase vol fraction
  integer const ip_gas = m_phaseOrder[CapillaryPressureBase::PhaseType::GAS];
  integer const ip_oil = m_phaseOrder[CapillaryPressureBase::PhaseType::OIL];

  GEOS_UNUSED_VAR( ip_gas );
  if( ip_oil >= 0 )
  {
    real64 const volFracScaled = (phaseVolFraction[ip_oil] - m_phaseMinVolumeFraction[ip_oil]) * volFracScaleInv;
    real64 const exponentInv   = m_phaseCapPressureExponentInv[ip_oil];
    real64 const entryPressure = -m_phaseEntryPressure[ip_oil]; // for gas capillary pressure, take the opposite of the
                                                                // BC function

    real64 const wettingVolFracScaled           = 1-volFracScaled;
    real64 const dWettingVolFracScaled_dVolFrac =  -volFracScaleInv;

    evaluateBrooksCoreyFunction( wettingVolFracScaled,
                                 dWettingVolFracScaled_dVolFrac,
                                 exponentInv,
                                 entryPressure,
                                 eps,
                                 phaseCapPres[ip_oil],
                                 dPhaseCapPres_dPhaseVolFrac[ip_oil][ip_oil] );
  }
}

GEOS_HOST_DEVICE
inline void
BrooksCoreyCapillaryPressureUpdate::
  computeInv( arraySlice1d< real64, compflow::USD_PHASE - 1 > const & phaseVolFraction,
              arraySlice1d< real64 const, cappres::USD_CAPPRES - 2 > const & phaseCapPres,
              arraySlice2d< real64, cappres::USD_CAPPRES_DS - 2 > const & dPhaseCapPres_dPhaseVolFrac ) const
{
  LvArray::forValuesInSlice( dPhaseCapPres_dPhaseVolFrac, []( real64 & val ){ val = 0.0; } );

  real64 const volFracScaleInv = 1.0 / m_volFracScale;

  // the Brooks-Corey model does not support volFracScaled = 0,
  // hence we need an epsilon value to avoid a division by zero
  // TODO: for S < epsilon, replace the original unbounded BC curve with a bounded power-law extension
  real64 const eps = m_capPressureEpsilon;


  // compute first water-oil capillary pressure as a function of water-phase vol fraction
  integer const ip_water = m_phaseOrder[CapillaryPressureBase::PhaseType::WATER];
  integer const ip_gas = m_phaseOrder[CapillaryPressureBase::PhaseType::GAS];
  if( ip_water >= 0 )
  {
    real64 const volFracScaled_eps = (eps - m_phaseMinVolumeFraction[ip_water]) * volFracScaleInv;
    real64 const exponentInv   = m_phaseCapPressureExponentInv[ip_water];
    real64 const entryPressure = m_phaseEntryPressure[ip_water];

    real64 const wettingVolFracScaled_eps           = volFracScaled_eps;
    real64 const dWettingVolFracScaled_dVolFrac = volFracScaleInv;

    real64 maxCapPres_eps = 0.0;
    real64 max_dpc_eps = 0.0;

    evaluateBrooksCoreyFunction( wettingVolFracScaled_eps,
                                 dWettingVolFracScaled_dVolFrac,
                                 exponentInv,
                                 entryPressure,
                                 eps,
                                 maxCapPres_eps,
                                 max_dpc_eps );

    evaluateBrooksCoreyFunctionInv( phaseCapPres[ip_water],
                                    ip_water,
                                    volFracScaleInv,
                                    exponentInv,
                                    entryPressure,
                                    maxCapPres_eps,
                                    m_phaseMinVolumeFraction[ip_water],
                                    m_phaseOrder,
                                    phaseVolFraction[ip_water],
                                    dPhaseCapPres_dPhaseVolFrac[ip_water][ip_water] );
    phaseVolFraction[ip_gas] = 1.0 - phaseVolFraction[ip_water];
  }


  // compute first gas-oil capillary pressure as a function of gas-phase vol fraction


  // if( ip_gas >= 0 )
  // {
  //   real64 const volFracScaled_eps = (eps - m_phaseMinVolumeFraction[ip_gas]) * volFracScaleInv;
  //   real64 const exponentInv   = m_phaseCapPressureExponentInv[ip_gas];
  //   real64 const entryPressure = -m_phaseEntryPressure[ip_gas]; // for gas capillary pressure, take the opposite of the
  //                                                               // BC function

  //   real64 const wettingVolFracScaled_eps           = 1-volFracScaled_eps;
  //   real64 const dWettingVolFracScaled_dVolFrac =  -volFracScaleInv;

  //   real64 maxCapPres_eps = 0.0;
  //   real64 max_dpc_eps = 0.0;

  //   evaluateBrooksCoreyFunction( wettingVolFracScaled_eps,
  //                                dWettingVolFracScaled_dVolFrac,
  //                                exponentInv,
  //                                entryPressure,
  //                                eps,
  //                                maxCapPres_eps,
  //                                max_dpc_eps );

  //   evaluateBrooksCoreyFunctionInv( phaseCapPres[ip_gas],
  //                                   ip_gas,
  //                                   volFracScaleInv,
  //                                   exponentInv,
  //                                   entryPressure,
  //                                   maxCapPres_eps,
  //                                   m_phaseMinVolumeFraction[ip_gas],
  //                                   m_phaseOrder,
  //                                   phaseVolFraction[ip_gas],
  //                                   dPhaseCapPres_dPhaseVolFrac[ip_gas][ip_gas] );
  //  phaseVolFraction[ip_water] = 1.0 - phaseVolFraction[ip_gas];
  // }
}

GEOS_HOST_DEVICE
inline void
BrooksCoreyCapillaryPressureUpdate::
  evaluateBrooksCoreyFunction( real64 const scaledWettingVolFrac,
                               real64 const dScaledWettingPhaseVolFrac_dVolFrac,
                               real64 const exponentInv,
                               real64 const entryPressure,
                               real64 const eps,
                               real64 & phaseCapPressure,
                               real64 & dPhaseCapPressure_dVolFrac )
{
  real64 const exponent = 1.0 / exponentInv; // div by 0 taken care of by initialization check

  phaseCapPressure           = 0.0;
  dPhaseCapPressure_dVolFrac = 0.0;

  if( scaledWettingVolFrac >= eps && scaledWettingVolFrac < 1.0 )
  {
    // intermediate value
    real64 const val = entryPressure / pow( scaledWettingVolFrac, exponent + 1 );

    phaseCapPressure           = val * scaledWettingVolFrac; // entryPressure * (S_w)^( - 1 / exponentInv )
    dPhaseCapPressure_dVolFrac = -dScaledWettingPhaseVolFrac_dVolFrac * val * exponent;
  }
  else // enforce a constant and bounded capillary pressure
  {
    phaseCapPressure = (scaledWettingVolFrac < eps)
                     ? entryPressure / pow( eps, exponent ) // div by 0 taken care of by initialization check
                     : entryPressure;
  }

}

GEOS_HOST_DEVICE
inline void
BrooksCoreyCapillaryPressureUpdate::
  evaluateBrooksCoreyFunctionInv( real64 const phaseCapPressure,
                                  int const ip,
                                  real64 const volFracScaleInv,
                                  real64 const exponentInv,
                                  real64 const entryPressure,
                                  real64 const maxCapPres_eps,
                                  real64 const phaseMinVolumeFraction,
                                  arrayView1d< integer const > const phaseOrder,
                                  real64 & phaseVolFraction,
                                  real64 & dPhaseCapPressure_dVolFrac )
{


  phaseVolFraction           = 0.0;
  real64 value               = 0.0;
  dPhaseCapPressure_dVolFrac = 0.0;
  integer const ip_oil = phaseOrder[CapillaryPressureBase::PhaseType::OIL];

  real64 const dScaledWettingPhaseVolFrac_dVolFrac = (ip == ip_oil)
  ? -volFracScaleInv : volFracScaleInv;

  if( phaseCapPressure <= maxCapPres_eps && phaseCapPressure >= entryPressure )
  {
    // intermediate value
    real64 const val =  pow( entryPressure, exponentInv ) / pow( phaseCapPressure, exponentInv + 1 );

    value           = (phaseCapPressure  * val) * volFracScaleInv + phaseMinVolumeFraction; // entryPressure * (S_w)^( - 1 / exponentInv )
    dPhaseCapPressure_dVolFrac = -dScaledWettingPhaseVolFrac_dVolFrac * val * exponentInv;
    phaseVolFraction = (ip == ip_oil) ? 1.0 - value : value;
  }
  else // enforce a constant and bounded capillary pressure
  {
    real64 const val = (phaseCapPressure > maxCapPres_eps)
                     ? pow( entryPressure, exponentInv ) / pow( maxCapPres_eps, exponentInv ) : 1.0;
    value           = val * volFracScaleInv + phaseMinVolumeFraction;
    phaseVolFraction = (ip == ip_oil) ? 1.0 - value : value;
  }

}


} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_BROOKSCOREYCAPILLARYPRESSURE_HPP
