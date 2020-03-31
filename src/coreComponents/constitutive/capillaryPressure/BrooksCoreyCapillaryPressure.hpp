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
 * @file BrooksCoreyCapillaryPressure.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_BROOKSCOREYCAPILLARYPRESSURE_HPP
#define GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_BROOKSCOREYCAPILLARYPRESSURE_HPP

#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"

namespace geosx
{

namespace constitutive
{

class BrooksCoreyCapillaryPressure : public CapillaryPressureBase
{
public:

  BrooksCoreyCapillaryPressure( std::string const & name,
                                dataRepository::Group * const parent );

  virtual ~BrooksCoreyCapillaryPressure() override;

  void DeliverClone( string const & name,
                     Group * const parent,
                     std::unique_ptr< ConstitutiveBase > & clone ) const override;

  static std::string CatalogName() { return "BrooksCoreyCapillaryPressure"; }

  virtual string GetCatalogName() override { return CatalogName(); }

  // CapillaryPressure-specific interface

  virtual void BatchUpdate( arrayView2d< real64 const > const & phaseVolumeFraction ) override;

  virtual void PointUpdate( arraySlice1d< real64 const > const & phaseVolFraction,
                            localIndex const k,
                            localIndex const q ) override;

  inline static void Compute( localIndex const NP,
                              arraySlice1d< real64 const > const & phaseVolFraction,
                              arraySlice1d< real64 > const & phaseCapPressure,
                              arraySlice2d< real64 > const & dPhaseCapPressure_dPhaseVolFrac,
                              arraySlice1d< integer const > const & phaseOrder,
                              arraySlice1d< real64 const > const & phaseMinVolumeFraction,
                              arraySlice1d< real64 const > const & phaseCapPressureExponentInv,
                              arraySlice1d< real64 const > const & phaseEntryPressure,
                              real64 const capPressureEpsilon,
                              real64 const volFracScale );

  struct viewKeyStruct : CapillaryPressureBase::viewKeyStruct
  {
    static constexpr auto phaseMinVolumeFractionString      = "phaseMinVolumeFraction";
    static constexpr auto phaseCapPressureExponentInvString = "phaseCapPressureExponentInv";
    static constexpr auto phaseEntryPressureString          = "phaseEntryPressure";
    static constexpr auto capPressureEpsilonString          = "capPressureEpsilon";

    using ViewKey = dataRepository::ViewKey;

    ViewKey phaseMinVolumeFraction      = { phaseMinVolumeFractionString };
    ViewKey phaseCapPressureExponentInv = { phaseCapPressureExponentInvString };
    ViewKey phaseEntryPressure          = { phaseEntryPressureString };
    ViewKey capPressureEpsilon          = { capPressureEpsilonString };

  } viewKeysBrooksCoreyCapillaryPressure;

protected:
  virtual void PostProcessInput() override;

  inline static void EvaluateBrooksCoreyFunction( real64 const scaledWettingVolFrac,
                                                  real64 const dScaledWettingPhaseVolFrac_dVolFrac,
                                                  real64 & phaseCapPressure,
                                                  real64 & dPhaseCapPressure_dVolFrac,
                                                  real64 const exponentInv,
                                                  real64 const entryPressure,
                                                  real64 const eps );

  array1d< real64 > m_phaseMinVolumeFraction;
  array1d< real64 > m_phaseCapPressureExponentInv;
  array1d< real64 > m_phaseEntryPressure;

  real64 m_capPressureEpsilon;
  real64 m_volFracScale;
};


inline void
BrooksCoreyCapillaryPressure::Compute( localIndex const NP,
                                       arraySlice1d< real64 const > const & phaseVolFraction,
                                       arraySlice1d< real64 > const & phaseCapPressure,
                                       arraySlice2d< real64 > const & dPhaseCapPressure_dVolFrac,
                                       arraySlice1d< integer const > const & phaseOrder,
                                       arraySlice1d< real64 const > const & phaseMinVolumeFraction,
                                       arraySlice1d< real64 const > const & phaseCapPressureExponentInv,
                                       arraySlice1d< real64 const > const & phaseEntryPressure,
                                       real64 const capPressureEpsilon,
                                       real64 const volFracScale )
{

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    for( localIndex jp = 0; jp < NP; ++jp )
    {
      dPhaseCapPressure_dVolFrac[ip][jp] = 0.0;
    }
  }

  real64 const volFracScaleInv = 1.0 / volFracScale;

  // the Brooks-Corey model does not support volFracScaled = 0,
  // hence we need an epsilon value to avoid a division by zero
  // TODO: for S < epsilon, replace the original unbounded BC curve with a bounded power-law extension
  real64 const eps = capPressureEpsilon;


  // compute first water-oil capillary pressure as a function of water-phase vol fraction
  integer const ip_water = phaseOrder[CapillaryPressureBase::PhaseType::WATER];
  if( ip_water >= 0 )
  {
    real64 const volFracScaled = (phaseVolFraction[ip_water] - phaseMinVolumeFraction[ip_water]) * volFracScaleInv;
    real64 const exponentInv   = phaseCapPressureExponentInv[ip_water];
    real64 const entryPressure = phaseEntryPressure[ip_water];

    real64 const wettingVolFracScaled           = volFracScaled;
    real64 const dWettingVolFracScaled_dVolFrac = volFracScaleInv;

    EvaluateBrooksCoreyFunction( wettingVolFracScaled,
                                 dWettingVolFracScaled_dVolFrac,
                                 phaseCapPressure[ip_water],
                                 dPhaseCapPressure_dVolFrac[ip_water][ip_water],
                                 exponentInv,
                                 entryPressure,
                                 eps );

  }


  // compute first gas-oil capillary pressure as a function of gas-phase vol fraction
  integer const ip_gas = phaseOrder[CapillaryPressureBase::PhaseType::GAS];
  if( ip_gas >= 0 )
  {
    real64 const volFracScaled = (phaseVolFraction[ip_gas] - phaseMinVolumeFraction[ip_gas]) * volFracScaleInv;
    real64 const exponentInv   = phaseCapPressureExponentInv[ip_gas];
    real64 const entryPressure = -phaseEntryPressure[ip_gas]; // for gas capillary pressure, take the opposite of the VG
                                                              // function

    real64 const wettingVolFracScaled           = 1-volFracScaled;
    real64 const dWettingVolFracScaled_dVolFrac =  -volFracScaleInv;

    EvaluateBrooksCoreyFunction( wettingVolFracScaled,
                                 dWettingVolFracScaled_dVolFrac,
                                 phaseCapPressure[ip_gas],
                                 dPhaseCapPressure_dVolFrac[ip_gas][ip_gas],
                                 exponentInv,
                                 entryPressure,
                                 eps );

  }
}


inline void
BrooksCoreyCapillaryPressure::EvaluateBrooksCoreyFunction( real64 const scaledWettingVolFrac,
                                                           real64 const dScaledWettingPhaseVolFrac_dVolFrac,
                                                           real64 & phaseCapPressure,
                                                           real64 & dPhaseCapPressure_dVolFrac,
                                                           real64 const exponentInv,
                                                           real64 const entryPressure,
                                                           real64 const eps )
{
  real64 const exponent = 1 / exponentInv; // div by 0 taken care of by initialization check

  phaseCapPressure           = 0.0;
  dPhaseCapPressure_dVolFrac = 0.0;

  if( scaledWettingVolFrac >= eps && scaledWettingVolFrac < 1.0 )
  {
    // intermediate value
    real64 const val = entryPressure / std::pow( scaledWettingVolFrac, exponent + 1 );

    phaseCapPressure           = val * scaledWettingVolFrac; // entryPressure * (S_w)^( - 1 / exponentInv )
    dPhaseCapPressure_dVolFrac = -dScaledWettingPhaseVolFrac_dVolFrac * val * exponent;
  }
  else // enforce a constant and bounded capillary pressure
  {
    phaseCapPressure = (scaledWettingVolFrac < eps)
                     ? entryPressure / std::pow( eps, exponent ) // div by 0 taken care of by initialization check
                     : entryPressure;
  }

}


} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_BROOKSCOREYCAPILLARYPRESSURE_HPP
