/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
  * @file VanGenuchtenCapillaryPressure.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_VANGENUCHTENCAPILLARYPRESSURE_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_VANGENUCHTENCAPILLARYPRESSURE_HPP

#include "constitutive/CapillaryPressure/CapillaryPressureBase.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const vanGenuchtenCapillaryPressure = "VanGenuchtenCapillaryPressure";
}
}

namespace constitutive
{

class VanGenuchtenCapillaryPressure : public CapillaryPressureBase
{
public:

  VanGenuchtenCapillaryPressure( std::string const & name,
                                 dataRepository::ManagedGroup * const parent );

  virtual ~VanGenuchtenCapillaryPressure() override;

  std::unique_ptr<ConstitutiveBase> DeliverClone( string const & name,
                                                  ManagedGroup * const parent ) const override;

  static std::string CatalogName() { return dataRepository::keys::vanGenuchtenCapillaryPressure; }

  virtual string GetCatalogName() override { return CatalogName(); }

  // CapillaryPressure-specific interface

  virtual void BatchUpdate( arrayView2d<real64 const> const & phaseVolumeFraction ) override;

  virtual void PointUpdate( arraySlice1d<real64 const> const & phaseVolFraction,
                            localIndex const k,
                            localIndex const q ) override;

  inline static void Compute( localIndex const NP,
                              arraySlice1d<real64  const> const & phaseVolFraction,
                              arraySlice1d<real64> const & phaseCapPressure,
                              arraySlice2d<real64> const & dPhaseCapPressure_dPhaseVolFrac,
                              arraySlice1d<integer const> const & phaseOrder,
                              arraySlice1d<real64  const> const & phaseMinVolumeFraction,
                              arraySlice1d<real64  const> const & phaseCapPressureExponentInv,
                              arraySlice1d<real64  const> const & phaseCapPressureMultiplier,
                              real64 const & capPressureEpsilon,
                              real64 const & volFracScale );

  struct viewKeyStruct : CapillaryPressureBase::viewKeyStruct
  {
    static constexpr auto phaseMinVolumeFractionString      = "phaseMinVolumeFraction";
    static constexpr auto phaseCapPressureExponentInvString = "phaseCapPressureExponentInv";
    static constexpr auto phaseCapPressureMultiplierString  = "phaseCapPressureMultiplier";
    static constexpr auto capPressureEpsilonString = "capPressureEpsilon";
    
    using ViewKey = dataRepository::ViewKey;

    ViewKey phaseMinVolumeFraction      = { phaseMinVolumeFractionString };
    ViewKey phaseCapPressureExponentInv = { phaseCapPressureExponentInvString };
    ViewKey phaseCapPressureMultiplier  = { phaseCapPressureMultiplierString };
    ViewKey capPressureEpsilon          = { capPressureEpsilonString };

  } viewKeysVanGenuchtenCapillaryPressure;

protected:
  virtual void PostProcessInput() override;
  
  inline static void
  EvaluateVanGenuchtenFunction( real64 const & scaledWettingVolFrac,
                                real64 const & dScaledWettingPhaseVolFrac_dVolFrac,
                                real64 & phaseCapPressure,
                                real64 & dPhaseCapPressure_dVolFrac,
                                real64 const & exponentInv,
                                real64 const & multiplier,
                                real64 const & eps );
  
  array1d<real64> m_phaseMinVolumeFraction;
  array1d<real64> m_phaseCapPressureExponentInv;
  array1d<real64> m_phaseCapPressureMultiplier;

  real64 m_capPressureEpsilon;
  real64 m_volFracScale;
};


inline void
VanGenuchtenCapillaryPressure::Compute( localIndex const NP,
                                        arraySlice1d<real64  const> const & phaseVolFraction,
                                        arraySlice1d<real64> const & phaseCapPressure,
                                        arraySlice2d<real64> const & dPhaseCapPressure_dVolFrac,
                                        arraySlice1d<integer const> const & phaseOrder,
                                        arraySlice1d<real64  const> const & phaseMinVolumeFraction,
                                        arraySlice1d<real64  const> const & phaseCapPressureExponentInv,
                                        arraySlice1d<real64  const> const & phaseCapPressureMultiplier,
                                        real64 const & capPressureEpsilon,
                                        real64 const & volFracScale )
{ 

  
  for (localIndex ip = 0; ip < NP; ++ip)
  {
    for (localIndex jp = 0; jp < NP; ++jp)
    {
      dPhaseCapPressure_dVolFrac[ip][jp] = 0.0;
    }
  }

  // the VanGenuchten model does not support volFracScaled = 0 and = 1
  // hence we need an epsilon value to avoid a division by zero
  // TODO: for S < epsilon and S > 1 - epsilon, replace the original unbounded VG curve with a bounded power-law extension
  real64 const eps = capPressureEpsilon;

  real64 const volFracScaleInv = 1.0 / volFracScale;

  
  // compute first water-oil capillary pressure as a function of water-phase vol fraction
  integer const ip_water = phaseOrder[CapillaryPressureBase::PhaseType::WATER];
  if (ip_water >= 0)
  {  

    real64 const volFracScaled = (phaseVolFraction[ip_water] - phaseMinVolumeFraction[ip_water]) * volFracScaleInv;
    real64 const exponentInv   = phaseCapPressureExponentInv[ip_water]; // div by 0 taken care of by initialization check
    real64 const multiplier    = phaseCapPressureMultiplier[ip_water];
      
    real64 const scaledWettingVolFrac                = volFracScaled;
    real64 const dScaledWettingPhaseVolFrac_dVolFrac = volFracScaleInv;
	
    EvaluateVanGenuchtenFunction( scaledWettingVolFrac,
                                  dScaledWettingPhaseVolFrac_dVolFrac,
                                  phaseCapPressure[ip_water],
                                  dPhaseCapPressure_dVolFrac[ip_water][ip_water],
                                  exponentInv,
                                  multiplier,
                                  eps );

  }

  
  // then compute the oil-gas capillary pressure as a function of gas-phase vol fraction
  integer const ip_gas = phaseOrder[CapillaryPressureBase::PhaseType::GAS];
  if (ip_gas >= 0)
  {  
    real64 const volFracScaled = (phaseVolFraction[ip_gas] - phaseMinVolumeFraction[ip_gas]) * volFracScaleInv;
    real64 const exponentInv   = phaseCapPressureExponentInv[ip_gas]; // div by 0 taken care of by initialization check
    real64 const multiplier    = -phaseCapPressureMultiplier[ip_gas]; // for gas capillary pressure, take the opposite of the VG function

    real64 const scaledWettingVolFrac                = 1-volFracScaled;
    real64 const dScaledWettingPhaseVolFrac_dVolFrac =  -volFracScaleInv;
    
    EvaluateVanGenuchtenFunction( scaledWettingVolFrac,
                                  dScaledWettingPhaseVolFrac_dVolFrac,
                                  phaseCapPressure[ip_gas],
                                  dPhaseCapPressure_dVolFrac[ip_gas][ip_gas],
                                  exponentInv,
                                  multiplier,
                                  eps );

  }
}

inline void
VanGenuchtenCapillaryPressure::EvaluateVanGenuchtenFunction( real64 const & scaledWettingVolFrac,
                                                             real64 const & dScaledWettingPhaseVolFrac_dVolFrac,
                                                             real64 & phaseCapPressure,
                                                             real64 & dPhaseCapPressure_dVolFrac,
                                                             real64 const & exponentInv,
                                                             real64 const & multiplier,
                                                             real64 const & eps )
{
  real64 const exponent = 1 / exponentInv; // div by 0 taken care of by initialization check
  
  phaseCapPressure           = 0.0;
  dPhaseCapPressure_dVolFrac = 0.0;

  if (scaledWettingVolFrac >= eps && scaledWettingVolFrac < 1.0-eps)
  {
    // intermediate value
    real64 const a = 1 / std::pow( scaledWettingVolFrac, exponent+1 ); 
    real64 const b = multiplier * std::pow( a * scaledWettingVolFrac - 1, 0.5*(1-exponentInv)-1 );

    phaseCapPressure           = b * ( a * scaledWettingVolFrac - 1 ); // multiplier * ( S_w^(-1/m) - 1 )^( (1-m)/2 )
    dPhaseCapPressure_dVolFrac = - dScaledWettingPhaseVolFrac_dVolFrac * 0.5 * exponent * ( 1 - exponentInv ) * a * b; 
  }
  else // enforce a constant and bounded capillary pressure
  {
    phaseCapPressure = (scaledWettingVolFrac < eps) // div by 0 taken care of by initialization check
                     ? multiplier * std::pow( 1 / std::pow( eps,   exponent ) - 1, 0.5*(1-exponentInv) )
                     : multiplier * std::pow( 1 / std::pow( 1-eps, exponent ) - 1, 0.5*(1-exponentInv) );
  }    
}


} // namespace constitutive

} // namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_VANGENUCHTENCAPILLARYPRESSURE_HPP

