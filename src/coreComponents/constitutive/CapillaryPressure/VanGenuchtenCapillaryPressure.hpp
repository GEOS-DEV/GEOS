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

  virtual void ProcessInputFile_PostProcess() override;

  // CapillaryPressure-specific interface

  virtual void BatchUpdate( arrayView2d<real64 const> const & phaseVolumeFraction ) override;

  virtual void PointUpdate( arraySlice1d<real64 const> const & phaseVolFraction,
                            localIndex const k,
                            localIndex const q ) override;

  inline static void Compute( localIndex const NP,
			      arraySlice1d<integer const> const & phaseTypes,
                              arraySlice1d<real64 const> const & phaseVolFraction,
                              arraySlice1d<real64> const & phaseCapPressure,
                              arraySlice2d<real64> const & dPhaseCapPressure_dPhaseVolFrac,
                              arraySlice1d<real64 const> const &  phaseMinVolumeFraction,
                              arraySlice1d<real64 const> const & phaseCapPressureExponentInv,
                              arraySlice1d<real64 const> const & phaseCapPressureMultiplier,
			      real64 const & capPressureEpsilon,
                              real64 const & satScale );

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

  array1d<real64> m_phaseMinVolumeFraction;
  array1d<real64> m_phaseCapPressureExponentInv;
  array1d<real64> m_phaseCapPressureMultiplier;

  real64 m_capPressureEpsilon;
  real64 m_satScale;
};


inline void
VanGenuchtenCapillaryPressure::Compute( localIndex const NP,
					arraySlice1d<integer const> const & phaseTypes,
                                        arraySlice1d<real64 const> const & phaseVolFraction,
                                        arraySlice1d<real64> const & capPressure,
                                        arraySlice2d<real64> const & dCapPressure_dVolFrac,
                                        arraySlice1d<real64 const> const & phaseMinVolumeFraction,
                                        arraySlice1d<real64 const> const & phaseCapPressureExponentInv,
                                        arraySlice1d<real64 const> const & phaseCapPressureMultiplier,
					real64 const & capPressureEpsilon,
                                        real64 const & satScale )
{ 

  for (localIndex ip = 0; ip < NP; ++ip)
  {
    for (localIndex jp = 0; jp < NP; ++jp)
    {
      dCapPressure_dVolFrac[ip][jp] = 0.0;
    }
  }

  // the VanGenuchten model does not support satScaled = 0 and = 1
  // hence we need an epsilon value to avoid a division by zero 
  real64 const eps = capPressureEpsilon;

  real64 const satScaleInv = 1.0 / satScale;

  for (localIndex ip = 0; ip < NP; ++ip)
  {
    real64 const satScaled   = (phaseVolFraction[ip] - phaseMinVolumeFraction[ip]) * satScaleInv;
    real64 const exponentInv = phaseCapPressureExponentInv[ip]; // div by 0 taken care of by initialization check
    real64 const exponent    = 1/phaseCapPressureExponentInv[ip]; 
    real64 const multiplier  = phaseCapPressureMultiplier[ip];

    switch (phaseTypes[ip])
    {

      // compute first water-oil capillary pressure as a function of water-phase vol fraction
      case CapillaryPressureBase::PhaseType::WATER:
      {
	// TODO: put this in a separate function
        if (satScaled >= eps && satScaled < 1.0-eps)
        {
	  // intermediate value
	  real64 const a = 1 / std::pow( satScaled, exponent+1 ); // S_c^(-1/m-1)
	  real64 const b = multiplier * std::pow( a * satScaled - 1, 0.5*(1-exponentInv)-1 );

	  capPressure[ip] = b * ( a * satScaled - 1 );
	  dCapPressure_dVolFrac[ip][ip] = - 0.5 * exponent * ( 1 - exponentInv ) * satScaleInv * a * b; 
        }
	else
	{
	  capPressure[ip] = (satScaled < eps)
	                  ? multiplier * std::pow( 1 / std::pow( eps,   exponent ) - 1, 0.5*(1-exponentInv) )
	                  : multiplier * std::pow( 1 / std::pow( 1-eps, exponent ) - 1, 0.5*(1-exponentInv) );
	}    
        break;
      }

      // no capillary pressure for the oil phase
      case CapillaryPressureBase::PhaseType::OIL:
      {
        capPressure[ip] = 0;
        break;
      }
      
      // then compute the oil-gas capillary pressure as a function of gas-phase vol fraction
      case CapillaryPressureBase::PhaseType::GAS:
      {
	// TODO: put this in a separate functions
        if (satScaled > eps && satScaled <= 1.0-eps)
        {
	  // intermediate value
	  real64 const a = 1 / std::pow( 1-satScaled, exponent+1 ); // S_c^(-1/m-1)
	  real64 const b = multiplier * std::pow( a * (1-satScaled) - 1, 0.5*(1-exponentInv)-1 );

	  capPressure[ip] = b * ( a * (1-satScaled) - 1 );
	  dCapPressure_dVolFrac[ip][ip] = 0.5 * exponent * ( 1 - exponentInv ) * satScaleInv * a * b; 
        }
	else
	{
	  capPressure[ip] = (satScaled <= eps)
	                  ? multiplier * std::pow( 1 / std::pow( 1-eps, exponent ) - 1, 0.5*(1-exponentInv) )
	                  : multiplier * std::pow( 1 / std::pow( eps,   exponent ) - 1, 0.5*(1-exponentInv) );
	}
        break;	
      }

      default:
        GEOS_ERROR("Unsupported phase type");
    }
  }
}

} // namespace constitutive

} // namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_VANGENUCHTENCAPILLARYPRESSURE_HPP

