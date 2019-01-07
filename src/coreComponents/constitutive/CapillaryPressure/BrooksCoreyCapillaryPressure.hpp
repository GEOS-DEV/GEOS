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
  * @file BrooksCoreyCapillaryPressure.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_BROOKSCOREYCAPILLARYPRESSURE_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_BROOKSCOREYCAPILLARYPRESSURE_HPP

#include "constitutive/CapillaryPressure/CapillaryPressureBase.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const brooksCoreyCapillaryPressure = "BrooksCoreyCapillaryPressure";
}
}

namespace constitutive
{

class BrooksCoreyCapillaryPressure : public CapillaryPressureBase
{
public:

  BrooksCoreyCapillaryPressure( std::string const & name,
				dataRepository::ManagedGroup * const parent );

  virtual ~BrooksCoreyCapillaryPressure() override;

  std::unique_ptr<ConstitutiveBase> DeliverClone( string const & name,
                                                   ManagedGroup * const parent ) const override;

  static std::string CatalogName() { return dataRepository::keys::brooksCoreyCapillaryPressure; }

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
                              arraySlice1d<real64 const> const & phaseMinVolumeFraction,
                              arraySlice1d<real64 const> const & phaseCapPressureExponentInv,
                              arraySlice1d<real64 const> const & phaseEntryPressure,
			      real64 const & capPressureEpsilon,
                              real64 const & satScale );

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

  array1d<real64> m_phaseMinVolumeFraction;
  array1d<real64> m_phaseCapPressureExponentInv;
  array1d<real64> m_phaseEntryPressure;

  real64 m_capPressureEpsilon;
  real64 m_satScale;
};


inline void
BrooksCoreyCapillaryPressure::Compute( localIndex const NP,
				       arraySlice1d<integer const> const & phaseTypes,
                                       arraySlice1d<real64 const> const & phaseVolFraction,
                                       arraySlice1d<real64> const & capPressure,
                                       arraySlice2d<real64> const & dCapPressure_dVolFrac,
                                       arraySlice1d<real64 const> const & phaseMinVolumeFraction,
                                       arraySlice1d<real64 const> const & phaseCapPressureExponentInv,
                                       arraySlice1d<real64 const> const & phaseEntryPressure,
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

  real64 const satScaleInv = 1.0 / satScale;
  
  // the Brooks-Corey model does not support satScaled = 0,
  // hence we need an epsilon value to avoid a division by zero 
  real64 const eps = capPressureEpsilon;

  for (localIndex ip = 0; ip < NP; ++ip)
  {
    real64 const satScaled     = (phaseVolFraction[ip] - phaseMinVolumeFraction[ip]) * satScaleInv;
    real64 const exponent      = 1/phaseCapPressureExponentInv[ip]; // div by 0 taken care of by initialization check
    real64 const entryPressure = phaseEntryPressure[ip];
    
    switch (phaseTypes[ip])
    {

      // compute first water-oil capillary pressure as a function of water-phase vol fraction
      case CapillaryPressureBase::PhaseType::WATER:
      {
	// TODO: put this in a separate function
        if (satScaled >= eps && satScaled < 1.0)
        {
          // intermediate value
          real64 const  val = entryPressure / std::pow( satScaled, exponent + 1);

          capPressure[ip] = val * satScaled;
          dCapPressure_dVolFrac[ip][ip] = - val * exponent * satScaleInv;
        }
        else
        {
          capPressure[ip] = (satScaled < eps)
	                  ? entryPressure / std::pow( eps, exponent ) 
	                  : entryPressure;
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
	// TODO: put this in a separate function
	if (satScaled > 0.0 && satScaled <= 1-eps)
        {
          // intermediate value
          real64 const  val = entryPressure / std::pow( 1-satScaled, exponent + 1);

          capPressure[ip] = val * (1-satScaled);
          dCapPressure_dVolFrac[ip][ip] = val * exponent * satScaleInv;
        }	
        else
        {
          capPressure[ip] = (satScaled > 1-eps)
	                  ? entryPressure / std::pow( eps, exponent ) 
	                  : entryPressure;
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

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_BROOKSCOREYCAPILLARYPRESSURE_HPP

