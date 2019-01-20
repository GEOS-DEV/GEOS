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
  * @file ThreePhaseBakerRelativePermeability.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_THREEPHASERELATIVEPERMEABILITY_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_THREEPHASERELATIVEPERMEABILITY_HPP

#include "constitutive/RelPerm/RelativePermeabilityBase.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const threePhaseBakerRelativePermeability = "ThreePhaseBakerRelativePermeability";
}
}

namespace constitutive
{

class ThreePhaseBakerRelativePermeability : public RelativePermeabilityBase
{
public:
  
  ThreePhaseBakerRelativePermeability( std::string const & name, dataRepository::ManagedGroup * const parent );

  virtual ~ThreePhaseBakerRelativePermeability() override;

  std::unique_ptr<ConstitutiveBase> DeliverClone( string const & name,
                                                  ManagedGroup * const parent ) const override;

  static std::string CatalogName() { return dataRepository::keys::threePhaseBakerRelativePermeability; }

  virtual string GetCatalogName() override { return CatalogName(); }

  // RelPerm-specific interface

  virtual void BatchUpdate( arrayView2d<real64 const> const & phaseVolumeFraction ) override;

  virtual void PointUpdate( arraySlice1d<real64 const> const & phaseVolFraction,
                            localIndex const k,
                            localIndex const q ) override;

  inline static void Compute( localIndex const NP,
                              arraySlice1d<real64 const> const & phaseVolFraction,
                              arraySlice1d<real64> const & phaseRelPerm,
                              arraySlice2d<real64> const & dPhaseRelPerm_dPhaseVolFrac,
			      arraySlice1d<integer const> const & phaseOrder,
                              arraySlice1d<real64  const> const & phaseMinVolumeFraction,
                              arraySlice1d<real64  const> const & waterOilRelPermExponent,
                              arraySlice1d<real64  const> const & waterOilRelPermMaxValue,
                              arraySlice1d<real64  const> const & gasOilRelPermExponent,
                              arraySlice1d<real64  const> const & gasOilRelPermMaxValue,
                              real64 const & volFracScale);
  
  struct viewKeyStruct : RelativePermeabilityBase::viewKeyStruct
  {
    static constexpr auto phaseMinVolumeFractionString  = "phaseMinVolumeFraction";

    static constexpr auto waterOilRelPermExponentString = "waterOilRelPermExponent";
    static constexpr auto waterOilRelPermMaxValueString = "waterOilRelPermMaxValue";
    
    static constexpr auto gasOilRelPermExponentString   = "gasOilRelPermExponent";
    static constexpr auto gasOilRelPermMaxValueString   = "gasOilRelPermMaxValue";

    using ViewKey = dataRepository::ViewKey;

    ViewKey phaseMinVolumeFraction  = { phaseMinVolumeFractionString };
    
    ViewKey waterOilRelPermMaxValue = { waterOilRelPermMaxValueString };
    ViewKey waterOilRelPermExponent = { waterOilRelPermExponentString };

    ViewKey gasOilRelPermMaxValue   = { gasOilRelPermMaxValueString };
    ViewKey gasOilRelPermExponent   = { gasOilRelPermExponentString };

  } viewKeysThreePhaseBakerRelativePermeability;

protected:
  virtual void PostProcessInput() override;

  static inline void EvaluateBrooksCoreyFunction( real64 const & scaledVolFrac,
                                                  real64 const & dScaledVolFrac_dVolFrac,
                                                  real64 & relPerm,
                                                  real64 & dRelPerm_dVolFrac,
                                                  real64 const & exponent,
                                                  real64 const & maxValue );
  
  // order of the phase properties in the water-oil data
  struct WaterOilPairPhaseType
  {
    static constexpr integer WATER = 0; // first water phase property
    static constexpr integer OIL   = 1; // second oil phase property
  };

  // order of the phase properties in the gas-oil data
  struct GasOilPairPhaseType
  {
    static constexpr integer GAS   = 0; // first gas phase property
    static constexpr integer OIL   = 1; // second oil phase property
  };
  
  array1d<real64> m_phaseMinVolumeFraction;

  // water-oil data
  array1d<real64> m_waterOilRelPermExponent;
  array1d<real64> m_waterOilRelPermMaxValue;

  // gas-oil data
  array1d<real64> m_gasOilRelPermExponent;
  array1d<real64> m_gasOilRelPermMaxValue;

  real64 m_volFracScale;
};


inline void
ThreePhaseBakerRelativePermeability::Compute( localIndex const NP,
                                              arraySlice1d<real64 const> const & phaseVolFraction,
                                              arraySlice1d<real64> const & relPerm,
                                              arraySlice2d<real64> const & dRelPerm_dVolFrac,
					      arraySlice1d<integer const> const & phaseOrder,
                                              arraySlice1d<real64  const> const & phaseMinVolumeFraction,
                                              arraySlice1d<real64  const> const & waterOilRelPermExponent,
					      arraySlice1d<real64  const> const & waterOilRelPermMaxValue,
                                              arraySlice1d<real64  const> const & gasOilRelPermExponent,
					      arraySlice1d<real64  const> const & gasOilRelPermMaxValue,
                                              real64 const & volFracScale )
{

  for (localIndex ip = 0; ip < NP; ++ip)
  {
    for (localIndex jp = 0; jp < NP; ++jp)
    {
      dRelPerm_dVolFrac[ip][jp] = 0.0;
    }
  }
  
  real64 const volFracScaleInv = 1.0 / volFracScale;
  
  // 1) Water phase relative permeability using water-oil data as a function of Sw
  integer const ip_water = phaseOrder[RelativePermeabilityBase::PhaseType::WATER];

  real64 const shiftedWaterVolFrac = (phaseVolFraction[ip_water] - phaseMinVolumeFraction[ip_water]);
  real64 const scaledWaterVolFrac  = shiftedWaterVolFrac * volFracScaleInv;
  real64 const waterExponent = waterOilRelPermExponent[WaterOilPairPhaseType::WATER];
  real64 const waterMaxValue = waterOilRelPermMaxValue[WaterOilPairPhaseType::WATER];
  
  EvaluateBrooksCoreyFunction( scaledWaterVolFrac,
			       volFracScaleInv,
			       relPerm[ip_water],
			       dRelPerm_dVolFrac[ip_water][ip_water],
			       waterExponent,
			       waterMaxValue );

  
  
  // 2) Gas phase relative permeability using gas-oil data as a function of Sg
  integer const ip_gas = phaseOrder[RelativePermeabilityBase::PhaseType::GAS];

  real64 const shiftedGasVolFrac = (phaseVolFraction[ip_gas] - phaseMinVolumeFraction[ip_gas]);
  real64 const scaledGasVolFrac  = shiftedGasVolFrac * volFracScaleInv;
  real64 const gasExponent = gasOilRelPermExponent[GasOilPairPhaseType::GAS];
  real64 const gasMaxValue = gasOilRelPermMaxValue[GasOilPairPhaseType::GAS];

  EvaluateBrooksCoreyFunction( scaledGasVolFrac,
			       volFracScaleInv,
			       relPerm[ip_gas],
			       dRelPerm_dVolFrac[ip_gas][ip_gas],
			       gasExponent,
			       gasMaxValue );

  // 3) Oil phase relative permeability as a function of Sw and Sg
  integer const ip_oil = phaseOrder[RelativePermeabilityBase::PhaseType::OIL];
  real64 const scaledOilVolFrac = (phaseVolFraction[ip_oil] - phaseMinVolumeFraction[ip_oil]) * volFracScaleInv;
  /*
   * Below we use the notations:
   *    "_wo" suffix refers to two-phase Water-Oil data
   *    "_go" suffix refers to two-phase Gas-Oil data
   */
  
  // 3.1) Compute oil phase relative permeability as a function of So using water-oil data
  real64 const oilExponent_wo = waterOilRelPermExponent[WaterOilPairPhaseType::OIL];
  real64 const oilMaxValue_wo = waterOilRelPermMaxValue[WaterOilPairPhaseType::OIL];

  real64 oilRelPerm_wo = 0; // oil rel perm using two-phase water-oil data
  real64 dOilRelPerm_wo_dOilVolFrac = 0; // derivative w.r.t to So

  EvaluateBrooksCoreyFunction( scaledOilVolFrac,
			       volFracScaleInv,
			       oilRelPerm_wo, 
			       dOilRelPerm_wo_dOilVolFrac,
			       oilExponent_wo,
			       oilMaxValue_wo );
  
  // 3.2) Compute oil phase relative permeability as a function of Sg using gas-oil data
  real64 const oilExponent_go = gasOilRelPermExponent[GasOilPairPhaseType::OIL];
  real64 const oilMaxValue_go = gasOilRelPermMaxValue[GasOilPairPhaseType::OIL];

  real64 oilRelPerm_go = 0; // oil rel perm using two-phase gas-oil data
  real64 dOilRelPerm_go_dOilVolFrac = 0; // derivative w.r.t to So
  
  EvaluateBrooksCoreyFunction( scaledOilVolFrac,
			       volFracScaleInv,
			       oilRelPerm_go,
			       dOilRelPerm_go_dOilVolFrac,
			       oilExponent_go,
			       oilMaxValue_go );

  
  // 3.3) Compute the saturation-weighted interpolation

  // if water phase is immobile, then use the two-phase gas-oil data only
  if (shiftedWaterVolFrac < std::numeric_limits<real64>::epsilon()) 
  {
    relPerm[ip_oil] = oilRelPerm_go;
    dRelPerm_dVolFrac[ip_oil][ip_oil] = dOilRelPerm_go_dOilVolFrac;
  }
  // if gas phase is immobile, then use the two-phase water-oil data only
  else if (shiftedGasVolFrac < std::numeric_limits<real64>::epsilon()) 
  {
    relPerm[ip_oil] = oilRelPerm_wo;
    dRelPerm_dVolFrac[ip_oil][ip_oil] = dOilRelPerm_wo_dOilVolFrac;
  }
  // if both the water phase and the gas phase are mobile,
  // then use a saturation-weighted interpolation of the two-phase oil rel perms
  else
  {
    real64 const sumOilRelPerm = (shiftedWaterVolFrac * oilRelPerm_wo
			        + shiftedGasVolFrac   * oilRelPerm_go);
    real64 const dSumOilRelPerm_dWaterVolFrac = oilRelPerm_wo;
    real64 const dSumOilRelPerm_dOilVolFrac   = shiftedWaterVolFrac * dOilRelPerm_wo_dOilVolFrac
                                              + shiftedGasVolFrac   * dOilRelPerm_go_dOilVolFrac;
    real64 const dSumOilRelPerm_dGasVolFrac   = oilRelPerm_go;

    
    real64 const sumVolFrac    = shiftedWaterVolFrac + shiftedGasVolFrac;
    real64 const sumVolFracInv = 1 / sumVolFrac; // div by 0 handled by the if statement above
    real64 const dSumVolFracInv_dWaterVolFrac = - sumVolFracInv * sumVolFracInv;
    real64 const dSumVolFracInv_dGasVolFrac   = dSumVolFracInv_dWaterVolFrac;

    relPerm[ip_oil] = sumOilRelPerm * sumVolFracInv; // three-phase oil rel perm
    dRelPerm_dVolFrac[ip_oil][ip_water] = dSumOilRelPerm_dWaterVolFrac * sumVolFracInv  // derivative w.r.t. Sw
					+ sumOilRelPerm                * dSumVolFracInv_dWaterVolFrac;
    dRelPerm_dVolFrac[ip_oil][ip_oil]   = dSumOilRelPerm_dOilVolFrac   * sumVolFracInv; // derivative w.r.t. So
    dRelPerm_dVolFrac[ip_oil][ip_gas]   = dSumOilRelPerm_dGasVolFrac   * sumVolFracInv  // derivative w.r.t. Sg
					+ sumOilRelPerm                * dSumVolFracInv_dGasVolFrac;
  }
}

inline void 
ThreePhaseBakerRelativePermeability::EvaluateBrooksCoreyFunction(real64 const & scaledVolFrac,
								 real64 const & dScaledVolFrac_dVolFrac,
								 real64 & relPerm,
								 real64 & dRelPerm_dVolFrac,
								 real64 const & exponent,
								 real64 const & maxValue )
{
  if (scaledVolFrac > 0.0 && scaledVolFrac < 1.0)
  {
    // intermediate value
    real64 const v = maxValue * std::pow( scaledVolFrac, exponent - 1.0 );

    relPerm = v * scaledVolFrac;
    dRelPerm_dVolFrac = v * exponent * dScaledVolFrac_dVolFrac;
  }
  else
  {
    relPerm = (scaledVolFrac < 0.0) ? 0.0 : maxValue;
  }
}
  
} // namespace constitutive

} // namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_BROOKSCOREYRELATIVEPERMEABILITY_HPP
