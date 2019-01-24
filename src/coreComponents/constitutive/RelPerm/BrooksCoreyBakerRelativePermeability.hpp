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

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_BROOKSCOREYBAKERRELATIVEPERMEABILITY_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_BROOKSCOREYBAKERRELATIVEPERMEABILITY_HPP

#include "constitutive/RelPerm/RelativePermeabilityBase.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const brooksCoreyBakerRelativePermeability = "BrooksCoreyBakerRelativePermeability";
}
}

namespace constitutive
{

class BrooksCoreyBakerRelativePermeability : public RelativePermeabilityBase
{
public:
  
  BrooksCoreyBakerRelativePermeability( std::string const & name, dataRepository::ManagedGroup * const parent );

  virtual ~BrooksCoreyBakerRelativePermeability() override;

  std::unique_ptr<ConstitutiveBase> DeliverClone( string const & name,
                                                  ManagedGroup * const parent ) const override;

  static std::string CatalogName() { return dataRepository::keys::brooksCoreyBakerRelativePermeability; }

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

  } viewKeysBrooksCoreyBakerRelativePermeability;

protected:
  virtual void PostProcessInput() override;

  static inline void EvaluateBrooksCoreyFunction( real64 const & scaledVolFrac,
                                                  real64 const & dScaledVolFrac_dVolFrac,
                                                  real64 & relPerm,
                                                  real64 & dRelPerm_dVolFrac,
                                                  real64 const & exponent,
                                                  real64 const & maxValue );

  static inline void InterpolateTwoPhaseRelPerms( real64 const & shiftedWaterVolFrac,
						  real64 const & shiftedGasVolFrac,
						  real64 & threePhaseRelPerm,
						  arraySlice1d<real64> const & dThreePhaseRelPerm_dVolFrac,
						  arraySlice1d<integer const> const & phaseOrder,
						  real64 const & relPerm_wo,
						  real64 const & dRelPerm_wo_dOilVolFrac,
						  real64 const & relPerm_go,
						  real64 const & dRelPerm_go_dOilVolFrac);
  
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
BrooksCoreyBakerRelativePermeability::Compute( localIndex const NP,
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
  integer const ip_water = phaseOrder[RelativePermeabilityBase::PhaseType::WATER];
  integer const ip_oil   = phaseOrder[RelativePermeabilityBase::PhaseType::OIL];
  integer const ip_gas   = phaseOrder[RelativePermeabilityBase::PhaseType::GAS];

  real64 oilRelPerm_wo = 0; // oil rel perm using two-phase gas-oil data
  real64 dOilRelPerm_wo_dOilVolFrac = 0; // derivative w.r.t to So
  real64 oilRelPerm_go = 0; // oil rel perm using two-phase gas-oil data
  real64 dOilRelPerm_go_dOilVolFrac = 0; // derivative w.r.t to So

  // this function assumes that the oil phase can always be present (i.e., ip_oil > 0)
  
  // 1) Water and oil phase relative permeabilities using water-oil data
  if (ip_water > 0)
  {
    real64 const scaledWaterVolFrac = (phaseVolFraction[ip_water] - phaseMinVolumeFraction[ip_water]) * volFracScaleInv;
    real64 const scaledOilVolFrac   = (phaseVolFraction[ip_oil]   - phaseMinVolumeFraction[ip_oil])   * volFracScaleInv;
    
    real64 const waterExponent = waterOilRelPermExponent[WaterOilPairPhaseType::WATER];
    real64 const waterMaxValue = waterOilRelPermMaxValue[WaterOilPairPhaseType::WATER];

    // water rel perm
    EvaluateBrooksCoreyFunction( scaledWaterVolFrac,
                                 volFracScaleInv,
			         relPerm[ip_water],
			         dRelPerm_dVolFrac[ip_water][ip_water],
			         waterExponent,
			         waterMaxValue );

    real64 const oilExponent_wo = waterOilRelPermExponent[WaterOilPairPhaseType::OIL];
    real64 const oilMaxValue_wo = waterOilRelPermMaxValue[WaterOilPairPhaseType::OIL];
    
    // oil rel perm
    EvaluateBrooksCoreyFunction( scaledOilVolFrac,
                                 volFracScaleInv,
                                 oilRelPerm_wo,
                                 dOilRelPerm_wo_dOilVolFrac,
                                 oilExponent_wo,
                                 oilMaxValue_wo );

  }
  
  
  // 2) Gas and oil phase relative permeabilities using gas-oil data
  if (ip_gas > 0)
  {
    real64 const scaledGasVolFrac = (phaseVolFraction[ip_gas] - phaseMinVolumeFraction[ip_gas]) * volFracScaleInv;
    real64 const scaledOilVolFrac = (phaseVolFraction[ip_oil] - phaseMinVolumeFraction[ip_oil]) * volFracScaleInv;
	
    real64 const gasExponent = gasOilRelPermExponent[GasOilPairPhaseType::GAS];
    real64 const gasMaxValue = gasOilRelPermMaxValue[GasOilPairPhaseType::GAS];

    // gas rel perm
    EvaluateBrooksCoreyFunction( scaledGasVolFrac,
			         volFracScaleInv,
			         relPerm[ip_gas],
			         dRelPerm_dVolFrac[ip_gas][ip_gas],
			         gasExponent,
			         gasMaxValue );

    real64 const oilExponent_go = gasOilRelPermExponent[GasOilPairPhaseType::OIL];
    real64 const oilMaxValue_go = gasOilRelPermMaxValue[GasOilPairPhaseType::OIL];
    
    // oil rel perm
    EvaluateBrooksCoreyFunction( scaledOilVolFrac,
                                 volFracScaleInv,
                                 oilRelPerm_go,
                                 dOilRelPerm_go_dOilVolFrac,
                                 oilExponent_go,
                                 oilMaxValue_go );

    
  }
  

  // 3) Compute the "three-phase" oil relperm

  // if no gas, use water-oil data
  if (ip_gas < 0)
  {
    relPerm[ip_oil] = oilRelPerm_wo;
    dRelPerm_dVolFrac[ip_oil][ip_oil] = dOilRelPerm_wo_dOilVolFrac;
  }
  // if no water, use gas-oil data
  else if (ip_water < 0)
  {
    relPerm[ip_oil] = oilRelPerm_go;
    dRelPerm_dVolFrac[ip_oil][ip_oil] = dOilRelPerm_go_dOilVolFrac;
  }
  // if water and oil and gas can be present, use saturation-weighted interpolation
  else
  {
    real64 const shiftedWaterVolFrac = (phaseVolFraction[ip_water] - phaseMinVolumeFraction[ip_water]);
    real64 const shiftedGasVolFrac   = (phaseVolFraction[ip_gas]   - phaseMinVolumeFraction[ip_gas]);
    
    InterpolateTwoPhaseRelPerms( shiftedWaterVolFrac,
				 shiftedGasVolFrac,
				 relPerm[ip_oil],
				 dRelPerm_dVolFrac[ip_oil],
				 phaseOrder,
				 oilRelPerm_wo,
				 dOilRelPerm_wo_dOilVolFrac,
				 oilRelPerm_go,
				 dOilRelPerm_go_dOilVolFrac );
  }
}

inline void 
BrooksCoreyBakerRelativePermeability::EvaluateBrooksCoreyFunction(real64 const & scaledVolFrac,
								 real64 const & dScaledVolFrac_dVolFrac,
								 real64 & relPerm,
								 real64 & dRelPerm_dVolFrac,
								 real64 const & exponent,
								 real64 const & maxValue )
{
  relPerm           = 0.0;
  dRelPerm_dVolFrac = 0.0;
  
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

inline void
BrooksCoreyBakerRelativePermeability::InterpolateTwoPhaseRelPerms(real64 const & shiftedWaterVolFrac,
								 real64 const & shiftedGasVolFrac,
								 real64 & threePhaseRelPerm,
								 arraySlice1d<real64> const & dThreePhaseRelPerm_dVolFrac,
								 arraySlice1d<integer const> const & phaseOrder,
								 real64 const & relPerm_wo,
								 real64 const & dRelPerm_wo_dOilVolFrac,
								 real64 const & relPerm_go,
								 real64 const & dRelPerm_go_dOilVolFrac)
{
  integer const ip_water = phaseOrder[RelativePermeabilityBase::PhaseType::WATER];
  integer const ip_oil   = phaseOrder[RelativePermeabilityBase::PhaseType::OIL];
  integer const ip_gas   = phaseOrder[RelativePermeabilityBase::PhaseType::GAS];
  
  // if water phase is immobile, then use the two-phase gas-oil data only
  if (shiftedWaterVolFrac < std::numeric_limits<real64>::epsilon()) 
  {
    threePhaseRelPerm = relPerm_go;
    dThreePhaseRelPerm_dVolFrac[ip_oil] = dRelPerm_go_dOilVolFrac;
  }
  // if gas phase is immobile, then use the two-phase water-oil data only
  else if (shiftedGasVolFrac < std::numeric_limits<real64>::epsilon()) 
  {
    threePhaseRelPerm = relPerm_wo;
    dThreePhaseRelPerm_dVolFrac[ip_oil] = dRelPerm_wo_dOilVolFrac;
  }
  // if both the water phase and the gas phase are mobile,
  // then use a saturation-weighted interpolation of the two-phase oil rel perms
  else
  {
    real64 const sumRelPerm = (shiftedWaterVolFrac * relPerm_wo
			     + shiftedGasVolFrac   * relPerm_go);
    real64 const dSumRelPerm_dWaterVolFrac = relPerm_wo;
    real64 const dSumRelPerm_dOilVolFrac   = shiftedWaterVolFrac * dRelPerm_wo_dOilVolFrac
                                           + shiftedGasVolFrac   * dRelPerm_go_dOilVolFrac;
    real64 const dSumRelPerm_dGasVolFrac   = relPerm_go;

    
    real64 const sumVolFrac    = shiftedWaterVolFrac + shiftedGasVolFrac;
    real64 const sumVolFracInv = 1 / sumVolFrac; // div by 0 handled by the if statement above
    real64 const dSumVolFracInv_dWaterVolFrac = - sumVolFracInv * sumVolFracInv;
    real64 const dSumVolFracInv_dGasVolFrac   = dSumVolFracInv_dWaterVolFrac;

    threePhaseRelPerm = sumRelPerm * sumVolFracInv; // three-phase oil rel perm
    dThreePhaseRelPerm_dVolFrac[ip_water] = dSumRelPerm_dWaterVolFrac * sumVolFracInv  // derivative w.r.t. Sw
                                          + sumRelPerm                * dSumVolFracInv_dWaterVolFrac;
    dThreePhaseRelPerm_dVolFrac[ip_oil]   = dSumRelPerm_dOilVolFrac   * sumVolFracInv; // derivative w.r.t. So
    dThreePhaseRelPerm_dVolFrac[ip_gas]   = dSumRelPerm_dGasVolFrac   * sumVolFracInv  // derivative w.r.t. Sg
				          + sumRelPerm                * dSumVolFracInv_dGasVolFrac;
  }
}

  
} // namespace constitutive

} // namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_BROOKSCOREYBAKERRELATIVEPERMEABILITY_HPP
