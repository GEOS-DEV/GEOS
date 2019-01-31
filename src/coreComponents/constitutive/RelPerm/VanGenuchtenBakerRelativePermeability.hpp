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
  * @file VanGenuchtenBakerRelativePermeability.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_VANGENUCHTENBAKERRELATIVEPERMEABILITY_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_VANGENUCHTENBAKERRELATIVEPERMEABILITY_HPP

#include "constitutive/RelPerm/RelativePermeabilityBase.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const vanGenuchtenBakerRelativePermeability = "VanGenuchtenBakerRelativePermeability";
}
}

namespace constitutive
{

class VanGenuchtenBakerRelativePermeability : public RelativePermeabilityBase
{
public:
  
  VanGenuchtenBakerRelativePermeability( std::string const & name, dataRepository::ManagedGroup * const parent );

  virtual ~VanGenuchtenBakerRelativePermeability() override;

  std::unique_ptr<ConstitutiveBase> DeliverClone( string const & name,
                                                  ManagedGroup * const parent ) const override;

  static std::string CatalogName() { return dataRepository::keys::vanGenuchtenBakerRelativePermeability; }

  virtual string GetCatalogName() override { return CatalogName(); }

  // RelPerm-specific interface

  /**
   * @brief Function to update state of a single material point.
   * @param[in] phaseVolFraction input phase volume fraction
   * @param[in] k the first index of the storage arrays (elem index)
   * @param[in] q the secound index of the storage arrays (quadrature index)
   *
   * @note This function performs a point update, but should not be called
   *       within a kernel since it is virtual, and the required data is not
   *       guaranteed to be in the target memory space.
   */
  virtual void BatchUpdate( arrayView2d<real64 const> const & phaseVolumeFraction ) override;

  /**
   * @brief Perform a batch constitutive update (all points).
   * @param[in] phaseVolFraction input phase volume fraction
   */
  virtual void PointUpdate( arraySlice1d<real64 const> const & phaseVolFraction,
                            localIndex const k,
                            localIndex const q ) override;

  
  /**
   * @brief Computes the phase relative permeabilities using the Van Genuchten-Baker method
   * @param NP phase index
   * @param[in] phaseVolFraction vector of phase volume fractions
   * @param[out] phaseRelPerm the computed relative permeability value vector for all phases
   * @param[out] dPhaseRelPerm_dPhaseVolFrac the computed partial derivative of the relative wrt to the volume fraction of the phases
   * @param[in] phaseOrder vector of phase orders
   * @param[in] phaseMinVolumeFraction vector of minimum phase volume fractions
   * @param[in] waterOilRelPermExponentInv vector of exponents used in the computation of the water-oil relative permeability
   * @param[in] waterOilRelPermMaxValue vector of water-oil permeability curve end-point values 
   * @param[in] gasOilRelPermExponentInv vector of exponents used in the computation of the gas-oil relative permeability
   * @param[in] gasOilRelPermMaxValue vector of gas-oil permeability curve end-point values
   * @param[in] volFracScale scaling factor to apply to the entire relative permeability curve
   * @return (void)
   *
   * This function computes an entire relative permeability curve based on the Van Genuchten-Baker method
   * Reference: Eclipse technical description
   */
  inline static void Compute( localIndex const NP,
                              arraySlice1d<real64 const> const & phaseVolFraction,
                              arraySlice1d<real64> const & phaseRelPerm,
                              arraySlice2d<real64> const & dPhaseRelPerm_dPhaseVolFrac,
                              arraySlice1d<integer const> const & phaseOrder,
                              arraySlice1d<real64  const> const & phaseMinVolumeFraction,
                              arraySlice1d<real64  const> const & waterOilRelPermExponentInv,
                              arraySlice1d<real64  const> const & waterOilRelPermMaxValue,
                              arraySlice1d<real64  const> const & gasOilRelPermExponentInv,
                              arraySlice1d<real64  const> const & gasOilRelPermMaxValue,
                              real64 const & volFracScale);
  
  struct viewKeyStruct : RelativePermeabilityBase::viewKeyStruct
  {
    static constexpr auto phaseMinVolumeFractionString  = "phaseMinVolumeFraction";

    static constexpr auto waterOilRelPermExponentInvString = "waterOilRelPermExponentInv";
    static constexpr auto waterOilRelPermMaxValueString    = "waterOilRelPermMaxValue";
    
    static constexpr auto gasOilRelPermExponentInvString = "gasOilRelPermExponentInv";
    static constexpr auto gasOilRelPermMaxValueString    = "gasOilRelPermMaxValue";

    using ViewKey = dataRepository::ViewKey;

    ViewKey phaseMinVolumeFraction = { phaseMinVolumeFractionString };
    
    ViewKey waterOilRelPermMaxValue    = { waterOilRelPermMaxValueString };
    ViewKey waterOilRelPermExponentInv = { waterOilRelPermExponentInvString };

    ViewKey gasOilRelPermMaxValue    = { gasOilRelPermMaxValueString };
    ViewKey gasOilRelPermExponentInv = { gasOilRelPermExponentInvString };

  } viewKeysVanGenuchtenBakerRelativePermeability;

protected:
  virtual void PostProcessInput() override;


  /**
   * @brief Evaluate the Van Genuchten relperm function for a given (scalar) phase saturation
   * @param[in] scaledVolFrac the scaled volume fraction for this phase
   * @param[in] dScaledVolFrac_dVolFrac the derivative of scaled volume fraction for this phase wrt to the volume fraction
   * @param[out] relperm the relative permeability for this phase
   * @param[out] dRelPerm_dVolFrac the derivative of the relative permeability wrt to the volume fraction of the phase
   * @param[in] exponentInv the inverse of the exponent used in the VG model
   * @param[in] maxValue the endpoint relative permeability value
   * @return (void)
   *
   * This function evaluates the relperm function and its derivative at a given phase saturation 
   * Reference: Eclipse technical description and Petrowiki
   */
  static inline void EvaluateVanGenuchtenFunction( real64 const & scaledVolFrac,
                                                   real64 const & dScaledVolFrac_dVolFrac,
                                                   real64 & relPerm,
                                                   real64 & dRelPerm_dVolFrac,
                                                   real64 const & exponentInv,
                                                   real64 const & maxValue );

  /**
   * @brief Interpolate the two-phase relperms to compute the three-phase relperm 
   * @param[in] shiftedWaterVolFrac
   * @param[in] gasVolFrac
   * @param[out] threePhaseRelPerm
   * @param[out] dThreePhaseRelPerm_dVolFrac
   * @param[in] relPerm_wo
   * @param[in] dRelPerm_wo_dOilVolFrac
   * @param[in] relPerm_go
   * @param[in] dRelPerm_go_dOilVolFrac
   * @return (void)
   *
   * This function interpolates the two-phase relperms to compute the three-phase relperm
   * The interpolation is based on the modified Baker method, also used as default in Eclipse
   * Reference: Eclipse technical description
   */
  static inline void InterpolateTwoPhaseRelPerms( real64 const & shiftedWaterVolFrac,
                                                  real64 const & gasVolFrac,
                                                  real64 & threePhaseRelPerm,
                                                  arraySlice1d<real64> const & dThreePhaseRelPerm_dVolFrac,
                                                  arraySlice1d<integer const> const & phaseOrder,
                                                  real64 const & relPerm_wo,
                                                  real64 const & dRelPerm_wo_dOilVolFrac,
                                                  real64 const & relPerm_go,
                                                  real64 const & dRelPerm_go_dOilVolFrac );
  
  array1d<real64> m_phaseMinVolumeFraction;

  // water-oil data
  array1d<real64> m_waterOilRelPermExponentInv;
  array1d<real64> m_waterOilRelPermMaxValue;

  // gas-oil data
  array1d<real64> m_gasOilRelPermExponentInv;
  array1d<real64> m_gasOilRelPermMaxValue;

  real64 m_volFracScale;
};


inline void
VanGenuchtenBakerRelativePermeability::Compute( localIndex const NP,
                                                arraySlice1d<real64 const> const & phaseVolFraction,
                                                arraySlice1d<real64> const & relPerm,
                                                arraySlice2d<real64> const & dRelPerm_dVolFrac,
                                                arraySlice1d<integer const> const & phaseOrder,
                                                arraySlice1d<real64  const> const & phaseMinVolumeFraction,
                                                arraySlice1d<real64  const> const & waterOilRelPermExponentInv,
                                                arraySlice1d<real64  const> const & waterOilRelPermMaxValue,
                                                arraySlice1d<real64  const> const & gasOilRelPermExponentInv,
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
  integer const ip_water = phaseOrder[PhaseType::WATER];
  integer const ip_oil   = phaseOrder[PhaseType::OIL];
  integer const ip_gas   = phaseOrder[PhaseType::GAS];

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
    
    real64 const waterExponentInv = waterOilRelPermExponentInv[WaterOilPairPhaseType::WATER];
    real64 const waterMaxValue = waterOilRelPermMaxValue[WaterOilPairPhaseType::WATER];

    // water rel perm
    EvaluateVanGenuchtenFunction( scaledWaterVolFrac,
                                  volFracScaleInv,
                                  relPerm[ip_water],
                                  dRelPerm_dVolFrac[ip_water][ip_water],
                                  waterExponentInv,
                                  waterMaxValue );

    real64 const oilExponentInv_wo = waterOilRelPermExponentInv[WaterOilPairPhaseType::OIL];
    real64 const oilMaxValue_wo = waterOilRelPermMaxValue[WaterOilPairPhaseType::OIL];
    
    // oil rel perm
    EvaluateVanGenuchtenFunction( scaledOilVolFrac,
                                  volFracScaleInv,
                                  oilRelPerm_wo,
                                  dOilRelPerm_wo_dOilVolFrac,
                                  oilExponentInv_wo,
                                  oilMaxValue_wo );

  }
  
  
  // 2) Gas and oil phase relative permeabilities using gas-oil data
  if (ip_gas > 0)
  {
    real64 const scaledGasVolFrac = (phaseVolFraction[ip_gas] - phaseMinVolumeFraction[ip_gas]) * volFracScaleInv;
    real64 const scaledOilVolFrac = (phaseVolFraction[ip_oil] - phaseMinVolumeFraction[ip_oil]) * volFracScaleInv;
	
    real64 const gasExponentInv = gasOilRelPermExponentInv[GasOilPairPhaseType::GAS];
    real64 const gasMaxValue = gasOilRelPermMaxValue[GasOilPairPhaseType::GAS];

    // gas rel perm
    EvaluateVanGenuchtenFunction( scaledGasVolFrac,
                                  volFracScaleInv,
                                  relPerm[ip_gas],
                                  dRelPerm_dVolFrac[ip_gas][ip_gas],
                                  gasExponentInv,
                                  gasMaxValue );

    real64 const oilExponentInv_go = gasOilRelPermExponentInv[GasOilPairPhaseType::OIL];
    real64 const oilMaxValue_go    = gasOilRelPermMaxValue[GasOilPairPhaseType::OIL];
    
    // oil rel perm
    EvaluateVanGenuchtenFunction( scaledOilVolFrac,
                                  volFracScaleInv,
                                  oilRelPerm_go,
                                  dOilRelPerm_go_dOilVolFrac,
                                  oilExponentInv_go,
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
    
    InterpolateTwoPhaseRelPerms( shiftedWaterVolFrac,
				 phaseVolFraction[ip_gas],
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
VanGenuchtenBakerRelativePermeability::EvaluateVanGenuchtenFunction( real64 const & scaledVolFrac,
                                                                     real64 const & dScaledVolFrac_dVolFrac,
                                                                     real64 & relPerm,
                                                                     real64 & dRelPerm_dVolFrac,
                                                                     real64 const & exponentInv,
                                                                     real64 const & maxValue )
{
  real64 const exponent = 1 / exponentInv;
  
  relPerm           = 0.0;
  dRelPerm_dVolFrac = 0.0;
  
  if (scaledVolFrac > 0.0 && scaledVolFrac < 1.0)
  {
    // intermediate values
    real64 const a = std::pow( scaledVolFrac, exponent-1 ); 
    real64 const b = std::pow( 1 - a * scaledVolFrac, exponentInv-1 );
    real64 const c = ( 1 - b * ( 1 - a * scaledVolFrac ) );
    real64 const volFracSquared = scaledVolFrac * scaledVolFrac;
    real64 const dVolFracSquared_dVolFrac = 2 * dScaledVolFrac_dVolFrac * scaledVolFrac;
    
    relPerm  = maxValue * volFracSquared * c;

    dRelPerm_dVolFrac  = dVolFracSquared_dVolFrac * c
                       + volFracSquared * dScaledVolFrac_dVolFrac * a * b;
    dRelPerm_dVolFrac *= maxValue;
  }
  else
  {
    relPerm = (scaledVolFrac < 0.0) ? 0.0 : maxValue;
  }
}

inline void
VanGenuchtenBakerRelativePermeability::InterpolateTwoPhaseRelPerms( real64 const & shiftedWaterVolFrac,
                                                                    real64 const & gasVolFrac,
                                                                    real64 & threePhaseRelPerm,
                                                                    arraySlice1d<real64> const & dThreePhaseRelPerm_dVolFrac,
                                                                    arraySlice1d<integer const> const & phaseOrder,
                                                                    real64 const & relPerm_wo,
                                                                    real64 const & dRelPerm_wo_dOilVolFrac,
                                                                    real64 const & relPerm_go,
                                                                    real64 const & dRelPerm_go_dOilVolFrac )
{
  integer const ip_water = phaseOrder[PhaseType::WATER];
  integer const ip_oil   = phaseOrder[PhaseType::OIL];
  integer const ip_gas   = phaseOrder[PhaseType::GAS];
  
  // if water phase is immobile, then use the two-phase gas-oil data only
  if (shiftedWaterVolFrac < std::numeric_limits<real64>::epsilon()) 
  {
    threePhaseRelPerm = relPerm_go;
    dThreePhaseRelPerm_dVolFrac[ip_oil] = dRelPerm_go_dOilVolFrac;
  }
  // if gas phase is immobile, then use the two-phase water-oil data only
  else if (gasVolFrac < std::numeric_limits<real64>::epsilon()) 
  {
    threePhaseRelPerm = relPerm_wo;
    dThreePhaseRelPerm_dVolFrac[ip_oil] = dRelPerm_wo_dOilVolFrac;
  }
  // if both the water phase and the gas phase are mobile,
  // then use a saturation-weighted interpolation of the two-phase oil rel perms
  else
  {
    real64 const sumRelPerm = (shiftedWaterVolFrac * relPerm_wo
			     + gasVolFrac   * relPerm_go);
    real64 const dSumRelPerm_dWaterVolFrac = relPerm_wo;
    real64 const dSumRelPerm_dOilVolFrac   = shiftedWaterVolFrac * dRelPerm_wo_dOilVolFrac
                                           + gasVolFrac   * dRelPerm_go_dOilVolFrac;
    real64 const dSumRelPerm_dGasVolFrac   = relPerm_go;

    
    real64 const sumVolFrac    = shiftedWaterVolFrac + gasVolFrac;
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

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_VANGENUCHTENBAKERRELATIVEPERMEABILITY_HPP
