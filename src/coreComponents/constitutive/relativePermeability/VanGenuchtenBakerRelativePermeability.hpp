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
 * @file VanGenuchtenBakerRelativePermeability.hpp
 */

#ifndef GEOS_CONSTITUTIVE_VANGENUCHTENBAKERRELATIVEPERMEABILITY_HPP
#define GEOS_CONSTITUTIVE_VANGENUCHTENBAKERRELATIVEPERMEABILITY_HPP

#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityInterpolators.hpp"

namespace geos
{
namespace constitutive
{

class VanGenuchtenBakerRelativePermeabilityUpdate final : public RelativePermeabilityBaseUpdate
{
public:

  VanGenuchtenBakerRelativePermeabilityUpdate( arrayView2d< real64 const > const & phaseMinVolumeFraction,
                                               arrayView2d< real64 const > const & waterOilRelPermExponentInv,
                                               arrayView2d< real64 const > const & waterOilRelPermMaxValue,
                                               arrayView2d< real64 const > const & gasOilRelPermExponentInv,
                                               arrayView2d< real64 const > const & gasOilRelPermMaxValue,
                                               arrayView1d< real64 const > const & volFracScale,
                                               arrayView1d< integer const > const & phaseTypes,
                                               arrayView1d< integer const > const & phaseOrder,
                                               arrayView4d< real64, relperm::USD_RELPERM > const & phaseRelPerm,
                                               arrayView5d< real64, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac,
                                               arrayView3d< real64, relperm::USD_PHASE > const & phaseTrappedVolFrac )
    : RelativePermeabilityBaseUpdate( phaseTypes,
                                      phaseOrder,
                                      phaseRelPerm,
                                      dPhaseRelPerm_dPhaseVolFrac,
                                      phaseTrappedVolFrac ),
    m_phaseMinVolumeFraction( phaseMinVolumeFraction ),
    m_waterOilRelPermExponentInv( waterOilRelPermExponentInv ),
    m_waterOilRelPermMaxValue( waterOilRelPermMaxValue ),
    m_gasOilRelPermExponentInv( gasOilRelPermExponentInv ),
    m_gasOilRelPermMaxValue( gasOilRelPermMaxValue ),
    m_volFracScale( volFracScale )
  {}

  GEOS_HOST_DEVICE
  void compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                arraySlice1d< real64, relperm::USD_PHASE - 2 > const & phaseTrappedVolFrac,
                arraySlice2d< real64, relperm::USD_RELPERM - 2 > const & phaseRelPerm,
                arraySlice3d< real64, relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const;

  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const override
  {
    compute( phaseVolFraction,
             m_phaseTrappedVolFrac[k][q],
             m_phaseRelPerm[k][q],
             m_dPhaseRelPerm_dPhaseVolFrac[k][q] );
  }

private:

  /**
   * @brief Evaluate the Van Genuchten relperm function for a given (scalar) phase saturation
   * @param[in] scaledVolFrac the scaled volume fraction for this phase
   * @param[in] dScaledVolFrac_dVolFrac the derivative of scaled volume fraction for this phase wrt to the volume
   * fraction
   * @param[out] relperm the relative permeability for this phase
   * @param[out] dRelPerm_dVolFrac the derivative of the relative permeability wrt to the volume fraction of the phase
   * @param[in] exponentInv the inverse of the exponent used in the VG model
   * @param[in] maxValue the endpoint relative permeability value
   * @return (void)
   *
   * This function evaluates the relperm function and its derivative at a given phase saturation
   * Reference: Eclipse technical description and Petrowiki
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void evaluateVanGenuchtenFunction( real64 const & scaledVolFrac,
                                            real64 const & dScaledVolFrac_dVolFrac,
                                            real64 const & exponentInv,
                                            real64 const & maxValue,
                                            real64 & relPerm,
                                            real64 & dRelPerm_dVolFrac );

  arrayView2d< real64 const > m_phaseMinVolumeFraction;

  arrayView2d< real64 const > m_waterOilRelPermExponentInv;
  arrayView2d< real64 const > m_waterOilRelPermMaxValue;

  arrayView2d< real64 const > m_gasOilRelPermExponentInv;
  arrayView2d< real64 const > m_gasOilRelPermMaxValue;

  arrayView1d< real64 const > m_volFracScale;

};

class VanGenuchtenBakerRelativePermeability : public RelativePermeabilityBase
{
public:

  VanGenuchtenBakerRelativePermeability( string const & name, dataRepository::Group * const parent );

  static string catalogName() { return "VanGenuchtenBakerRelativePermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = VanGenuchtenBakerRelativePermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  struct viewKeyStruct : RelativePermeabilityBase::viewKeyStruct
  {
    static constexpr char const * phaseMinVolumeFractionString() { return "phaseMinVolumeFraction"; }
    static constexpr char const * waterOilRelPermExponentInvString() { return "waterOilRelPermExponentInv"; }
    static constexpr char const * waterOilRelPermMaxValueString() { return "waterOilRelPermMaxValue"; }
    static constexpr char const * gasOilRelPermExponentInvString() { return "gasOilRelPermExponentInv"; }
    static constexpr char const * gasOilRelPermMaxValueString() { return "gasOilRelPermMaxValue"; }
    static constexpr char const * volFracScaleString() { return "volFracScale"; }
  };

  arrayView2d< real64 const > getPhaseMinVolumeFraction() const override { return m_phaseMinVolumeFraction; };

protected:

  virtual void postProcessInput() override;

  array2d< real64 > m_phaseMinVolumeFraction;

  // water-oil data
  array2d< real64 > m_waterOilRelPermExponentInv;
  array2d< real64 > m_waterOilRelPermMaxValue;

  // gas-oil data
  array2d< real64 > m_gasOilRelPermExponentInv;
  array2d< real64 > m_gasOilRelPermMaxValue;

  array1d< real64 > m_volFracScale;
};


GEOS_HOST_DEVICE
inline void
VanGenuchtenBakerRelativePermeabilityUpdate::
  compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
           arraySlice1d< real64, relperm::USD_PHASE - 2 > const & phaseTrappedVolFrac,
           arraySlice2d< real64, relperm::USD_RELPERM - 2 > const & phaseRelPerm,
           arraySlice3d< real64, relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const
{
  LvArray::forValuesInSlice( dPhaseRelPerm_dPhaseVolFrac, []( real64 & val ){ val = 0.0; } );

  using PT = RelativePermeabilityBase::PhaseType;
  integer const ipWater = m_phaseOrder[PT::WATER];
  integer const ipOil   = m_phaseOrder[PT::OIL];
  integer const ipGas  = m_phaseOrder[PT::GAS];

  for( int dir = 0; dir < 3; ++dir )
  {
    real64 const volFracScaleInv = 1.0 / m_volFracScale[dir];

    real64 oilRelPerm_wo = 0.0; // oil rel perm using two-phase gas-oil data
    real64 dOilRelPerm_wo_dOilVolFrac = 0.0; // derivative w.r.t to So
    real64 oilRelPerm_go = 0.0; // oil rel perm using two-phase gas-oil data
    real64 dOilRelPerm_go_dOilVolFrac = 0.0; // derivative w.r.t to So

    // this function assumes that the oil phase can always be present (i.e., ipOil > 0)

    // 1) Water and oil phase relative permeabilities using water-oil data
    if( ipWater >= 0 )
    {
      real64 const scaledWaterVolFrac = (phaseVolFraction[ipWater] - m_phaseMinVolumeFraction[dir][ipWater]) * volFracScaleInv;
      real64 const scaledOilVolFrac   = (phaseVolFraction[ipOil]   - m_phaseMinVolumeFraction[dir][ipOil])   * volFracScaleInv;

      using WOPT = RelativePermeabilityBase::WaterOilPairPhaseType;
      real64 const waterExponentInv = m_waterOilRelPermExponentInv[WOPT::WATER][dir];
      real64 const waterMaxValue = m_waterOilRelPermMaxValue[WOPT::WATER][dir];

      // water rel perm
      evaluateVanGenuchtenFunction( scaledWaterVolFrac,
                                    volFracScaleInv,
                                    waterExponentInv,
                                    waterMaxValue,
                                    phaseRelPerm[ipWater][dir],
                                    dPhaseRelPerm_dPhaseVolFrac[ipWater][ipWater][dir] );

      real64 const oilExponentInv_wo = m_waterOilRelPermExponentInv[WOPT::OIL][dir];
      real64 const oilMaxValue_wo = m_waterOilRelPermMaxValue[WOPT::OIL][dir];

      // oil rel perm
      evaluateVanGenuchtenFunction( scaledOilVolFrac,
                                    volFracScaleInv,
                                    oilExponentInv_wo,
                                    oilMaxValue_wo,
                                    oilRelPerm_wo,
                                    dOilRelPerm_wo_dOilVolFrac );

    }

    // 2) Gas and oil phase relative permeabilities using gas-oil data
    if( ipGas >= 0 )
    {
      real64 const scaledGasVolFrac = (phaseVolFraction[ipGas] - m_phaseMinVolumeFraction[dir][ipGas]) * volFracScaleInv;
      real64 const scaledOilVolFrac = (phaseVolFraction[ipOil] - m_phaseMinVolumeFraction[dir][ipOil]) * volFracScaleInv;

      using GOPT = RelativePermeabilityBase::GasOilPairPhaseType;
      real64 const gasExponentInv = m_gasOilRelPermExponentInv[GOPT::GAS][dir];
      real64 const gasMaxValue = m_gasOilRelPermMaxValue[GOPT::GAS][dir];
      // gas rel perm
      evaluateVanGenuchtenFunction( scaledGasVolFrac,
                                    volFracScaleInv,
                                    gasExponentInv,
                                    gasMaxValue,
                                    phaseRelPerm[ipGas][dir],
                                    dPhaseRelPerm_dPhaseVolFrac[ipGas][ipGas][dir] );

      real64 const oilExponentInv_go = m_gasOilRelPermExponentInv[GOPT::OIL][dir];
      real64 const oilMaxValue_go    = m_gasOilRelPermMaxValue[GOPT::OIL][dir];

      // oil rel perm
      evaluateVanGenuchtenFunction( scaledOilVolFrac,
                                    volFracScaleInv,
                                    oilExponentInv_go,
                                    oilMaxValue_go,
                                    oilRelPerm_go,
                                    dOilRelPerm_go_dOilVolFrac );
    }

    // 3) Compute the "three-phase" oil relperm

    // if no gas, use water-oil data
    if( ipGas < 0 )
    {
      phaseRelPerm[ipOil][dir] = oilRelPerm_wo;
      dPhaseRelPerm_dPhaseVolFrac[ipOil][ipOil][dir] = dOilRelPerm_wo_dOilVolFrac;
    }
    // if no water, use gas-oil data
    else if( ipWater < 0 )
    {
      phaseRelPerm[ipOil][dir] = oilRelPerm_go;
      dPhaseRelPerm_dPhaseVolFrac[ipOil][ipOil][dir] = dOilRelPerm_go_dOilVolFrac;
    }
    // if water and oil and gas can be present, use saturation-weighted interpolation
    else
    {
      real64 const shiftedWaterVolFrac = (phaseVolFraction[ipWater] - m_phaseMinVolumeFraction[dir][ipWater]);

      // TODO: change name of the class and add template to choose interpolation
      relpermInterpolators::Baker::compute( shiftedWaterVolFrac,
                                            phaseVolFraction[ipGas],
                                            m_phaseOrder,
                                            oilRelPerm_wo,
                                            dOilRelPerm_wo_dOilVolFrac,
                                            oilRelPerm_go,
                                            dOilRelPerm_go_dOilVolFrac,
                                            phaseRelPerm[ipOil][dir],
                                            dPhaseRelPerm_dPhaseVolFrac[ipOil][dir] );
    }

    // update trapped phase volume fraction
    if( ipWater >= 0 )
    {
      phaseTrappedVolFrac[ipWater] = LvArray::math::min( phaseVolFraction[ipWater], m_phaseMinVolumeFraction[dir][ipWater] );
    }
    if( ipGas >= 0 )
    {
      phaseTrappedVolFrac[ipGas] = LvArray::math::min( phaseVolFraction[ipGas], m_phaseMinVolumeFraction[dir][ipGas] );
    }
    if( ipOil >= 0 )
    {
      phaseTrappedVolFrac[ipOil] = LvArray::math::min( phaseVolFraction[ipOil], m_phaseMinVolumeFraction[dir][ipOil] );
    }
  }//end dir loop

}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
VanGenuchtenBakerRelativePermeabilityUpdate::
  evaluateVanGenuchtenFunction( real64 const & scaledVolFrac,
                                real64 const & dScaledVolFrac_dVolFrac,
                                real64 const & exponentInv,
                                real64 const & maxValue,
                                real64 & relPerm,
                                real64 & dRelPerm_dVolFrac )
{
  real64 const exponent = 1.0 / exponentInv;

  relPerm           = 0.0;
  dRelPerm_dVolFrac = 0.0;

  if( scaledVolFrac > 0.0 && scaledVolFrac < 1.0 )
  {
    // intermediate values
    real64 const a = pow( scaledVolFrac, exponent-1 );
    real64 const b = pow( 1 - a * scaledVolFrac, exponentInv-1 );
    real64 const c = ( 1 - b * ( 1 - a * scaledVolFrac ) );
    real64 const volFracSquared = scaledVolFrac * scaledVolFrac;
    real64 const dVolFracSquared_dVolFrac = 2 * dScaledVolFrac_dVolFrac * scaledVolFrac;

    relPerm = maxValue * volFracSquared * c;
    dRelPerm_dVolFrac = dVolFracSquared_dVolFrac * c + volFracSquared * dScaledVolFrac_dVolFrac * a * b;
    dRelPerm_dVolFrac *= maxValue;
  }
  else
  {
    relPerm = (scaledVolFrac <= 0.0) ? 0.0 : maxValue;
  }
}

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_VANGENUCHTENBAKERRELATIVEPERMEABILITY_HPP
