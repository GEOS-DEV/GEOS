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
 * @file BrooksCoreyRelativePermeability.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_BROOKSCOREYRELATIVEPERMEABILITY_HPP
#define GEOSX_CONSTITUTIVE_BROOKSCOREYRELATIVEPERMEABILITY_HPP

#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"

namespace geosx
{
namespace constitutive
{

class BrooksCoreyRelativePermeabilityUpdate final : public RelativePermeabilityBaseUpdate
{
public:

  BrooksCoreyRelativePermeabilityUpdate( arrayView1d< real64 const > const & phaseMinVolumeFraction,
                                         arrayView1d< real64 const > const & phaseRelPermExponent,
                                         arrayView1d< real64 const > const & phaseRelPermMaxValue,
                                         real64 const volFracScale,
                                         arrayView1d< integer const > const & phaseTypes,
                                         arrayView1d< integer const > const & phaseOrder,
                                         arrayView3d< real64, relperm::USD_RELPERM > const & phaseRelPerm,
                                         arrayView4d< real64, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac )
    : RelativePermeabilityBaseUpdate( phaseTypes,
                                      phaseOrder,
                                      phaseRelPerm,
                                      dPhaseRelPerm_dPhaseVolFrac ),
    m_phaseMinVolumeFraction( phaseMinVolumeFraction ),
    m_phaseRelPermExponent( phaseRelPermExponent ),
    m_phaseRelPermMaxValue( phaseRelPermMaxValue ),
    m_volFracScale( volFracScale )
  {}

  GEOSX_HOST_DEVICE
  virtual void compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                        arraySlice1d< real64, relperm::USD_RELPERM - 2 > const & phaseRelPerm,
                        arraySlice2d< real64, relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const override;

  GEOSX_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const override
  {
    compute( phaseVolFraction,
             m_phaseRelPerm[k][q],
             m_dPhaseRelPerm_dPhaseVolFrac[k][q] );
  }

private:

  arrayView1d< real64 const > m_phaseMinVolumeFraction;
  arrayView1d< real64 const > m_phaseRelPermExponent;
  arrayView1d< real64 const > m_phaseRelPermMaxValue;
  real64 m_volFracScale;
};

class BrooksCoreyRelativePermeability : public RelativePermeabilityBase
{
public:

  BrooksCoreyRelativePermeability( string const & name, dataRepository::Group * const parent );

//START_SPHINX_INCLUDE_00
  static string catalogName() { return "BrooksCoreyRelativePermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = BrooksCoreyRelativePermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

//START_SPHINX_INCLUDE_01
  struct viewKeyStruct : RelativePermeabilityBase::viewKeyStruct
  {
    static constexpr char const * phaseMinVolumeFractionString() { return "phaseMinVolumeFraction"; }
    static constexpr char const * phaseRelPermExponentString() { return "phaseRelPermExponent"; }
    static constexpr char const * phaseRelPermMaxValueString() { return "phaseRelPermMaxValue"; }
    static constexpr char const * volFracScaleString() { return "volFracScale"; }
  } vieKeysBrooksCoreyRelativePermeability;
//END_SPHINX_INCLUDE_01

protected:

  virtual void postProcessInput() override;

//START_SPHINX_INCLUDE_02
  array1d< real64 > m_phaseMinVolumeFraction;
  array1d< real64 > m_phaseRelPermExponent;
  array1d< real64 > m_phaseRelPermMaxValue;

  real64 m_volFracScale;
//END_SPHINX_INCLUDE_02
};

GEOSX_HOST_DEVICE
inline void
BrooksCoreyRelativePermeabilityUpdate::
  compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
           arraySlice1d< real64, relperm::USD_RELPERM - 2 > const & phaseRelPerm,
           arraySlice2d< real64, relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const
{
  LvArray::forValuesInSlice( dPhaseRelPerm_dPhaseVolFrac, []( real64 & val ){ val = 0.0; } );

  real64 const satScaleInv = 1.0 / m_volFracScale;

  for( localIndex ip = 0; ip < numPhases(); ++ip )
  {
    real64 const satScaled = (phaseVolFraction[ip] - m_phaseMinVolumeFraction[ip]) * satScaleInv;
    real64 const exponent  = m_phaseRelPermExponent[ip];
    real64 const scale     = m_phaseRelPermMaxValue[ip];

    if( satScaled > 0.0 && satScaled < 1.0 )
    {
      // intermediate value
      real64 const v = scale * pow( satScaled, exponent - 1.0 );

      phaseRelPerm[ip] = v * satScaled;
      dPhaseRelPerm_dPhaseVolFrac[ip][ip] = v * exponent * satScaleInv;
    }
    else
    {
      phaseRelPerm[ip] = (satScaled <= 0.0) ? 0.0 : scale;
    }
  }
}

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_BROOKSCOREYRELATIVEPERMEABILITY_HPP
