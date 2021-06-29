/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
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

  BrooksCoreyRelativePermeabilityUpdate( arrayView2d< real64 const > const & phaseMinVolumeFraction,
                                         arrayView2d< real64 const > const & phaseRelPermExponent,
                                         arrayView2d< real64 const > const & phaseRelPermMaxValue,
                                         arrayView1d< real64 const > const & volFracScale,
                                         arrayView1d< integer const > const & phaseTypes,
                                         arrayView1d< integer const > const & phaseOrder,
                                         arrayView3d< real64 > const & phaseRelPerm,
                                         arrayView4d< real64 > const & dPhaseRelPerm_dPhaseVolFrac )
    : RelativePermeabilityBaseUpdate( phaseTypes,
                                      phaseOrder,
                                      phaseRelPerm,
                                      dPhaseRelPerm_dPhaseVolFrac ),
    m_phaseMinVolumeFraction( phaseMinVolumeFraction ),
    m_phaseRelPermExponent( phaseRelPermExponent ),
    m_phaseRelPermMaxValue( phaseRelPermMaxValue ),
    m_volFracScale( volFracScale )
  {}

  /// Default copy constructor
  BrooksCoreyRelativePermeabilityUpdate( BrooksCoreyRelativePermeabilityUpdate const & ) = default;

  /// Default move constructor
  BrooksCoreyRelativePermeabilityUpdate( BrooksCoreyRelativePermeabilityUpdate && ) = default;

  /// Deleted copy assignment operator
  BrooksCoreyRelativePermeabilityUpdate & operator=( BrooksCoreyRelativePermeabilityUpdate const & ) = delete;

  /// Deleted move assignment operator
  BrooksCoreyRelativePermeabilityUpdate & operator=( BrooksCoreyRelativePermeabilityUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void compute( arraySlice1d< real64 const > const & phaseVolFraction,
                arraySlice1d< real64 > const & phaseRelPerm,
                arraySlice2d< real64 > const & dPhaseRelPerm_dPhaseVolFrac,
                arraySlice1d< real64 const > const & phaseMinVolumeFraction,
                arraySlice1d< real64 const > const & phaseRelPermExponent,
                arraySlice1d< real64 const > const & phaseRelPermMaxValue,
                real64 const & volFracScale ) const;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void update( localIndex const k,
                       localIndex const q,
                       arraySlice1d< real64 const > const & phaseVolFraction ) const override
  {
    compute( phaseVolFraction,
             m_phaseRelPerm[k][q],
             m_dPhaseRelPerm_dPhaseVolFrac[k][q],
             m_phaseMinVolumeFraction[k],
             m_phaseRelPermExponent[k],
             m_phaseRelPermMaxValue[k],
             m_volFracScale[k] );
  }

private:

  arrayView2d< real64 const > m_phaseMinVolumeFraction;
  arrayView2d< real64 const > m_phaseRelPermExponent;
  arrayView2d< real64 const > m_phaseRelPermMaxValue;
  arrayView1d< real64 const > m_volFracScale;
};

class BrooksCoreyRelativePermeability : public RelativePermeabilityBase
{
public:

  BrooksCoreyRelativePermeability( string const & name, dataRepository::Group * const parent );

  virtual ~BrooksCoreyRelativePermeability() override;

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
    static constexpr char const * defaultPhaseMinVolumeFractionString() { return "defaultPhaseMinVolumeFraction"; }
    static constexpr char const * defaultPhaseRelPermExponentString() { return "defaultPhaseRelPermExponent"; }
    static constexpr char const * defaultPhaseRelPermMaxValueString() { return "defaultPhaseRelPermMaxValue"; }

    static constexpr char const * phaseMinVolumeFractionString() { return "phaseMinVolumeFraction"; }
    static constexpr char const * phaseRelPermExponentString() { return "phaseRelPermExponent"; }
    static constexpr char const * phaseRelPermMaxValueString() { return "phaseRelPermMaxValue"; }
    static constexpr char const * volFracScaleString() { return "volFracScale"; }

    dataRepository::ViewKey defaultPhaseMinVolumeFraction = { defaultPhaseMinVolumeFractionString() };
    dataRepository::ViewKey defaultPhaseRelPermExponent   = { defaultPhaseRelPermExponentString() };
    dataRepository::ViewKey defaultPhaseRelPermMaxValue   = { defaultPhaseRelPermMaxValueString() };
    dataRepository::ViewKey phaseMinVolumeFraction = { phaseMinVolumeFractionString() };
    dataRepository::ViewKey phaseRelPermExponent   = { phaseRelPermExponentString() };
    dataRepository::ViewKey phaseRelPermMaxValue   = { phaseRelPermMaxValueString() };
  } vieKeysBrooksCoreyRelativePermeability;
//END_SPHINX_INCLUDE_01

protected:

  virtual void postProcessInput() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  virtual void resizeFields( localIndex const size, localIndex const numPts ) override;

//START_SPHINX_INCLUDE_02
  array1d< real64 > m_defaultPhaseMinVolumeFraction;
  array1d< real64 > m_defaultPhaseRelPermExponent;
  array1d< real64 > m_defaultPhaseRelPermMaxValue;

  // TODO: add the quadrature point dimension
  array2d< real64 > m_phaseMinVolumeFraction;
  array2d< real64 > m_phaseRelPermExponent;
  array2d< real64 > m_phaseRelPermMaxValue;

  array1d< real64 > m_volFracScale;
//END_SPHINX_INCLUDE_02
};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
BrooksCoreyRelativePermeabilityUpdate::
  compute( arraySlice1d< real64 const > const & phaseVolFraction,
           arraySlice1d< real64 > const & phaseRelPerm,
           arraySlice2d< real64 > const & dPhaseRelPerm_dPhaseVolFrac,
           arraySlice1d< real64 const > const & phaseMinVolumeFraction,
           arraySlice1d< real64 const > const & phaseRelPermExponent,
           arraySlice1d< real64 const > const & phaseRelPermMaxValue,
           real64 const & volFracScale ) const
{
  for( real64 & val : dPhaseRelPerm_dPhaseVolFrac )
  {
    val = 0.0;
  }

  real64 const satScaleInv = 1.0 / volFracScale;

  for( localIndex ip = 0; ip < numPhases(); ++ip )
  {
    real64 const satScaled = (phaseVolFraction[ip] - phaseMinVolumeFraction[ip]) * satScaleInv;
    real64 const exponent  = phaseRelPermExponent[ip];
    real64 const scale     = phaseRelPermMaxValue[ip];

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
