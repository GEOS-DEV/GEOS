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
 * @file VolumeWeightedThermalConductivity.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_THERMALCONDUCTIVITY_VOLUMEWEIGHTEDTHERMALCONDUCTIVITY_HPP_
#define GEOSX_CONSTITUTIVE_THERMALCONDUCTIVITY_VOLUMEWEIGHTEDTHERMALCONDUCTIVITY_HPP_

#include "constitutive/thermalConductivity/ThermalConductivityBase.hpp"


namespace geosx
{
namespace constitutive
{

/**
 * @brief The update class for volume-weighted thermal conductivity
 */
class VolumeWeightedThermalConductivityUpdate : public ThermalConductivityBaseUpdate
{
public:

  /**
   * @brief Constructor for the class performing the thermal conductivity updates
   * @param effectiveConductivity the array of cell-wise effective conductivities in the subregion
   * @param dEffectiveConductivity_dPhaseVolFrac the array of cell-wise derivatives of effective conductivities wrt phase vol fractions in
   * the subregion
   * @param rockThermalConductivity the array of cell-wise rock thermal conductivities
   * @param phaseThermalConductivity the array of fluid phase thermal conductivities
   */
  VolumeWeightedThermalConductivityUpdate( arrayView3d< real64 > const & effectiveConductivity,
                                           arrayView4d< real64 > const & dEffectiveConductivity_dPhaseVolFrac,
                                           arrayView3d< real64 const > const & rockThermalConductivity,
                                           arrayView1d< real64 const > const & phaseThermalConductivity )
    : ThermalConductivityBaseUpdate( effectiveConductivity, dEffectiveConductivity_dPhaseVolFrac ),
    m_rockThermalConductivity( rockThermalConductivity ),
    m_phaseThermalConductivity( phaseThermalConductivity )
  {}

  /**
   * @brief Cell-wise compute function for the volume-weighted effective conductivity
   * @param[in] laggedPorosity the cell lagged porosity
   * @param[in] phaseVolFrac the cell phase volume fraction
   * @param[in] rockThermalConductivity the cell rock thermal conductivity
   * @param[out] effectiveConductivity the cell effective conductivity
   * @param[out] dEffectiveConductivity_dPhaseVolFrac the cell derivative of effective conductivity wrt phase volume fractions
   */
  GEOSX_HOST_DEVICE
  void compute( real64 const & laggedPorosity,
                arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFrac,
                arraySlice1d< real64 const > const rockThermalConductivity,
                arraySlice1d< real64 > const effectiveConductivity,
                arraySlice2d< real64 > const dEffectiveConductivity_dPhaseVolFrac ) const
  {
    LvArray::forValuesInSlice( effectiveConductivity, []( real64 & val ){ val = 0.0; } );
    for( integer ip = 0; ip < numPhases(); ++ip )
    {
      real64 const phaseVolFracTimesConductivity = phaseVolFrac[ip] * m_phaseThermalConductivity[ip];
      real64 const porosityTimesConductivity = m_phaseThermalConductivity[ip] * laggedPorosity;
      for( integer dir = 0; dir < 3; ++dir )
      {
        effectiveConductivity[dir] += phaseVolFracTimesConductivity;
        dEffectiveConductivity_dPhaseVolFrac[dir][ip] = porosityTimesConductivity;
      }
    }
    for( integer dir = 0; dir < 3; ++dir )
    {
      effectiveConductivity[dir] *= laggedPorosity;
      effectiveConductivity[dir] += ( 1.0 - laggedPorosity ) * rockThermalConductivity[dir];
    }
  }

  GEOSX_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const & laggedPorosity,
                       arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFrac ) const override
  {
    compute( laggedPorosity,
             phaseVolFrac,
             m_rockThermalConductivity[k][q],
             m_effectiveConductivity[k][q],
             m_dEffectiveConductivity_dPhaseVolFrac[k][q] );
  }

private:

  /// View on the cell-wise rock thermal conductivity
  arrayView3d< real64 const > m_rockThermalConductivity;

  /// View on the phase thermal conductivity
  arrayView1d< real64 const > m_phaseThermalConductivity;

};

/**
 * @brief The class for volume-weighted thermal conductivity
 */
class VolumeWeightedThermalConductivity : public ThermalConductivityBase
{
public:

  /**
   * @brief Constructor for the class implementing volume-weighted conductivity
   * @param[in] name the name of the class
   * @param[in] parent pointer to the parent Group
   */
  VolumeWeightedThermalConductivity( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "VolumeWeightedThermalConductivity"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = VolumeWeightedThermalConductivityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_effectiveConductivity,
                          m_dEffectiveConductivity_dPhaseVolFrac,
                          m_rockThermalConductivity,
                          m_phaseThermalConductivity );
  }

  struct viewKeyStruct : public ThermalConductivityBase::viewKeyStruct
  {
    static constexpr char const * rockThermalConductivityComponentsString() { return "rockThermalConductivityComponents"; }
    static constexpr char const * phaseThermalConductivityString() { return "phaseThermalConductivity"; }
  } viewKeys;

protected:

  virtual void postProcessInput() override;

private:

  /// default rock thermal conductivity in the subRegion
  R1Tensor m_rockThermalConductivityComponents;

  /// cell-wise rock thermal conductivity in the subRegion
  array3d< real64 > m_rockThermalConductivity;

  /// fluid phase thermal conductivity
  array1d< real64 > m_phaseThermalConductivity;

};

} // namespace constitutive

} // namespace geosx


#endif //GEOSX_CONSTITUTIVE_THERMALCONDUCTIVITY_VOLUMEWEIGHTEDTHERMALCONDUCTIVITY_HPP_
