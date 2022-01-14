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
 * @file ConstantThermalConductivity.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_THERMALCONDUCTIVITY_CONSTANTTHERMALCONDUCTIVITY_HPP_
#define GEOSX_CONSTITUTIVE_THERMALCONDUCTIVITY_CONSTANTTHERMALCONDUCTIVITY_HPP_

#include "constitutive/thermalConductivity/ThermalConductivityBase.hpp"


namespace geosx
{
namespace constitutive
{

/**
 * @brief The update class for constant thermal conductivity (does not do anything)
 */
class ConstantThermalConductivityUpdate : public ThermalConductivityBaseUpdate
{
public:

  /**
   * @brief Constructor for the class performing the thermal conductivity updates
   * @param effectiveConductivity the array of cell-wise effective conductivities in the subregion
   * @param dEffectiveConductivity_dPhaseVolFrac the array of cell-wise derivatives of effective conductivities wrt phase vol fractions in
   * the subregion
   */
  ConstantThermalConductivityUpdate( arrayView2d< real64 > const & effectiveConductivity,
                                     arrayView3d< real64 > const & dEffectiveConductivity_dPhaseVolFrac )
    : ThermalConductivityBaseUpdate( effectiveConductivity, dEffectiveConductivity_dPhaseVolFrac )
  {}

  GEOSX_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const & laggedPorosity,
                       arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const override
  { GEOSX_UNUSED_VAR( k, q, laggedPorosity, phaseVolFraction ); }

};

/**
 * @brief The class for constant thermal conductivity
 */
class ConstantThermalConductivity : public ThermalConductivityBase
{
public:

  /**
   * @brief Constructor for the class storing constant conductivity
   * @param[in] name the name of the class
   * @param[in] parent pointer to the parent Group
   */
  ConstantThermalConductivity( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "ConstantThermalConductivity"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ConstantThermalConductivityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_effectiveConductivity,
                          m_dEffectiveConductivity_dPhaseVolFrac );
  }

  struct viewKeyStruct : public ThermalConductivityBase::viewKeyStruct
  {
    static constexpr char const * defaultThermalConductivityString() { return "defaultThermalConductivity"; }
  } viewKeys;

protected:

  virtual void postProcessInput() override;

private:

  /// default thermal conductivity in the subRegion
  real64 m_defaultThermalConductivity;

};

} // namespace constitutive

} // namespace geosx


#endif //GEOSX_CONSTITUTIVE_THERMALCONDUCTIVITY_CONSTANTTHERMALCONDUCTIVITY_HPP_
