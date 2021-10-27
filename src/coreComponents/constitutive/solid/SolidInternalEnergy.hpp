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
 * @file SolidInternalEnergy.hpp
 */


#ifndef GEOSX_CONSTITUTIVE_SOLID_SOLIDINTERNALENERGY_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_SOLIDINTERNALENERGY_HPP_

#include "constitutive/ConstitutiveBase.hpp"

namespace geosx
{

namespace constitutive
{


class SolidInternalEnergyUpdates
{
public:

  SolidInternalEnergyUpdates( arrayView2d< real64 > const & internalEnergy,
                              arrayView2d< real64 > const & dInternalEnergy_dTemperature,
                              real64 const & heatCapacity,
                              real64 const & referenceTemperature,
                              real64 const & referenceInternalEnergy ):
    m_internalEnergy( internalEnergy ),
    m_dInternalEnergy_dTemperature( dInternalEnergy_dTemperature ),
    m_heatCapacity( heatCapacity ),
    m_referenceTemperature( referenceTemperature ),
    m_referenceInternalEnergy( referenceInternalEnergy )
  {}

  GEOSX_DEVICE
  void update( localIndex const k,
               real64 const & temperature ) const
  {
    compute( temperature,
             m_internalEnergy[k][0],
             m_dInternalEnergy_dTemperature[k][0] );
  }

  GEOSX_DEVICE
  void compute( real64 const & temperature,
                real64 & internalEnergy,
                real64 & dInternalEnergy_dTemperature ) const
  {
    internalEnergy = m_referenceInternalEnergy + m_heatCapacity * ( temperature - m_referenceTemperature );
    dInternalEnergy_dTemperature =  m_heatCapacity;
  }

private:

  /// Solid internal energy
  arrayView2d< real64 > m_internalEnergy;

  /// Derivative of the solid internal Energy w.r.t. to the temperature
  arrayView2d< real64 > m_dInternalEnergy_dTemperature;

  /// solid specific heat capacity
  real64 m_heatCapacity;

  /// Reference temperature
  real64 m_referenceTemperature;

  /// Reference internal energy
  real64 m_referenceInternalEnergy;
};

class SolidInternalEnergy : public ConstitutiveBase
{
public:

  SolidInternalEnergy( string const & name, Group * const parent );

  static string catalogName() { return "SolidInternalEnergy"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * internalEnergyString() { return "internalEnergy"; }
    static constexpr char const * oldInternalEnergyString() { return "oldInternalEnergy"; }
    static constexpr char const * dInternalEnergy_dTemperatureString() { return "dInternalEnergy_dTemperature"; }
    static constexpr char const * heatCapacityString() { return "heatCapacity"; }
    static constexpr char const * referenceTemperatureString() { return "referenceTemperature"; }
    static constexpr char const * referenceInternalEnergyString() { return "referenceInternalEnergy"; }
  } viewKeys;

  using KernelWrapper = SolidInternalEnergyUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelUpdates()
  {
    return KernelWrapper( m_internalEnergy,
                          m_dInternalEnergy_dTemperature,
                          m_heatCapacity,
                          m_referenceTemperature,
                          m_referenceInternalEnergy );
  }


  /**
   * @brief Const accessor for internalEnergy.
   * @return Accessor
   */
  arrayView2d< real64 const > const  getInternalEnergy() const { return m_internalEnergy; }

  /**
   * @brief Const/non-mutable accessor for oldInternalEnergy.
   * @return Accessor
   */
  arrayView2d< real64 const > const  getOldInternalEnergy() const { return m_oldInternalEnergy; }

  /**
   * @brief Const/non-mutable accessor for oldInternalEnergy.
   * @return Accessor
   */
  arrayView2d< real64 const > const  getDinternalEnergy_dTemperature() const { return m_dInternalEnergy_dTemperature; }

  /// Save state data in preparation for next timestep
  virtual void saveConvergedState() const override;

private:

  /// Solid internal energy
  array2d< real64 > m_internalEnergy;

  /// Old solid internal energy
  array2d< real64 > m_oldInternalEnergy;

  /// Derivative of the solid internal Energy w.r.t. to the temperature
  array2d< real64 > m_dInternalEnergy_dTemperature;

  /// solid specific heat capacity
  real64 m_heatCapacity;

  /// Reference temperature
  real64 m_referenceTemperature;

  /// Reference internal energy
  real64 m_referenceInternalEnergy;

};

} /* namespace constitutive */

} /* namespace geosx */

#endif
