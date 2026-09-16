/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DamageCloggingPermeability.hpp
 */

#ifndef GEOS_CONSTITUTIVE_PERMEABILITY_DAMAGECLOGGINGPERMEABILITY_HPP_
#define GEOS_CONSTITUTIVE_PERMEABILITY_DAMAGECLOGGINGPERMEABILITY_HPP_

#include "constitutive/permeability/DamagePermeability.hpp"


namespace geos
{
namespace constitutive
{

class DamageCloggingPermeabilityUpdate : public DamagePermeabilityUpdate
{
public:

  DamageCloggingPermeabilityUpdate( arrayView3d< real64 > const & permeability,
                                    arrayView3d< real64 > const & dPerm_dPressure,
                                    arrayView1d< real64 > const & cloggedPoreFraction,
                                    real64 const & bulkPermeability,
                                    real64 const & damageDependenceConstant,
                                    real64 const & cloggingExponent,
                                    real64 const & minCloggingMultiplier )
    : DamagePermeabilityUpdate( permeability, dPerm_dPressure, bulkPermeability, damageDependenceConstant ),
    m_cloggedPoreFraction( cloggedPoreFraction ),
    m_cloggingExponent( cloggingExponent ),
    m_minCloggingMultiplier( minCloggingMultiplier )
  {}

  /**
   * @brief Damage-enhanced permeability reduced by mineral clogging, k = k_b exp(c d) max(f_min, (1-theta)^n)
   * @param k the element index
   * @param damage the element-averaged damage
   * @param cloggedPoreFraction fraction theta of the pre-precipitation pore space filled by precipitate
   */
  GEOS_HOST_DEVICE
  void updateDamageCloggingPermeability( localIndex const k,
                                         real64 const & damage,
                                         real64 const & cloggedPoreFraction ) const
  {
    real64 const openFraction = LvArray::math::max( 1.0 - cloggedPoreFraction, 0.0 );
    real64 const cloggingMultiplier = LvArray::math::max( m_minCloggingMultiplier,
                                                          pow( openFraction, m_cloggingExponent ) );
    real64 const matrixPermeability = m_bulkPermeability * LvArray::math::exp( m_damageDependenceConstant * damage ) * cloggingMultiplier;

    m_cloggedPoreFraction[k] = cloggedPoreFraction;
    for( localIndex dim=0; dim<3; ++dim )
    {
      m_permeability[k][0][dim] = matrixPermeability;
    }
  }

  GEOS_HOST_DEVICE
  real64 getCloggedPoreFraction( localIndex const k ) const { return m_cloggedPoreFraction[k]; }

private:

  /// Fraction of the pre-precipitation pore space filled by precipitate
  arrayView1d< real64 > m_cloggedPoreFraction;

  /// Exponent n of the clogging multiplier (1-theta)^n
  real64 m_cloggingExponent;

  /// Lower bound of the clogging multiplier
  real64 m_minCloggingMultiplier;

};


class DamageCloggingPermeability : public DamagePermeability
{
public:

  DamageCloggingPermeability( string const & name, Group * const parent );

  static string catalogName() { return "DamageCloggingPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = DamageCloggingPermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dPressure,
                          m_cloggedPoreFraction,
                          m_bulkPermeability,
                          m_damageDependenceConstant,
                          m_cloggingExponent,
                          m_minCloggingMultiplier );
  }

  struct viewKeyStruct : public DamagePermeability::viewKeyStruct
  {
    static constexpr char const * cloggingExponentString() { return "cloggingExponent"; }
    static constexpr char const * minCloggingMultiplierString() { return "minCloggingMultiplier"; }
    static constexpr char const * cloggedPoreFractionString() { return "cloggedPoreFraction"; }
  };

private:

  /// Fraction of the pre-precipitation pore space filled by precipitate
  array1d< real64 > m_cloggedPoreFraction;

  /// Exponent n of the clogging multiplier (1-theta)^n
  real64 m_cloggingExponent;

  /// Lower bound of the clogging multiplier
  real64 m_minCloggingMultiplier;

};

}/* namespace constitutive */

} /* namespace geos */


#endif //GEOS_CONSTITUTIVE_PERMEABILITY_DAMAGECLOGGINGPERMEABILITY_HPP_
