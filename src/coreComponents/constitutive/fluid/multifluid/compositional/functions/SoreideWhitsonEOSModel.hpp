/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SoreideWhitsonEOSModel.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_SOREIDEWHITSONEOSMODEL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_SOREIDEWHITSONEOSMODEL_HPP_

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ComponentProperties.hpp"
#include "CubicEOSPhaseModel.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

enum class SoreideWhitsonPhaseType : integer
{
  Vapour,
  Aqueous
};

enum class EquationOfStateType : integer
{
  PengRobinson,
  SoaveRedlichKwong
};

namespace
{

template< EquationOfStateType EOS_TYPE >
struct CubicModel {};

template<>
struct CubicModel< EquationOfStateType::PengRobinson >
{
  using type = PengRobinsonEOS;
};

template<>
struct CubicModel< EquationOfStateType::SoaveRedlichKwong >
{
  using type = SoaveRedlichKwongEOS;
};

}

template< SoreideWhitsonPhaseType PHASE_TYPE, EquationOfStateType EOS_TYPE >
struct SoreideWhitsonEOSModel
{
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  /**
   * @brief Generate a catalog name
   */
  static constexpr char const * catalogName(){ return "SoreideWhitson"; }

  /**
   * @brief Calculate the pure coefficients
   * @details Computes the pure coefficients
   * @param[in] ic Component index
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] componentProperties The compositional component properties
   * @param[out] aCoefficient pure coefficient (A)
   * @param[out] bCoefficient pure coefficient (B)
   */
  GEOS_HOST_DEVICE
  static void
  computePureCoefficients( integer const ic,
                           real64 const & pressure,
                           real64 const & temperature,
                           ComponentProperties::KernelWrapper const & componentProperties,
                           real64 & aCoefficient,
                           real64 & bCoefficient );

  /**
   * @brief Calculate the pure coefficients derivatives
   * @details Computes the pure coefficients derivatives
   * @param[in] ic Component index
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] componentProperties The compositional component properties
   * @param[out] aCoefficient pure coefficient (A)
   * @param[out] bCoefficient pure coefficient (B)
   * @param[out] daCoefficient_dp pure coefficient (A) derivative w.r.t. pressure
   * @param[out] dbCoefficient_dp pure coefficient (B) derivative w.r.t. pressure
   * @param[out] daCoefficient_dt pure coefficient (A) derivative w.r.t. temperature
   * @param[out] dbCoefficient_dt pure coefficient (B) derivative w.r.t. temperature
   */
  GEOS_HOST_DEVICE
  static void
  computePureCoefficients( integer const ic,
                           real64 const & pressure,
                           real64 const & temperature,
                           ComponentProperties::KernelWrapper const & componentProperties,
                           real64 & aCoefficient,
                           real64 & bCoefficient,
                           real64 & daCoefficient_dp,
                           real64 & dbCoefficient_dp,
                           real64 & daCoefficient_dt,
                           real64 & dbCoefficient_dt );
};

namespace
{
template< SoreideWhitsonPhaseType PHASE_TYPE, EquationOfStateType EOS_TYPE >
struct PureCoefficientCalculator {};

template< EquationOfStateType EOS_TYPE >
struct PureCoefficientCalculator< SoreideWhitsonPhaseType::Aqueous, EOS_TYPE >
{
  GEOS_HOST_DEVICE
  static void
  calculate( integer const ic,
             real64 const & pressure,
             real64 const & temperature,
             ComponentProperties::KernelWrapper const & componentProperties,
             real64 & aCoefficient,
             real64 & bCoefficient,
             real64 & daCoefficient_dp,
             real64 & dbCoefficient_dp,
             real64 & daCoefficient_dt,
             real64 & dbCoefficient_dt )
  {
    using CubicEOS = typename CubicModel< EOS_TYPE >::type;
    arraySlice1d< real64 const > const & criticalPressure = componentProperties.m_componentCriticalPressure;
    arraySlice1d< real64 const > const & criticalTemperature = componentProperties.m_componentCriticalTemperature;

    real64 const pr = pressure / criticalPressure[ic];
    real64 const tr = temperature / criticalTemperature[ic];

    real64 const trToMinus2 = 1.0/(tr*tr);
    real64 const trToMinus3 = trToMinus2/tr;

    real64 constexpr salinity = 0.1;
    real64 const csw = salinity < MultiFluidConstants::minForSpeciesPresence ? 0.0 : LvArray::math::exp( 1.1 * LvArray::math::log( salinity ) );

    real64 const sqrtAlpha = 1.0 - 0.4530*(1.0 - tr*(1.0 - 0.0103*csw)) + 0.0034*(trToMinus3 - 1.0);
    real64 const dSqrtAlpha_dtr = 0.4530*(1.0 - 0.0103*csw) - 0.0102*trToMinus3/tr;
    real64 const alpha = sqrtAlpha * sqrtAlpha;
    real64 const dAlpha_dt = 2.0 * sqrtAlpha * dSqrtAlpha_dtr / criticalTemperature[ic];

    aCoefficient = CubicEOS::omegaA * pr * trToMinus2 * alpha;
    bCoefficient = CubicEOS::omegaB * pr / tr;

    daCoefficient_dp = aCoefficient / pressure;
    dbCoefficient_dp = bCoefficient / pressure;

    daCoefficient_dt = CubicEOS::omegaA * pr * (-2.0*trToMinus3 * alpha / criticalTemperature[ic] + trToMinus2 * dAlpha_dt);
    dbCoefficient_dt = -bCoefficient / temperature;
  }
};

template< EquationOfStateType EOS_TYPE >
struct PureCoefficientCalculator< SoreideWhitsonPhaseType::Vapour, EOS_TYPE >
{
  GEOS_HOST_DEVICE
  static void calculate( integer const ic,
                         real64 const & pressure,
                         real64 const & temperature,
                         ComponentProperties::KernelWrapper const & componentProperties,
                         real64 & aCoefficient,
                         real64 & bCoefficient,
                         real64 & daCoefficient_dp,
                         real64 & dbCoefficient_dp,
                         real64 & daCoefficient_dt,
                         real64 & dbCoefficient_dt )
  {
    using CubicEOS = CubicEOSPhaseModel< typename CubicModel< EOS_TYPE >::type >;
    CubicEOS::computePureCoefficients( ic,
                                       pressure,
                                       temperature,
                                       componentProperties,
                                       aCoefficient,
                                       bCoefficient,
                                       daCoefficient_dp,
                                       dbCoefficient_dp,
                                       daCoefficient_dt,
                                       dbCoefficient_dt );
  }
};
}

template< SoreideWhitsonPhaseType PHASE_TYPE, EquationOfStateType EOS_TYPE >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSModel< PHASE_TYPE, EOS_TYPE >::
computePureCoefficients( integer const ic,
                         real64 const & pressure,
                         real64 const & temperature,
                         ComponentProperties::KernelWrapper const & componentProperties,
                         real64 & aCoefficient,
                         real64 & bCoefficient )
{
  real64 daCoefficient_dp = 0.0;
  real64 dbCoefficient_dp = 0.0;
  real64 daCoefficient_dt = 0.0;
  real64 dbCoefficient_dt = 0.0;
  computePureCoefficients( ic,
                           pressure,
                           temperature,
                           componentProperties,
                           aCoefficient,
                           bCoefficient,
                           daCoefficient_dp,
                           dbCoefficient_dp,
                           daCoefficient_dt,
                           dbCoefficient_dt );
}

template< SoreideWhitsonPhaseType PHASE_TYPE, EquationOfStateType EOS_TYPE >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSModel< PHASE_TYPE, EOS_TYPE >::
computePureCoefficients( integer const ic,
                         real64 const & pressure,
                         real64 const & temperature,
                         ComponentProperties::KernelWrapper const & componentProperties,
                         real64 & aCoefficient,
                         real64 & bCoefficient,
                         real64 & daCoefficient_dp,
                         real64 & dbCoefficient_dp,
                         real64 & daCoefficient_dt,
                         real64 & dbCoefficient_dt )
{
  PureCoefficientCalculator< PHASE_TYPE, EOS_TYPE >::calculate(
    ic,
    pressure,
    temperature,
    componentProperties,
    aCoefficient,
    bCoefficient,
    daCoefficient_dp,
    dbCoefficient_dp,
    daCoefficient_dt,
    dbCoefficient_dt );
}

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_SOREIDEWHITSONEOSMODEL_HPP_
