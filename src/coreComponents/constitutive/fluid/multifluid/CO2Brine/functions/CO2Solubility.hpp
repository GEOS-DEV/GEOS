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
 * @file CO2Solubility.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_CO2SOLUBILITY_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_CO2SOLUBILITY_HPP_

#include "FlashModelBase.hpp"

#include "constitutive/fluid/multifluid/CO2Brine/functions/PVTFunctionHelpers.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{

namespace constitutive
{

namespace PVTProps
{

constexpr real64 minForDivision = 1e-10;

class CO2SolubilityUpdate final : public FlashModelBaseUpdate
{
public:

  using PhaseProp = MultiFluidVar< real64, 3 >;
  using PhaseComp = MultiFluidVar< real64, 4 >;

  CO2SolubilityUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                       TableFunction const & CO2SolubilityTable,
                       integer const CO2Index,
                       integer const waterIndex,
                       integer const phaseGasIndex,
                       integer const phaseLiquidIndex )
    : FlashModelBaseUpdate( componentMolarWeight ),
    m_CO2SolubilityTable( CO2SolubilityTable.createKernelWrapper() ),
    m_CO2Index( CO2Index ),
    m_waterIndex( waterIndex ),
    m_phaseGasIndex( phaseGasIndex ),
    m_phaseLiquidIndex( phaseLiquidIndex )
  {}

  template< int USD1 >
  GEOS_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & compFraction,
                PhaseProp::SliceType const phaseFraction,
                PhaseComp::SliceType const phaseCompFraction ) const;

  virtual void move( LvArray::MemorySpace const space, bool const touch ) override
  {
    FlashModelBaseUpdate::move( space, touch );
    m_CO2SolubilityTable.move( space, touch );
  }

protected:

  /// Table with CO2 solubility tabulated as a function (P,T)
  TableFunction::KernelWrapper m_CO2SolubilityTable;

  /// Index of the CO2 phase
  integer m_CO2Index;

  /// Index of the water phase
  integer m_waterIndex;

  /// Index of the gas phase
  integer m_phaseGasIndex;

  /// Index of the liquid phase
  integer m_phaseLiquidIndex;

};

class CO2Solubility : public FlashModelBase
{
public:

  CO2Solubility( string const & name,
                 string_array const & inputParams,
                 string_array const & phaseNames,
                 string_array const & componentNames,
                 array1d< real64 > const & componentMolarWeight,
                 bool const printTable );

  static string catalogName() { return "CO2Solubility"; }

  virtual string getCatalogName() const final { return catalogName(); }

  /**
   * @copydoc FlashModelBase::checkTablesParameters( real64 pressure, real64 temperature )
   */
  void checkTablesParameters( real64 pressure, real64 temperature ) const override final;

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = CO2SolubilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

private:

  /// Table to compute solubility as a function of pressure and temperature
  TableFunction const * m_CO2SolubilityTable;

  /// Index of the CO2 component
  integer m_CO2Index;

  /// Index of the water component
  integer m_waterIndex;

  /// Index of the gas phase
  integer m_phaseGasIndex;

  /// Index of the liquid phase
  integer m_phaseLiquidIndex;
};

template< int USD1 >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
CO2SolubilityUpdate::compute( real64 const & pressure,
                              real64 const & temperature,
                              arraySlice1d< real64 const, USD1 > const & compFraction,
                              PhaseProp::SliceType const phaseFraction,
                              PhaseComp::SliceType const phaseCompFraction ) const
{
  using Deriv = multifluid::DerivativeOffset;

  // solubility mol/kg(water)  X = Csat/W
  real64 const input[2] = { pressure, temperature };
  real64 solubilityDeriv[2]{};
  real64 solubility = m_CO2SolubilityTable.compute( input, solubilityDeriv );

  solubility *= m_componentMolarWeight[m_waterIndex];
  for( integer ic = 0; ic < 2; ++ic )
  {
    solubilityDeriv[ic] *= m_componentMolarWeight[m_waterIndex];
  }

  // Y = C/W = z/(1-z)
  real64 Y = 0.0;
  real64 dY_dCompFrac[2]{};

  if( compFraction[m_CO2Index] > 1.0 - minForDivision )
  {
    Y = compFraction[m_CO2Index] / minForDivision;
    dY_dCompFrac[m_CO2Index] = 1.0 / minForDivision;
    dY_dCompFrac[m_waterIndex] = 0.0;
  }
  else
  {
    real64 const oneMinusCompFracInv = 1.0 / (1.0 - compFraction[m_CO2Index]);
    Y = compFraction[m_CO2Index] * oneMinusCompFracInv;
    dY_dCompFrac[m_CO2Index] = oneMinusCompFracInv * oneMinusCompFracInv;
    dY_dCompFrac[m_waterIndex] = 0.0;
  }

  auto setZero = []( real64 & val ){ val = 0.0; };
  LvArray::forValuesInSlice( phaseFraction.derivs, setZero );
  LvArray::forValuesInSlice( phaseCompFraction.derivs, setZero );

  if( Y < solubility )
  {
    // liquid phase only

    // 1) Compute phase fractions

    phaseFraction.value[m_phaseLiquidIndex] = 1.0;
    phaseFraction.value[m_phaseGasIndex] = 0.0;

    // 2) Compute phase component fractions

    phaseCompFraction.value[m_phaseGasIndex][m_CO2Index] = 1.0;
    phaseCompFraction.value[m_phaseGasIndex][m_waterIndex] = 0.0;
    for( localIndex ic = 0; ic < 2; ++ic )
    {
      phaseCompFraction.value[m_phaseLiquidIndex][ic] = compFraction[ic];

      for( localIndex jc = 0; jc < 2; ++jc )
      {
        phaseCompFraction.derivs[m_phaseLiquidIndex][ic][Deriv::dC+jc] = (ic == jc ) ? 1.0 : 0.0;
        phaseCompFraction.derivs[m_phaseGasIndex][ic][Deriv::dC+jc] = 0.0;
      }
    }
  }
  else
  {
    // two-phase flow

    // 1) Compute phase fractions

    // liquid phase fraction = (Csat + W) / (C + W) = (Csat/W + 1) / (C/W + 1)
    real64 const onePlusYInv = 1.0 / ( 1.0 + Y );
    phaseFraction.value[m_phaseLiquidIndex] = (solubility + 1.0) * onePlusYInv;

    phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dP] = solubilityDeriv[0] * onePlusYInv;
    phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dT] = solubilityDeriv[1] * onePlusYInv;
    phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dC+m_CO2Index] =
      -dY_dCompFrac[m_CO2Index] * phaseFraction.value[m_phaseLiquidIndex] * onePlusYInv;
    phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dC+m_waterIndex] =
      -dY_dCompFrac[m_waterIndex] * phaseFraction.value[m_phaseLiquidIndex] * onePlusYInv;

    phaseFraction.value[m_phaseGasIndex] = 1.0 - phaseFraction.value[m_phaseLiquidIndex];

    phaseFraction.derivs[m_phaseGasIndex][Deriv::dP] = -phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dP];
    phaseFraction.derivs[m_phaseGasIndex][Deriv::dT] = -phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dT];
    phaseFraction.derivs[m_phaseGasIndex][Deriv::dC+m_CO2Index] = -phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dC+m_CO2Index];
    phaseFraction.derivs[m_phaseGasIndex][Deriv::dC+m_waterIndex] = -phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dC+m_waterIndex];

    // 2) Compute phase component fractions

    // liquid phase composition  CO2 = Csat / (Csat + W) = (Csat/W) / (Csat/W + 1)
    real64 const onePlusSolubilityInv = 1.0 / ( 1.0 + solubility );
    phaseCompFraction.value[m_phaseLiquidIndex][m_CO2Index] = solubility * onePlusSolubilityInv;

    phaseCompFraction.derivs[m_phaseLiquidIndex][m_CO2Index][Deriv::dP] = solubilityDeriv[0] * (onePlusSolubilityInv*onePlusSolubilityInv);
    phaseCompFraction.derivs[m_phaseLiquidIndex][m_CO2Index][Deriv::dT] = solubilityDeriv[1] * (onePlusSolubilityInv*onePlusSolubilityInv);

    phaseCompFraction.value[m_phaseLiquidIndex][m_waterIndex] = 1.0 - phaseCompFraction.value[m_phaseLiquidIndex][m_CO2Index];

    phaseCompFraction.derivs[m_phaseLiquidIndex][m_waterIndex][Deriv::dP] = -phaseCompFraction.derivs[m_phaseLiquidIndex][m_CO2Index][Deriv::dP];
    phaseCompFraction.derivs[m_phaseLiquidIndex][m_waterIndex][Deriv::dT] = -phaseCompFraction.derivs[m_phaseLiquidIndex][m_CO2Index][Deriv::dT];

    // gas phase composition  CO2 = 1.0

    phaseCompFraction.value[m_phaseGasIndex][m_CO2Index]   = 1.0;
    phaseCompFraction.value[m_phaseGasIndex][m_waterIndex] = 0.0;

    phaseCompFraction.derivs[m_phaseGasIndex][m_CO2Index][Deriv::dP]   = 0.0;
    phaseCompFraction.derivs[m_phaseGasIndex][m_waterIndex][Deriv::dT] = 0.0;
    phaseCompFraction.derivs[m_phaseGasIndex][m_CO2Index][Deriv::dP]   = 0.0;
    phaseCompFraction.derivs[m_phaseGasIndex][m_waterIndex][Deriv::dT] = 0.0;
    // phaseCompFraction does not depend on globalComponentFraction

  }
}

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2SOLUBILITY_HPP_
