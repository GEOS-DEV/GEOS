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
 * @file CO2Solubility.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2SOLUBILITY_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2SOLUBILITY_HPP_

#include "FlashModelBase.hpp"

#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"
#include "managers/Functions/TableFunction.hpp"

namespace geosx
{

namespace constitutive
{

namespace PVTProps
{

constexpr real64 minForDivision = 1e-10;

class CO2SolubilityUpdate final : public FlashModelBaseUpdate
{
public:

  CO2SolubilityUpdate( arrayView1d< string const > const & componentNames,
                       arrayView1d< real64 const > const & componentMolarWeight,
                       TableFunction * CO2SolubilityTable,
                       localIndex const CO2Index,
                       localIndex const waterIndex,
                       localIndex const phaseGasIndex,
                       localIndex const phaseLiquidIndex )
    : FlashModelBaseUpdate( componentNames,
                            componentMolarWeight ),
    m_CO2SolubilityTable( CO2SolubilityTable->createKernelWrapper() ),
    m_CO2Index( CO2Index ),
    m_waterIndex( waterIndex ),
    m_phaseGasIndex( phaseGasIndex ),
    m_phaseLiquidIndex( phaseLiquidIndex )
  {}

  /// Default copy constructor
  CO2SolubilityUpdate( CO2SolubilityUpdate const & ) = default;

  /// Default move constructor
  CO2SolubilityUpdate( CO2SolubilityUpdate && ) = default;

  /// Deleted copy assignment operator
  CO2SolubilityUpdate & operator=( CO2SolubilityUpdate const & ) = delete;

  /// Deleted move assignment operator
  CO2SolubilityUpdate & operator=( CO2SolubilityUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void compute( real64 const & pressure,
                        real64 const & temperature,
                        arraySlice1d< real64 const > const & compFraction,
                        arraySlice1d< real64 > const & phaseFraction,
                        arraySlice1d< real64 > const & dPhaseFraction_dPressure,
                        arraySlice1d< real64 > const & dPhaseFraction_dTemperature,
                        arraySlice2d< real64 > const & dPhaseFraction_dCompFraction,
                        arraySlice2d< real64 > const & phaseCompFraction,
                        arraySlice2d< real64 > const & dPhaseCompFraction_dPressure,
                        arraySlice2d< real64 > const & dPhaseCompFraction_dTemperature,
                        arraySlice3d< real64 > const & dPhaseCompFraction_dCompFraction ) const;

protected:

  /// Table with CO2 solubility tabulated as a function (P,T)
  TableFunction::KernelWrapper const m_CO2SolubilityTable;

  /// Index of the CO2 phase
  localIndex const m_CO2Index;

  /// Index of the water phase
  localIndex const m_waterIndex;

  /// Index of the gas phase
  localIndex const m_phaseGasIndex;

  /// Index of the liquid phase
  localIndex const m_phaseLiquidIndex;

};

class CO2Solubility : public FlashModelBase
{
public:

  CO2Solubility( string_array const & inputPara,
                 string_array const & phaseNames,
                 string_array const & componentNames,
                 array1d< real64 > const & componentMolarWeight );

  ~CO2Solubility() override {}

  static string catalogName() { return "CO2Solubility"; }

  virtual string getCatalogName() const override final { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = CO2SolubilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

private:

  void makeTable( string_array const & inputPara );

  void calculateCO2Solubility( real64 const & tolerance,
                               PVTProps::PTTableCoordinates const & tableCoords,
                               real64 const & salinity,
                               array1d< real64 > const & values );

  void CO2SolubilityFunction( real64 const & tolerance,
                              real64 const & T,
                              real64 const & P,
                              real64 & V_r,
                              real64 (*f)( real64 const & x1, real64 const & x2, real64 const & x3 ) );

  /// Table to compute solubility as a function of pressure and temperature
  TableFunction * m_CO2SolubilityTable;

  /// Index of the CO2 component
  localIndex m_CO2Index;

  /// Index of the water component
  localIndex m_waterIndex;

  /// Index of the gas phase
  localIndex m_phaseGasIndex;

  /// Index of the liquid phase
  localIndex m_phaseLiquidIndex;
};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void CO2SolubilityUpdate::compute( real64 const & pressure,
                                   real64 const & temperature,
                                   arraySlice1d< real64 const > const & compFraction,
                                   arraySlice1d< real64 > const & phaseFraction,
                                   arraySlice1d< real64 > const & dPhaseFraction_dPressure,
                                   arraySlice1d< real64 > const & dPhaseFraction_dTemperature,
                                   arraySlice2d< real64 > const & dPhaseFraction_dCompFraction,
                                   arraySlice2d< real64 > const & phaseCompFraction,
                                   arraySlice2d< real64 > const & dPhaseCompFraction_dPressure,
                                   arraySlice2d< real64 > const & dPhaseCompFraction_dTemperature,
                                   arraySlice3d< real64 > const & dPhaseCompFraction_dCompFraction ) const
{
  // solubility mol/kg(water)  X = Csat/W
  real64 const input[2] = { pressure, temperature };
  real64 solubility = 0.0;
  real64 solubilityDeriv[2]{};
  m_CO2SolubilityTable.compute( input, solubility, solubilityDeriv );

  solubility *= m_componentMolarWeight[m_waterIndex];
  for( localIndex i = 0; i < 2; ++i )
  {
    solubilityDeriv[i] *= m_componentMolarWeight[m_waterIndex];
  }

  real64 Y = 0.0;
  real64 dY_dCompFrac[2]{};

  // Y = C/W = z/(1-z)

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

  if( Y < solubility )
  {
    // liquid phase only

    // 1) Compute phase fractions

    phaseFraction[m_phaseLiquidIndex] = 1.0;
    phaseFraction[m_phaseGasIndex] = 0.0;
    for( real64 & val : dPhaseFraction_dPressure )
    {
      val = 0.0;
    }
    for( real64 & val : dPhaseFraction_dTemperature )
    {
      val = 0.0;
    }
    for( real64 & val : dPhaseFraction_dCompFraction )
    {
      val = 0.0;
    }

    // 2) Compute phase component fractions

    for( localIndex ic = 0; ic < 2; ++ic )
    {
      phaseCompFraction[m_phaseLiquidIndex][ic] = compFraction[ic];
      // the two following lines are not present in Yue's code, unclear if this will have some consequences
      phaseCompFraction[m_phaseGasIndex][m_CO2Index] = 1.0;
      phaseCompFraction[m_phaseGasIndex][m_waterIndex] = 0.0;
      for( localIndex jc = 0; jc < 2; ++jc )
      {
        dPhaseCompFraction_dCompFraction[m_phaseLiquidIndex][ic][jc] = (ic == jc ) ? 1.0 : 0.0;
        dPhaseCompFraction_dCompFraction[m_phaseGasIndex][ic][jc] = 0.0;
      }
    }
    for( real64 & val : dPhaseCompFraction_dPressure )
    {
      val = 0.0;
    }
    for( real64 & val : dPhaseCompFraction_dTemperature )
    {
      val = 0.0;
    }

  }
  else
  {
    // two-phase flow

    // 1) Compute phase fractions

    // liquid phase fraction = (Csat + W) / (C + W) = (Csat/W + 1) / (C/W + 1)
    real64 const onePlusYInv = 1.0 / ( 1.0 + Y );
    phaseFraction[m_phaseLiquidIndex] = (solubility + 1.0) * onePlusYInv;
    dPhaseFraction_dPressure[m_phaseLiquidIndex] = solubilityDeriv[0] * onePlusYInv;
    dPhaseFraction_dTemperature[m_phaseLiquidIndex] = solubilityDeriv[1] * onePlusYInv;
    dPhaseFraction_dCompFraction[m_phaseLiquidIndex][m_CO2Index] =
      -dY_dCompFrac[m_CO2Index] * phaseFraction[m_phaseLiquidIndex] * onePlusYInv;
    dPhaseFraction_dCompFraction[m_phaseLiquidIndex][m_waterIndex] =
      -dY_dCompFrac[m_waterIndex] * phaseFraction[m_phaseLiquidIndex] * onePlusYInv;

    phaseFraction[m_phaseGasIndex] = 1.0 - phaseFraction[m_phaseLiquidIndex];
    dPhaseFraction_dPressure[m_phaseGasIndex] = -dPhaseFraction_dPressure[m_phaseLiquidIndex];
    dPhaseFraction_dTemperature[m_phaseGasIndex] = -dPhaseFraction_dTemperature[m_phaseLiquidIndex];
    dPhaseFraction_dCompFraction[m_phaseGasIndex][m_CO2Index] = -dPhaseFraction_dCompFraction[m_phaseLiquidIndex][m_CO2Index];
    dPhaseFraction_dCompFraction[m_phaseGasIndex][m_waterIndex] = -dPhaseFraction_dCompFraction[m_phaseLiquidIndex][m_waterIndex];

    // 2) Compute phase component fractions

    // liquid phase composition  CO2 = Csat / (Csat + W) = (Csat/W) / (Csat/W + 1)
    real64 const onePlusSolubilityInv = 1.0 / ( 1.0 + solubility );
    phaseCompFraction[m_phaseLiquidIndex][m_CO2Index] = solubility * onePlusSolubilityInv;
    dPhaseCompFraction_dPressure[m_phaseLiquidIndex][m_CO2Index] = solubilityDeriv[0] * (onePlusSolubilityInv*onePlusSolubilityInv);
    dPhaseCompFraction_dTemperature[m_phaseLiquidIndex][m_CO2Index] = solubilityDeriv[1] * (onePlusSolubilityInv*onePlusSolubilityInv);

    phaseCompFraction[m_phaseLiquidIndex][m_waterIndex] = 1.0 - phaseCompFraction[m_phaseLiquidIndex][m_CO2Index];
    dPhaseCompFraction_dPressure[m_phaseLiquidIndex][m_waterIndex] = -dPhaseCompFraction_dPressure[m_phaseLiquidIndex][m_CO2Index];
    dPhaseCompFraction_dTemperature[m_phaseLiquidIndex][m_waterIndex] = -dPhaseCompFraction_dTemperature[m_phaseLiquidIndex][m_CO2Index];

    // gas phase composition  CO2 = 1.0

    phaseCompFraction[m_phaseGasIndex][m_CO2Index] = 1.0;
    phaseCompFraction[m_phaseGasIndex][m_waterIndex] = 0.0;
    dPhaseCompFraction_dPressure[m_phaseGasIndex][m_CO2Index] = 0.0;
    dPhaseCompFraction_dPressure[m_phaseGasIndex][m_waterIndex] = 0.0;
    dPhaseCompFraction_dTemperature[m_phaseGasIndex][m_CO2Index] = 0.0;
    dPhaseCompFraction_dTemperature[m_phaseGasIndex][m_waterIndex] = 0.0;

    // phaseCompFraction does not depend on globalComponentFraction
    for( real64 & val : dPhaseCompFraction_dCompFraction )
    {
      val = 0.0;
    }
  }
}

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2SOLUBILITY_HPP_
