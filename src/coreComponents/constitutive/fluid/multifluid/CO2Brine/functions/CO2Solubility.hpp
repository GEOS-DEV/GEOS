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

  using PhaseProp = MultiFluidVar< real64, 3, multifluid::LAYOUT_PHASE, multifluid::LAYOUT_PHASE_DC >;
  using PhaseComp = MultiFluidVar< real64, 4, multifluid::LAYOUT_PHASE_COMP, multifluid::LAYOUT_PHASE_COMP_DC >;

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

  template< int USD1, int USD2, int USD3 >
  GEOS_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & compFraction,
                arraySlice1d< real64, USD2 > const & phaseFraction,
                arraySlice2d< real64, USD3 > const & phaseCompFraction ) const;

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

  /// Expected number of components
  static constexpr integer numComps = 2;

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
                 array1d< real64 > const & componentMolarWeight );

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

template< int USD1, int USD2, int USD3 >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
CO2SolubilityUpdate::compute( real64 const & pressure,
                              real64 const & temperature,
                              arraySlice1d< real64 const, USD1 > const & compFraction,
                              arraySlice1d< real64, USD2 > const & phaseFraction,
                              arraySlice2d< real64, USD3 > const & phaseCompFraction ) const
{
  // solubility mol/kg(water)  X = Csat/W
  real64 const input[2] = { pressure, temperature };
  real64 solubility = m_CO2SolubilityTable.compute( input );
  solubility *= m_componentMolarWeight[m_waterIndex];

  // Y = C/W = z/(1-z)
  real64 Y;

  if( compFraction[m_CO2Index] > 1.0 - minForDivision )
  {
    Y = compFraction[m_CO2Index] / minForDivision;
  }
  else
  {
    real64 const oneMinusCompFracInv = 1.0 / (1.0 - compFraction[m_CO2Index]);
    Y = compFraction[m_CO2Index] * oneMinusCompFracInv;
  }

  if( Y < solubility )
  {
    // liquid phase only

    // 1) Compute phase fractions

    phaseFraction[m_phaseLiquidIndex] = 1.0;
    phaseFraction[m_phaseGasIndex] = 0.0;

    // 2) Compute phase component fractions

    for( localIndex ic = 0; ic < 2; ++ic )
    {
      phaseCompFraction[m_phaseLiquidIndex][ic] = compFraction[ic];
      // the two following lines are not present in Yue's code, unclear if this will have some consequences
      phaseCompFraction[m_phaseGasIndex][m_CO2Index] = 1.0;
      phaseCompFraction[m_phaseGasIndex][m_waterIndex] = 0.0;
    }
  }
  else
  {
    // two-phase flow

    // 1) Compute phase fractions

    // liquid phase fraction = (Csat + W) / (C + W) = (Csat/W + 1) / (C/W + 1)
    real64 const onePlusYInv = 1.0 / ( 1.0 + Y );
    phaseFraction[m_phaseLiquidIndex] = (solubility + 1.0) * onePlusYInv;
    phaseFraction[m_phaseGasIndex] = 1.0 - phaseFraction[m_phaseLiquidIndex];

    // 2) Compute phase component fractions

    // liquid phase composition  CO2 = Csat / (Csat + W) = (Csat/W) / (Csat/W + 1)
    real64 const onePlusSolubilityInv = 1.0 / ( 1.0 + solubility );
    phaseCompFraction[m_phaseLiquidIndex][m_CO2Index] = solubility * onePlusSolubilityInv;

    phaseCompFraction[m_phaseLiquidIndex][m_waterIndex] = 1.0 - phaseCompFraction[m_phaseLiquidIndex][m_CO2Index];

    // gas phase composition  CO2 = 1.0
    phaseCompFraction[m_phaseGasIndex][m_CO2Index] = 1.0;
    phaseCompFraction[m_phaseGasIndex][m_waterIndex] = 0.0;
  }
}

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

  // Solubility of CO2 is read from the tables in the form of moles of CO2 per kg of water
  // Solubility of water is read from the tables in the form of moles of water per kg of CO2
  real64 const input[2] = { pressure, temperature };

  real64 co2SolubilityDeriv[2]{};
  real64 watSolubilityDeriv[] = { 0.0, 0.0 };
  real64 co2Solubility = m_CO2SolubilityTable.compute( input, co2SolubilityDeriv );
  real64 watSolubility = 0.0;

  // Convert the solubility to mole/mole
  co2Solubility *= m_componentMolarWeight[m_waterIndex];
  watSolubility *= m_componentMolarWeight[m_CO2Index];
  for( integer ic = 0; ic < 2; ++ic )
  {
    co2SolubilityDeriv[ic] *= m_componentMolarWeight[m_waterIndex];
    watSolubilityDeriv[ic] *= m_componentMolarWeight[m_CO2Index];
  }

  real64 const z_co2 = compFraction[m_CO2Index];
  real64 const z_wat = compFraction[m_waterIndex];

  real64 const determinant = 1.0 - co2Solubility*watSolubility;
  real64 invDeterminant = 0.0;
  real64 invDeterminantDeriv[] = { 0.0, 0.0 };
  if( minForDivision < fabs( determinant ) )
  {
    invDeterminant = 1.0 / determinant;
    for( integer ic = 0; ic < 2; ic++ )
    {
      invDeterminantDeriv[ic] = invDeterminant*invDeterminant*(co2Solubility*watSolubilityDeriv[ic] + watSolubility*co2SolubilityDeriv[ic]);
    }
  }
  else
  {
    invDeterminant = 1.0 / minForDivision;
  }

  real64 x_co2 = co2Solubility * (z_wat - watSolubility * z_co2) * invDeterminant;
  real64 x_co2Deriv[4]{ 0.0, 0.0, 0.0, 0.0 };
  if( minForDivision < x_co2 )
  {
    // Pressure and temperature derivatives
    for( integer ic = 0; ic < 2; ic++ )
    {
      x_co2Deriv[ic] = co2SolubilityDeriv[ic] * (z_wat - watSolubility * z_co2) * invDeterminant
                       - co2Solubility * watSolubilityDeriv[ic] * z_co2 * invDeterminant
                       + co2Solubility * (z_wat - watSolubility * z_co2) * invDeterminantDeriv[ic];
    }
    // Composition derivatives
    x_co2Deriv[2+m_CO2Index] = -co2Solubility * watSolubility * invDeterminant;
    x_co2Deriv[2+m_waterIndex] = co2Solubility * invDeterminant;
  }
  else
  {
    x_co2 = 0.0;
  }

  real64 y_wat = watSolubility * (z_co2 - x_co2);
  real64 y_watDeriv[4]{ 0.0, 0.0, 0.0, 0.0 };
  if( minForDivision < y_wat )
  {
    // Pressure and temperature derivatives
    for( integer ic = 0; ic < 2; ic++ )
    {
      y_watDeriv[ic] = watSolubilityDeriv[ic] * (z_co2 - x_co2)
                       - watSolubility * x_co2Deriv[ic];
    }
    // Composition derivatives
    y_watDeriv[2+m_CO2Index] = watSolubility*(1.0 - x_co2Deriv[2+m_CO2Index]);
    y_watDeriv[2+m_waterIndex] = -watSolubility * x_co2Deriv[2+m_waterIndex];
  }
  else
  {
    y_wat = 0.0;
  }

  // Liquid and vapour phase fractions
  real64 const L = x_co2 + z_wat - y_wat;
  real64 const V = y_wat + z_co2 - x_co2; // = 1.0 - L;

  if( minForDivision < L && minForDivision < V )
  {
    // Two phases

    // 1) Compute phase fractions and derivatives

    real64 const dL_dP = x_co2Deriv[0] - y_watDeriv[0];
    real64 const dL_dT = x_co2Deriv[1] - y_watDeriv[1];
    real64 const dL_dzco2 = x_co2Deriv[2+m_CO2Index] - y_watDeriv[2+m_CO2Index];
    real64 const dL_dzwat = 1.0 + x_co2Deriv[2+m_waterIndex] - y_watDeriv[2+m_waterIndex];

    real64 const dV_dP = -dL_dP;
    real64 const dV_dT = -dL_dT;
    real64 const dV_dzco2 = -dL_dzco2;
    real64 const dV_dzwat = -dL_dzwat;

    phaseFraction.value[m_phaseLiquidIndex] = L;

    phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dP] = dL_dP;
    phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dT] = dL_dT;
    phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dC+m_CO2Index] = dL_dzco2;
    phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dC+m_waterIndex] = dL_dzwat;

    phaseFraction.value[m_phaseGasIndex] = V;

    phaseFraction.derivs[m_phaseGasIndex][Deriv::dP] = dV_dP;
    phaseFraction.derivs[m_phaseGasIndex][Deriv::dT] = dV_dT;
    phaseFraction.derivs[m_phaseGasIndex][Deriv::dC+m_CO2Index] = dV_dzco2;
    phaseFraction.derivs[m_phaseGasIndex][Deriv::dC+m_waterIndex] = dV_dzwat;

    // 2) Compute phase component fractions and derivatives
    real64 const invL2 = 1.0 / (L*L);
    real64 const invV2 = 1.0 / (V*V);

    phaseCompFraction.value[m_phaseLiquidIndex][m_CO2Index] = x_co2 / L;
    phaseCompFraction.derivs[m_phaseLiquidIndex][m_CO2Index][Deriv::dP] = (x_co2Deriv[0]*L - x_co2*dL_dP)*invL2;
    phaseCompFraction.derivs[m_phaseLiquidIndex][m_CO2Index][Deriv::dT] = (x_co2Deriv[1]*L - x_co2*dL_dT)*invL2;
    phaseCompFraction.derivs[m_phaseLiquidIndex][m_CO2Index][Deriv::dC+m_CO2Index] = (x_co2Deriv[2+m_CO2Index]*L - x_co2*dL_dzco2)*invL2;
    phaseCompFraction.derivs[m_phaseLiquidIndex][m_CO2Index][Deriv::dC+m_waterIndex] = (x_co2Deriv[2+m_waterIndex]*L - x_co2*dL_dzwat)*invL2;

    phaseCompFraction.value[m_phaseLiquidIndex][m_waterIndex] = (z_wat - y_wat) / L;
    phaseCompFraction.derivs[m_phaseLiquidIndex][m_waterIndex][Deriv::dP] = (-y_watDeriv[0]*L - (z_wat - y_wat)*dL_dP)*invL2;
    phaseCompFraction.derivs[m_phaseLiquidIndex][m_waterIndex][Deriv::dT] = (-y_watDeriv[1]*L - (z_wat - y_wat)*dL_dT)*invL2;
    phaseCompFraction.derivs[m_phaseLiquidIndex][m_waterIndex][Deriv::dC+m_CO2Index] = (-y_watDeriv[2+m_CO2Index]*L - (z_wat - y_wat)*dL_dzco2)*invL2;
    phaseCompFraction.derivs[m_phaseLiquidIndex][m_waterIndex][Deriv::dC+m_waterIndex] = ( (1.0-y_watDeriv[2+m_waterIndex])*L - (z_wat - y_wat)*dL_dzwat)*invL2;

    phaseCompFraction.value[m_phaseGasIndex][m_CO2Index] = (z_co2 - x_co2) / V;
    phaseCompFraction.derivs[m_phaseGasIndex][m_CO2Index][Deriv::dP] = (-x_co2Deriv[0]*V - (z_co2 - x_co2)*dV_dP)*invV2;
    phaseCompFraction.derivs[m_phaseGasIndex][m_CO2Index][Deriv::dT] = (-x_co2Deriv[1]*V - (z_co2 - x_co2)*dV_dT)*invV2;
    phaseCompFraction.derivs[m_phaseGasIndex][m_CO2Index][Deriv::dC+m_CO2Index] = ((1.0-x_co2Deriv[2+m_CO2Index])*V - (z_co2 - x_co2)*dV_dzco2)*invV2;
    phaseCompFraction.derivs[m_phaseGasIndex][m_CO2Index][Deriv::dC+m_waterIndex] = (-x_co2Deriv[2+m_waterIndex]*V - (z_co2 - x_co2)*dV_dzwat)*invV2;

    phaseCompFraction.value[m_phaseGasIndex][m_waterIndex] = y_wat / V;
    phaseCompFraction.derivs[m_phaseGasIndex][m_waterIndex][Deriv::dP] = (y_watDeriv[0]*V - y_wat*dV_dP)*invV2;
    phaseCompFraction.derivs[m_phaseGasIndex][m_waterIndex][Deriv::dT] = (y_watDeriv[1]*V - y_wat*dV_dT)*invV2;
    phaseCompFraction.derivs[m_phaseGasIndex][m_waterIndex][Deriv::dC+m_CO2Index] = (y_watDeriv[2+m_CO2Index]*V - y_wat*dV_dzco2)*invV2;
    phaseCompFraction.derivs[m_phaseGasIndex][m_waterIndex][Deriv::dC+m_waterIndex] = (y_watDeriv[2+m_waterIndex]*V - y_wat*dV_dzwat)*invV2;
  }
  else
  {
    // Single phase: Select the present phase
    integer const activePhase = minForDivision < L ? m_phaseLiquidIndex : m_phaseGasIndex;

    // Zero out everything to start
    auto setZero = []( real64 & val ){ val = 0.0; };
    LvArray::forValuesInSlice( phaseFraction.value, setZero );
    LvArray::forValuesInSlice( phaseCompFraction.value, setZero );
    LvArray::forValuesInSlice( phaseFraction.derivs, setZero );
    LvArray::forValuesInSlice( phaseCompFraction.derivs, setZero );

    // 1) Compute phase fractions

    phaseFraction.value[activePhase] = 1.0;

    // 2) Compute phase component fractions
    // Setup default values which will be overridden for the active phase
    phaseCompFraction.value[m_phaseGasIndex][m_CO2Index]   = 1.0;
    phaseCompFraction.value[m_phaseGasIndex][m_waterIndex] = 0.0;
    phaseCompFraction.value[m_phaseLiquidIndex][m_CO2Index]   = 0.0;
    phaseCompFraction.value[m_phaseLiquidIndex][m_waterIndex] = 1.0;
    // Set the global composition as the composition of the active phase
    for( integer ic = 0; ic < numComps; ++ic )
    {
      phaseCompFraction.value[activePhase][ic] = compFraction[ic];
      phaseCompFraction.derivs[activePhase][ic][Deriv::dC+ic] = 1.0;
    }
  }
}

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2SOLUBILITY_HPP_
