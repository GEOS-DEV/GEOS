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
 * @file BrineCO2Density.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINECO2DENSITY_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINECO2DENSITY_HPP_

#include "NewPVTFunctionBase.hpp"

#include "managers/Functions/TableFunction.hpp"

namespace geosx
{

namespace PVTProps
{

class BrineCO2DensityUpdate final : public PVTFunctionBaseUpdate
{
public:

  BrineCO2DensityUpdate( arrayView1d< string const > const & componentNames,
                         arrayView1d< real64 const > const & componentMolarWeight,
                         TableFunction * brineDensityTable,
                         localIndex const CO2Index,
                         localIndex const waterIndex )
    : PVTFunctionBaseUpdate( componentNames,
                             componentMolarWeight ),
    m_brineDensityTable( brineDensityTable->createKernelWrapper() ),
    m_CO2Index( CO2Index ),
    m_waterIndex( waterIndex )
  {}

  /// Default copy constructor
  BrineCO2DensityUpdate( BrineCO2DensityUpdate const & ) = default;

  /// Default move constructor
  BrineCO2DensityUpdate( BrineCO2DensityUpdate && ) = default;

  /// Deleted copy assignment operator
  BrineCO2DensityUpdate & operator=( BrineCO2DensityUpdate const & ) = delete;

  /// Deleted move assignment operator
  BrineCO2DensityUpdate & operator=( BrineCO2DensityUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void compute( real64 const & pressure,
                        real64 const & temperature,
                        arraySlice1d< real64 const > const & phaseComposition,
                        real64 & value,
                        real64 & dValue_dPressure,
                        real64 & dValue_dTemperature,
                        arraySlice1d< real64 > const & dValue_dPhaseComposition,
                        bool useMass = 0 ) const override;

protected:

  /// Table with brine density tabulated as a function (P,T,sal)
  TableFunction::KernelWrapper const m_brineDensityTable;

  /// Index of the CO2 component
  localIndex const m_CO2Index;

  /// Index of the water component
  localIndex const m_waterIndex;

};

class BrineCO2Density : public PVTFunctionBase
{
public:

  BrineCO2Density( array1d< string > const & inputPara,
                   array1d< string > const & componentNames,
                   array1d< real64 > const & componentMolarWeight );

  ~BrineCO2Density() override {}

  static string catalogName() { return "NewBrineCO2Density"; }

  virtual string getCatalogName() const override final { return catalogName(); }

  virtual PVTFunctionType functionType() const override
  {
    return PVTFunctionType::DENSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = BrineCO2DensityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();


private:

  void makeTable( array1d< string > const & inputPara );

  void calculateBrineDensity( array1d< array1d< real64 > > const & coordinates,
                              real64 const & salinity,
                              array1d< real64 > const & values );

  /// Table with brine density tabulated as a function of (P,T,sal)
  TableFunction * m_brineDensityTable;

  /// Index of the CO2 phase
  localIndex m_CO2Index;

  /// Index of the water phase
  localIndex m_waterIndex;

};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void BrineCO2DensityUpdate::compute( real64 const & pressure,
                                     real64 const & temperature,
                                     arraySlice1d< real64 const > const & phaseComposition,
                                     real64 & value,
                                     real64 & dValue_dPressure,
                                     real64 & dValue_dTemperature,
                                     arraySlice1d< real64 > const & dValue_dPhaseComposition,
                                     bool useMass ) const
{
  constexpr real64 a = 37.51;
  constexpr real64 b = -9.585e-2;
  constexpr real64 c = 8.740e-4;
  constexpr real64 d = -5.044e-7;

  real64 const input[2] = { pressure, temperature };
  real64 density = 0.0;
  real64 densityDeriv[2]{};
  m_brineDensityTable.compute( input, density, densityDeriv );

  real64 const squaredTemp = temperature * temperature;
  real64 const V = (  a
                      + b * temperature
                      + c * squaredTemp
                      + d * squaredTemp * temperature ) * 1e-6;
  real64 const dV_dTemp = ( b
                            + 2 * c * temperature
                            + 3 * d * squaredTemp ) * 1e-6;

  // CO2 concentration
  // C = X * den / (waterMW * (1.0 - X));
  real64 const denom = ( m_componentMolarWeight[m_waterIndex] * ( 1.0 - phaseComposition[m_CO2Index] ) );
  real64 const coef = phaseComposition[m_CO2Index] / denom;
  real64 dCoef_dComp[2]{};
  dCoef_dComp[m_CO2Index] = 1.0 / denom;
  dCoef_dComp[m_waterIndex] = -m_componentMolarWeight[m_waterIndex] * phaseComposition[m_CO2Index] / (denom * denom);

  real64 const conc = coef * density;
  real64 const dConc_dPres = coef * densityDeriv[0];
  real64 const dConc_dTemp = coef * densityDeriv[1];
  real64 dConc_dComp[2]{};
  dConc_dComp[m_CO2Index] = dCoef_dComp[m_CO2Index] * density;
  dConc_dComp[m_waterIndex] = dCoef_dComp[m_waterIndex] * density;

  // CO2 concentration times density times vol
  real64 const concDensVol = conc * density * V;
  real64 const dConcDensVol_dPres = ( dConc_dPres * density + conc * densityDeriv[0] ) * V;
  real64 const dConcDensVol_dTemp = ( dConc_dTemp * density + conc * densityDeriv[1] ) * V
                                    + conc * density * dV_dTemp;
  real64 dConcDensVol_dComp[2]{};
  dConcDensVol_dComp[m_CO2Index] = dConc_dComp[m_CO2Index] * density * V;
  dConcDensVol_dComp[m_waterIndex] = dConc_dComp[m_waterIndex] * density * V;

  // Brine density
  // value = den + CO2MW * C - C * den * V;
  if( useMass )
  {
    value = density
            + m_componentMolarWeight[m_CO2Index] * conc
            - concDensVol;
    dValue_dPressure = densityDeriv[0]
                       + m_componentMolarWeight[m_CO2Index] * dConc_dPres
                       - dConcDensVol_dPres;
    dValue_dTemperature = densityDeriv[1]
                          + m_componentMolarWeight[m_CO2Index] * dConc_dTemp
                          - dConcDensVol_dTemp;
    dValue_dPhaseComposition[m_CO2Index] = m_componentMolarWeight[m_CO2Index] * dConc_dComp[m_CO2Index]
                                           - dConcDensVol_dComp[m_CO2Index];
    dValue_dPhaseComposition[m_waterIndex] = m_componentMolarWeight[m_CO2Index] * dConc_dComp[m_waterIndex]
                                             - dConcDensVol_dComp[m_waterIndex];
  }
  else // value = den / waterMW + C - C * den * V / waterMW;
  {
    value = density / m_componentMolarWeight[m_waterIndex]
            + conc
            - concDensVol / m_componentMolarWeight[m_waterIndex];
    dValue_dPressure = densityDeriv[0]  / m_componentMolarWeight[m_waterIndex]
                       + dConc_dPres
                       - dConcDensVol_dPres / m_componentMolarWeight[m_waterIndex];
    dValue_dTemperature = densityDeriv[1]  / m_componentMolarWeight[m_waterIndex]
                          + dConc_dTemp
                          - dConcDensVol_dTemp  / m_componentMolarWeight[m_waterIndex];
    dValue_dPhaseComposition[m_CO2Index] = dConc_dComp[m_CO2Index]
                                           - dConcDensVol_dComp[m_CO2Index] / m_componentMolarWeight[m_waterIndex];
    dValue_dPhaseComposition[m_waterIndex] = dConc_dComp[m_waterIndex]
                                             - dConcDensVol_dComp[m_waterIndex] / m_componentMolarWeight[m_waterIndex];
  }
}

} // end namespace PVTProps

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINECO2DENSITY_HPP_
