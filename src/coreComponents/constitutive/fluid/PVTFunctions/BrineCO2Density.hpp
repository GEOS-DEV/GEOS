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

#include "PVTFunctionBase.hpp"

#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"
#include "functions/TableFunction.hpp"

namespace geosx
{

namespace constitutive
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
                        arraySlice1d< real64 const > const & dPhaseComposition_dPressure,
                        arraySlice1d< real64 const > const & dPhaseComposition_dTemperature,
                        arraySlice2d< real64 const > const & dPhaseComposition_dGlobalCompFraction,
                        real64 & value,
                        real64 & dValue_dPressure,
                        real64 & dValue_dTemperature,
                        arraySlice1d< real64 > const & dValue_dGlobalCompFraction,
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

  BrineCO2Density( string_array const & inputPara,
                   string_array const & componentNames,
                   array1d< real64 > const & componentMolarWeight );

  ~BrineCO2Density() override {}

  static string catalogName() { return "BrineCO2Density"; }

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

  void makeTable( string_array const & inputPara );

  void calculateBrineDensity( PVTProps::PTTableCoordinates const & tableCoords,
                              real64 const & salinity,
                              array1d< real64 > const & densities );

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
                                     arraySlice1d< real64 const > const & dPhaseComposition_dPressure,
                                     arraySlice1d< real64 const > const & dPhaseComposition_dTemperature,
                                     arraySlice2d< real64 const > const & dPhaseComposition_dGlobalCompFraction,
                                     real64 & value,
                                     real64 & dValue_dPressure,
                                     real64 & dValue_dTemperature,
                                     arraySlice1d< real64 > const & dValue_dGlobalCompFraction,
                                     bool useMass ) const
{
  // this method implements the method proposed by E. Garcia (2001)

  // these coefficients come from equation (2) from Garcia (2001)
  constexpr real64 a = 37.51;
  constexpr real64 b = -9.585e-2;
  constexpr real64 c = 8.740e-4;
  constexpr real64 d = -5.044e-7;

  real64 const input[2] = { pressure, temperature };
  real64 density = 0.0;
  real64 densityDeriv[2]{};
  m_brineDensityTable.compute( input, density, densityDeriv );

  // equation (2) from Garcia (2001)
  real64 const squaredTemp = temperature * temperature;
  real64 const V = (  a
                      + b * temperature
                      + c * squaredTemp
                      + d * squaredTemp * temperature ) * 1e-6;
  real64 const dV_dTemp = ( b
                            + 2 * c * temperature
                            + 3 * d * squaredTemp ) * 1e-6;

  // CO2 concentration
  real64 const wMwInv = 1.0 / m_componentMolarWeight[m_waterIndex];
  real64 const oneMinusCO2PhaseCompInv = 1.0 / ( 1.0 - phaseComposition[m_CO2Index] );
  real64 const oneMinusCO2PhaseCompInvSquared = oneMinusCO2PhaseCompInv * oneMinusCO2PhaseCompInv;
  real64 const coef = wMwInv * phaseComposition[m_CO2Index] * oneMinusCO2PhaseCompInv;
  real64 const dCoef_dPres = wMwInv * dPhaseComposition_dPressure[m_CO2Index] * oneMinusCO2PhaseCompInvSquared;
  real64 const dCoef_dTemp = wMwInv * dPhaseComposition_dTemperature[m_CO2Index] * oneMinusCO2PhaseCompInvSquared;
  real64 dCoef_dComp[2]{};
  dCoef_dComp[m_CO2Index] = wMwInv * dPhaseComposition_dGlobalCompFraction[m_CO2Index][m_CO2Index] * oneMinusCO2PhaseCompInvSquared;
  dCoef_dComp[m_waterIndex] = wMwInv * dPhaseComposition_dGlobalCompFraction[m_CO2Index][m_waterIndex] * oneMinusCO2PhaseCompInvSquared;

  real64 const conc = coef * density;
  real64 const dConc_dPres = dCoef_dPres * density + coef * densityDeriv[0];
  real64 const dConc_dTemp = dCoef_dTemp * density + coef * densityDeriv[1];
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
  // equation (1) from Garcia (2001)
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
    dValue_dGlobalCompFraction[m_CO2Index] = m_componentMolarWeight[m_CO2Index] * dConc_dComp[m_CO2Index]
                                             - dConcDensVol_dComp[m_CO2Index];
    dValue_dGlobalCompFraction[m_waterIndex] = m_componentMolarWeight[m_CO2Index] * dConc_dComp[m_waterIndex]
                                               - dConcDensVol_dComp[m_waterIndex];
  }
  else
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
    dValue_dGlobalCompFraction[m_CO2Index] = dConc_dComp[m_CO2Index]
                                             - dConcDensVol_dComp[m_CO2Index] / m_componentMolarWeight[m_waterIndex];
    dValue_dGlobalCompFraction[m_waterIndex] = dConc_dComp[m_waterIndex]
                                               - dConcDensVol_dComp[m_waterIndex] / m_componentMolarWeight[m_waterIndex];
  }
}

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINECO2DENSITY_HPP_
