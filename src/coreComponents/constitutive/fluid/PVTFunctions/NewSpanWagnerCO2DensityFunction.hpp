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
 * @file SpanWagnerCO2Density.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_SPANWAGNERCO2DENSITY_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_SPANWAGNERCO2DENSITY_HPP_

#include "NewPVTFunctionBase.hpp"

#include "managers/Functions/TableFunction.hpp"

namespace geosx
{

namespace PVTProps
{

class SpanWagnerCO2DensityUpdate final : public PVTFunctionBaseUpdate
{
public:

  SpanWagnerCO2DensityUpdate( arrayView1d< string const > const & componentNames,
                              arrayView1d< real64 const > const & componentMolarWeight,
                              TableFunction * CO2DensityTable,
                              localIndex const CO2Index )
    : PVTFunctionBaseUpdate( componentNames,
                             componentMolarWeight ),
    m_CO2DensityTable( CO2DensityTable->createKernelWrapper() ),
    m_CO2Index( CO2Index )
  {}

  /// Default copy constructor
  SpanWagnerCO2DensityUpdate( SpanWagnerCO2DensityUpdate const & ) = default;

  /// Default move constructor
  SpanWagnerCO2DensityUpdate( SpanWagnerCO2DensityUpdate && ) = default;

  /// Deleted copy assignment operator
  SpanWagnerCO2DensityUpdate & operator=( SpanWagnerCO2DensityUpdate const & ) = delete;

  /// Deleted move assignment operator
  SpanWagnerCO2DensityUpdate & operator=( SpanWagnerCO2DensityUpdate && ) = delete;

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
  TableFunction::KernelWrapper const m_CO2DensityTable;

  /// Index of the CO2 component
  localIndex const m_CO2Index;

};

class SpanWagnerCO2Density : public PVTFunctionBase
{
public:

  SpanWagnerCO2Density( array1d< string > const & inputPara,
                        array1d< string > const & componentNames,
                        array1d< real64 > const & componentMolarWeight );

  ~SpanWagnerCO2Density() override {}

  static string catalogName() { return "NewSpanWagnerCO2Density"; }

  virtual string getCatalogName() const override final { return catalogName(); }

  virtual PVTFunctionType functionType() const override
  {
    return PVTFunctionType::DENSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = SpanWagnerCO2DensityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  static
  void calculateCO2Density( array1d< array1d< real64 > > const & coordinates,
                            array1d< real64 > const & values );

private:

  void makeTable( array1d< string > const & inputPara );

  static
  void spanWagnerCO2DensityFunction( real64 const & T,
                                     real64 const & P,
                                     real64 & rho,
                                     real64 (*f)( real64 const & x1, real64 const & x2, real64 const & x3 ) );

  /// Table with CO2 density tabulated as a function of (P,T)
  TableFunction * m_CO2DensityTable;

  /// Index of the CO2 component
  localIndex m_CO2Index;

};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SpanWagnerCO2DensityUpdate::compute( real64 const & pressure,
                                          real64 const & temperature,
                                          arraySlice1d< real64 const > const & phaseComposition,
                                          real64 & value,
                                          real64 & dValue_dPressure,
                                          real64 & dValue_dTemperature,
                                          arraySlice1d< real64 > const & dValue_dPhaseComposition,
                                          bool useMass ) const
{
  GEOSX_UNUSED_VAR( phaseComposition );

  real64 const input[2] = { pressure, temperature };
  real64 densityDeriv[2]{};
  m_CO2DensityTable.compute( input, value, densityDeriv );
  dValue_dPressure = densityDeriv[0];
  dValue_dTemperature = densityDeriv[1];

  if( !useMass )
  {
    value /= m_componentMolarWeight[m_CO2Index];
    dValue_dPressure /= m_componentMolarWeight[m_CO2Index];
    dValue_dTemperature /= m_componentMolarWeight[m_CO2Index];
  }

  for( real64 & val : dValue_dPhaseComposition )
  {
    val = 0.0;
  }
}

} // end namespace PVTProps

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_SPANWAGNERCO2DENSITY_HPP_
