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
 * @file SpanWagnerCO2Density.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_SPANWAGNERCO2DENSITY_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_SPANWAGNERCO2DENSITY_HPP_

#include "PVTFunctionBase.hpp"

#include "constitutive/fluid/layouts.hpp"
#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{

namespace constitutive
{

namespace PVTProps
{

class SpanWagnerCO2DensityUpdate final : public PVTFunctionBaseUpdate
{
public:

  SpanWagnerCO2DensityUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                              TableFunction const & CO2DensityTable,
                              localIndex const CO2Index )
    : PVTFunctionBaseUpdate( componentMolarWeight ),
    m_CO2DensityTable( CO2DensityTable.createKernelWrapper() ),
    m_CO2Index( CO2Index )
  {}

  template< int USD1 >
  GEOS_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                real64 & value,
                bool useMass ) const;

  template< int USD1, int USD2, int USD3 >
  GEOS_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                arraySlice2d< real64 const, USD2 > const & dPhaseComposition,
                real64 & value,
                arraySlice1d< real64, USD3 > const & dValue,
                bool useMass ) const;

  virtual void move( LvArray::MemorySpace const space, bool const touch ) override
  {
    PVTFunctionBaseUpdate::move( space, touch );
    m_CO2DensityTable.move( space, touch );
  }

protected:

  /// Table with brine density tabulated as a function (P,T,sal)
  TableFunction::KernelWrapper m_CO2DensityTable;

  /// Index of the CO2 component
  localIndex m_CO2Index;

};

class SpanWagnerCO2Density : public PVTFunctionBase
{
public:

  SpanWagnerCO2Density( string const &,
                        string_array const & inputParams,
                        string_array const & componentNames,
                        array1d< real64 > const & componentMolarWeight );

  static string catalogName() { return "SpanWagnerCO2Density"; }

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
  KernelWrapper createKernelWrapper() const;

  static
  void calculateCO2Density( string const & functionName,
                            real64 const & tolerance,
                            PVTProps::PTTableCoordinates const & tableCoords,
                            array1d< real64 > const & densities );

private:

  /// Table with CO2 density tabulated as a function of (P,T)
  TableFunction const * m_CO2DensityTable;

  /// Index of the CO2 component
  localIndex m_CO2Index;

};

template< int USD1 >
GEOS_HOST_DEVICE
void SpanWagnerCO2DensityUpdate::compute( real64 const & pressure,
                                          real64 const & temperature,
                                          arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                          real64 & value,
                                          bool useMass ) const
{
  GEOS_UNUSED_VAR( phaseComposition );

  real64 const input[2] = { pressure, temperature };
  value = m_CO2DensityTable.compute( input );

  if( !useMass )
  {
    value /= m_componentMolarWeight[m_CO2Index];
  }
}

template< int USD1, int USD2, int USD3 >
GEOS_HOST_DEVICE
void SpanWagnerCO2DensityUpdate::compute( real64 const & pressure,
                                          real64 const & temperature,
                                          arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                          arraySlice2d< real64 const, USD2 > const & dPhaseComposition,
                                          real64 & value,
                                          arraySlice1d< real64, USD3 > const & dValue,
                                          bool useMass ) const
{
  GEOS_UNUSED_VAR( phaseComposition, dPhaseComposition );

  using Deriv = multifluid::DerivativeOffset;

  real64 const input[2] = { pressure, temperature };
  real64 densityDeriv[2]{};
  value = m_CO2DensityTable.compute( input, densityDeriv );

  LvArray::forValuesInSlice( dValue, []( real64 & val ){ val = 0.0; } );
  dValue[Deriv::dP] = densityDeriv[0];
  dValue[Deriv::dT] = densityDeriv[1];

  if( !useMass )
  {
    value /= m_componentMolarWeight[m_CO2Index];
    dValue[Deriv::dP] /= m_componentMolarWeight[m_CO2Index];
    dValue[Deriv::dT] /= m_componentMolarWeight[m_CO2Index];
  }

}

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geos

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_SPANWAGNERCO2DENSITY_HPP_
