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
 * @file FenghourCO2Viscosity.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_FENGHOURCO2VISCOSITY_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_FENGHOURCO2VISCOSITY_HPP_

#include "PVTFunctionBase.hpp"

#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"
#include "functions/TableFunction.hpp"

namespace geosx
{

namespace constitutive
{

namespace PVTProps
{

class FenghourCO2ViscosityUpdate final : public PVTFunctionBaseUpdate
{
public:

  FenghourCO2ViscosityUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                              TableFunction const & CO2ViscosityTable )
    : PVTFunctionBaseUpdate( componentMolarWeight ),
    m_CO2ViscosityTable( CO2ViscosityTable.createKernelWrapper() )
  {}

  template< int USD1 >
  GEOSX_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                real64 & value,
                bool useMass ) const;

  template< int USD1, int USD2, int USD3, int USD4 >
  GEOSX_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                arraySlice1d< real64 const, USD2 > const & dPhaseComposition_dPressure,
                arraySlice1d< real64 const, USD2 > const & dPhaseComposition_dTemperature,
                arraySlice2d< real64 const, USD3 > const & dPhaseComposition_dGlobalCompFraction,
                real64 & value,
                real64 & dValue_dPressure,
                real64 & dValue_dTemperature,
                arraySlice1d< real64, USD4 > const & dValue_dGlobalCompFraction,
                bool useMass ) const;

  virtual void move( LvArray::MemorySpace const space, bool const touch ) override
  {
    PVTFunctionBaseUpdate::move( space, touch );
    m_CO2ViscosityTable.move( space, touch );
  }

protected:

  /// Table with viscosity tabulated as a function (P,T)
  TableFunction::KernelWrapper m_CO2ViscosityTable;

};

class FenghourCO2Viscosity : public PVTFunctionBase
{
public:

  FenghourCO2Viscosity( string_array const & inputParams,
                        string_array const & componentNames,
                        array1d< real64 > const & componentMolarWeight );

  virtual ~FenghourCO2Viscosity() override = default;

  static string catalogName() { return "FenghourCO2Viscosity"; }

  virtual string getCatalogName() const override final { return catalogName(); }

  virtual PVTFunctionType functionType() const override
  {
    return PVTFunctionType::VISCOSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = FenghourCO2ViscosityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

private:

  /// Table with CO2 viscosity tabulated as a function of (P,T)
  TableFunction const * m_CO2ViscosityTable;

};

template< int USD1 >
GEOSX_HOST_DEVICE
void FenghourCO2ViscosityUpdate::compute( real64 const & pressure,
                                          real64 const & temperature,
                                          arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                          real64 & value,
                                          bool useMass ) const
{
  GEOSX_UNUSED_VAR( phaseComposition, useMass )

  real64 const input[2] = { pressure, temperature };
  real64 densityDeriv[2];
  m_CO2ViscosityTable.compute( input, value, densityDeriv );
}

template< int USD1, int USD2, int USD3, int USD4 >
GEOSX_HOST_DEVICE
void FenghourCO2ViscosityUpdate::compute( real64 const & pressure,
                                          real64 const & temperature,
                                          arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                          arraySlice1d< real64 const, USD2 > const & dPhaseComposition_dPressure,
                                          arraySlice1d< real64 const, USD2 > const & dPhaseComposition_dTemperature,
                                          arraySlice2d< real64 const, USD3 > const & dPhaseComposition_dGlobalCompFraction,
                                          real64 & value,
                                          real64 & dValue_dPressure,
                                          real64 & dValue_dTemperature,
                                          arraySlice1d< real64, USD4 > const & dValue_dGlobalCompFraction,
                                          bool useMass ) const
{
  GEOSX_UNUSED_VAR( phaseComposition,
                    dPhaseComposition_dPressure,
                    dPhaseComposition_dTemperature,
                    dPhaseComposition_dGlobalCompFraction,
                    useMass )

  real64 const input[2] = { pressure, temperature };
  real64 densityDeriv[2];
  m_CO2ViscosityTable.compute( input, value, densityDeriv );

  dValue_dPressure = densityDeriv[0];
  dValue_dTemperature = densityDeriv[1];
  LvArray::forValuesInSlice( dValue_dGlobalCompFraction, []( real64 & val ){ val = 0.0; } );
}

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_FENGHOURCO2VISCOSITY_HPP_
