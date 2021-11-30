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
 * @file CO2Enthalpy.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2ENTHALPY_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2ENTHALPY_HPP_

#include "PVTFunctionBase.hpp"

#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"
#include "functions/TableFunction.hpp"

namespace geosx
{

namespace constitutive
{

namespace PVTProps
{

class CO2EnthalpyUpdate final : public PVTFunctionBaseUpdate
{
public:

  CO2EnthalpyUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                     TableFunction const & CO2EnthalpyTable,
                     integer const CO2Index,
                     integer const waterIndex )
    : PVTFunctionBaseUpdate( componentMolarWeight ),
    m_CO2EnthalpyTable( CO2EnthalpyTable.createKernelWrapper() ),
    m_CO2Index( CO2Index ),
    m_waterIndex( waterIndex )
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
    m_CO2EnthalpyTable.move( space, touch );
  }

protected:

  /// Table with CO2 enthalpy tabulated as a function of (P,T)
  TableFunction::KernelWrapper m_CO2EnthalpyTable;

  /// Index of the CO2 component
  integer m_CO2Index;

  /// Index of the water component
  integer m_waterIndex;

};

class CO2Enthalpy : public PVTFunctionBase
{
public:

  CO2Enthalpy( string const & name,
               string_array const & inputParams,
               string_array const & componentNames,
               array1d< real64 > const & componentMolarWeight );

  static string catalogName() { return "CO2Enthalpy"; }

  virtual string getCatalogName() const final { return catalogName(); }

  virtual PVTFunctionType functionType() const override
  {
    return PVTFunctionType::ENTHALPY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = CO2EnthalpyUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  static void calculateCO2Enthalpy( PTTableCoordinates const & tableCoords,
                                    array1d< real64 > const & densities,
                                    array1d< real64 > const & enthalpies );


private:

  /// Table with CO2 enthalpy tabulated as a function of (P,T)
  TableFunction const * m_CO2EnthalpyTable;

  /// Index of the CO2 phase
  integer m_CO2Index;

  /// Index of the water phase
  integer m_waterIndex;

};

template< int USD1 >
GEOSX_HOST_DEVICE
void CO2EnthalpyUpdate::compute( real64 const & pressure,
                                 real64 const & temperature,
                                 arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                 real64 & value,
                                 bool useMass ) const
{
  GEOSX_UNUSED_VAR( pressure );

  real64 const enthalpy = m_CO2EnthalpyTable.compute( &temperature );
  real64 const CO2Enthalpy = m_CO2EnthalpyTable.compute( &temperature );


  //assume there are only CO2 and CO2 here.

  real64 const C = phaseComposition[m_waterIndex];

  real64 const waterMW = m_componentMolarWeight[m_waterIndex];
  real64 const CO2MW = m_componentMolarWeight[m_CO2Index];

  if( useMass )
  {
    real64 const X = C * waterMW / (C * waterMW + (1.0 - C) * CO2MW);
    value = (1.0 - X ) * CO2Enthalpy + X * enthalpy;
  }
  else
  {
    value = (1.0 - C ) * CO2Enthalpy / CO2MW + C * enthalpy / waterMW;
  }
}

template< int USD1, int USD2, int USD3, int USD4 >
GEOSX_HOST_DEVICE
void CO2EnthalpyUpdate::compute( real64 const & pressure,
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
  GEOSX_UNUSED_VAR( pressure );

  real64 const enthalpy = m_CO2EnthalpyTable.compute( &temperature );
  real64 const CO2Enthalpy = m_CO2EnthalpyTable.compute( &temperature );


  //assume there are only CO2 and CO2 here.

  real64 const C = phaseComposition[m_waterIndex];

  real64 const waterMW = m_componentMolarWeight[m_waterIndex];
  real64 const CO2MW = m_componentMolarWeight[m_CO2Index];

  if( useMass )
  {
    real64 const X = C * waterMW / (C * waterMW + (1.0 - C) * CO2MW);
    value = (1.0 - X ) * CO2Enthalpy + X * enthalpy;
  }
  else
  {
    value = (1.0 - C ) * CO2Enthalpy / CO2MW + C * enthalpy / waterMW;
  }
}

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2CO2DENSITY_HPP_
