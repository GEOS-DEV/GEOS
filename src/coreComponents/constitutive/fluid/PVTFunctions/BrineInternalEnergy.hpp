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
 * @file BrineInternalEnergy.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINEINTERNALENERGY_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINEINTERNALENERGY_HPP_

#include "PVTFunctionBase.hpp"

#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"
#include "functions/TableFunction.hpp"

namespace geosx
{

namespace constitutive
{

namespace PVTProps
{

class BrineInternalEnergyUpdate final : public PVTFunctionBaseUpdate
{
public:

  BrineInternalEnergyUpdate( arrayView1d< real64 const > const & componentMolarWeight )
    : PVTFunctionBaseUpdate( componentMolarWeight )
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
  }

protected:

};

class BrineInternalEnergy : public PVTFunctionBase
{
public:

  BrineInternalEnergy( string const & name,
                       string_array const & inputParams,
                       string_array const & componentNames,
                       array1d< real64 > const & componentMolarWeight )
    : PVTFunctionBase( name,
                       componentNames,
                       componentMolarWeight )
  {}

  static string catalogName() { return "BrineInternalEnergy"; }

  virtual string getCatalogName() const final { return catalogName(); }

  virtual PVTFunctionType functionType() const override
  {
    return PVTFunctionType::INTERNAL_ENERGY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = BrineInternalEnergyUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

private:

};

template< int USD1 >
GEOSX_HOST_DEVICE
void BrineInternalEnergyUpdate::compute( real64 const & pressure,
                                         real64 const & temperature,
                                         arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                         real64 & value,
                                         bool useMass ) const
{
  GEOSX_UNUSED_VAR( phaseComposition );
  GEOSX_UNUSED_VAR( useMass );

  value =  0.001 * pressure + 1.0 * temperature;
}

template< int USD1, int USD2, int USD3, int USD4 >
GEOSX_HOST_DEVICE
void BrineInternalEnergyUpdate::compute( real64 const & pressure,
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
  GEOSX_UNUSED_VAR( phaseComposition );
  GEOSX_UNUSED_VAR( dPhaseComposition_dPressure );
  GEOSX_UNUSED_VAR( dPhaseComposition_dTemperature );
  GEOSX_UNUSED_VAR( dPhaseComposition_dGlobalCompFraction );
  GEOSX_UNUSED_VAR( useMass );

  value =  0.001 * pressure + 1.0 * temperature;
  dValue_dPressure = 0.001;
  dValue_dTemperature = 1.0;

  LvArray::forValuesInSlice( dValue_dGlobalCompFraction, []( real64 & val ){ val = 0.0; } );
}

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINEINTERNALENERGY_HPP_
