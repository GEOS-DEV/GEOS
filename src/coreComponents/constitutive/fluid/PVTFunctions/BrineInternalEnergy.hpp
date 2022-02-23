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

#include "constitutive/fluid/layouts.hpp"
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

  template< int USD1, int USD2, int USD3 >
  GEOSX_HOST_DEVICE
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
  {
    // reserve for future: more accurate internal energy model should probably have some parameters
    GEOSX_UNUSED_VAR ( inputParams );
  }

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
  GEOSX_UNUSED_VAR( phaseComposition, useMass );

  value =  0.001 * pressure + 1.0 * temperature;
}

template< int USD1, int USD2, int USD3 >
GEOSX_HOST_DEVICE
void BrineInternalEnergyUpdate::compute( real64 const & pressure,
                                         real64 const & temperature,
                                         arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                         arraySlice2d< real64 const, USD2 > const & dPhaseComposition,
                                         real64 & value,
                                         arraySlice1d< real64, USD3 > const & dValue,
                                         bool useMass ) const
{
  GEOSX_UNUSED_VAR( phaseComposition, dPhaseComposition, useMass );

  using Deriv = multifluid::DerivativeOffset;

  value =  0.001 * pressure + 1.0 * temperature;
  LvArray::forValuesInSlice( dValue, []( real64 & val ){ val = 0.0; } );
  dValue[Deriv::dP] = 0.001;
  dValue[Deriv::dT] = 1.0;

}

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINEINTERNALENERGY_HPP_
