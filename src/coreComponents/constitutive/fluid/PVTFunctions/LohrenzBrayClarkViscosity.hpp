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
 * @file LohrenzBrayClarkViscosity.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_LOHRENZBRAYCLARKVISCOSITY_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_LOHRENZBRAYCLARKVISCOSITY_HPP_

#include "PVTCompositionalFunctionBase.hpp"

namespace geosx
{

namespace constitutive
{

namespace PVTProps
{

class LohrenzBrayClarkViscosityUpdate final : public PVTCompositionalFunctionBaseUpdate
{
public:

  LohrenzBrayClarkViscosityUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                                   arrayView1d< real64 const > const & componentCriticalPressure,
                                   arrayView1d< real64 const > const & componentCriticalTemperature,
                                   arrayView1d< real64 const > const & componentCriticalVolume,
                                   arrayView1d< real64 const > const & componentAcentricFactor,
                                   arrayView1d< real64 const > const & componentVolumeShift,
                                   arrayView2d< real64 const > const & componentBinaryCoeff )
    : PVTCompositionalFunctionBaseUpdate( componentMolarWeight,
                                          componentCriticalPressure,
                                          componentCriticalTemperature,
                                          componentCriticalVolume,
                                          componentAcentricFactor,
                                          componentVolumeShift,
                                          componentBinaryCoeff )
  {}

  // compute all phases at once, no derivatives version

  template< int USD1, int USD2 >
  GEOSX_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseDensity,
                arraySlice2d< real64 const, USD2 > const & phaseComposition,
                arraySlice1d< real64,       USD1 > const & phaseViscosity ) const
  {
    GEOSX_UNUSED_VAR( pressure, temperature, phaseDensity, phaseComposition );

    auto setConstant = []( real64 & val ){ val = 1.234; };
    LvArray::forValuesInSlice( phaseViscosity, setConstant );
  }

  /*
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
  */

protected:

};


class LohrenzBrayClarkViscosity : public PVTCompositionalFunctionBase
{
public:

  LohrenzBrayClarkViscosity( string const & name,
                             string_array const & componentNames,
                             array1d< real64 > const & componentMolarWeight,
                             array1d< real64 > const & componentCriticalPressure,
                             array1d< real64 > const & componentCriticalTemperature,
                             array1d< real64 > const & componentCriticalVolume,
                             array1d< real64 > const & componentAcentricFactor,
                             array1d< real64 > const & componentVolumeShift,
                             array2d< real64 > const & componentBinaryCoeff );

  virtual ~LohrenzBrayClarkViscosity() override = default;

  static string catalogName() { return "LohrenzBrayClarkViscosity"; }

  virtual string getCatalogName() const override final { return catalogName(); }

  virtual PVTCompositionalFunctionType functionType() const override
  {
    return PVTCompositionalFunctionType::VISCOSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = LohrenzBrayClarkViscosityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

private:

};

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_LOHRENZBRAYCLARKVISCOSITY_HPP_
