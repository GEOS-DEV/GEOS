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
 * @file BrineViscosity.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINEVISCOSITY_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINEVISCOSITY_HPP_

#include "PVTFunctionBase.hpp"

namespace geosx
{

namespace constitutive
{

namespace PVTProps
{

class BrineViscosityUpdate final : public PVTFunctionBaseUpdate
{
public:

  BrineViscosityUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                        real64 const coef0,
                        real64 const coef1 )
    : PVTFunctionBaseUpdate( componentMolarWeight ),
    m_coef0( coef0 ),
    m_coef1( coef1 )
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

  real64 m_coef0;

  real64 m_coef1;

};


class BrineViscosity : public PVTFunctionBase
{
public:

  BrineViscosity( string_array const & inputPara,
                  string_array const & componentNames,
                  array1d< real64 > const & componentMolarWeight );

  virtual ~BrineViscosity() override = default;

  static string catalogName() { return "BrineViscosity"; }

  virtual string getCatalogName() const override final { return catalogName(); }

  virtual PVTFunctionType functionType() const override
  {
    return PVTFunctionType::VISCOSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = BrineViscosityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

private:

  void makeCoefficients( string_array const & inputPara );

  real64 m_coef0;

  real64 m_coef1;

};

template< int USD1 >
GEOSX_HOST_DEVICE
void BrineViscosityUpdate::compute( real64 const & pressure,
                                    real64 const & temperature,
                                    arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                    real64 & value,
                                    bool useMass ) const
{
  GEOSX_UNUSED_VAR( pressure, phaseComposition, useMass )
  value = m_coef0 + m_coef1 * temperature;
}

template< int USD1, int USD2, int USD3, int USD4 >
GEOSX_HOST_DEVICE
void BrineViscosityUpdate::compute( real64 const & pressure,
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
  GEOSX_UNUSED_VAR( pressure,
                    phaseComposition,
                    dPhaseComposition_dPressure,
                    dPhaseComposition_dTemperature,
                    dPhaseComposition_dGlobalCompFraction,
                    useMass )

  value = m_coef0 + m_coef1 * temperature;
  compute( pressure, temperature, phaseComposition, value, useMass );
  dValue_dPressure = 0.0;
  dValue_dTemperature = m_coef1;
  LvArray::forValuesInSlice( dValue_dGlobalCompFraction, []( real64 & val ){ val = 0.0; } );
}

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINEVISCOSITY_HPP_
