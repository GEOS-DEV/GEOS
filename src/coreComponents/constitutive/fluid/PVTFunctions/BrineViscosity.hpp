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

  BrineViscosityUpdate( arrayView1d< string const > const & componentNames,
                        arrayView1d< real64 const > const & componentMolarWeight,
                        real64 const & coef0,
                        real64 const & coef1 )
    : PVTFunctionBaseUpdate( componentNames,
                             componentMolarWeight ),
    m_coef0( coef0 ),
    m_coef1( coef1 )
  {}

  /// Default copy constructor
  BrineViscosityUpdate( BrineViscosityUpdate const & ) = default;

  /// Default move constructor
  BrineViscosityUpdate( BrineViscosityUpdate && ) = default;

  /// Deleted copy assignment operator
  BrineViscosityUpdate & operator=( BrineViscosityUpdate const & ) = delete;

  /// Deleted move assignment operator
  BrineViscosityUpdate & operator=( BrineViscosityUpdate && ) = delete;

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

  real64 const m_coef0;

  real64 const m_coef1;

};


class BrineViscosity : public PVTFunctionBase
{
public:

  BrineViscosity( string_array const & inputPara,
                  string_array const & componentNames,
                  array1d< real64 > const & componentMolarWeight );

  ~BrineViscosity() override {}

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

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void BrineViscosityUpdate::compute( real64 const & pressure,
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
  GEOSX_UNUSED_VAR( pressure );
  GEOSX_UNUSED_VAR( phaseComposition );
  GEOSX_UNUSED_VAR( dPhaseComposition_dPressure );
  GEOSX_UNUSED_VAR( dPhaseComposition_dTemperature );
  GEOSX_UNUSED_VAR( dPhaseComposition_dGlobalCompFraction );

  GEOSX_UNUSED_VAR( useMass );

  value = m_coef0 + m_coef1 * temperature;
  dValue_dPressure = 0.0;
  dValue_dTemperature = m_coef1;
  for( real64 & val : dValue_dGlobalCompFraction )
  {
    val = 0.0;
  }
}

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINEVISCOSITY_HPP_
