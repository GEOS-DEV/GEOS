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
 * @file NoOpPVTFunction.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_NOOPPVTFUNCTION_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_NOOPPVTFUNCTION_HPP_

#include "PVTFunctionBase.hpp"

namespace geosx
{

namespace constitutive
{

namespace PVTProps
{

class NoOpPVTFunctionUpdate final : public PVTFunctionBaseUpdate
{
public:

  NoOpPVTFunctionUpdate( arrayView1d< real64 const > const & componentMolarWeight )
    : PVTFunctionBaseUpdate( componentMolarWeight )
  {}

  template< int USD1 >
  GEOSX_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                real64 & value,
                bool useMass ) const
  {
    GEOSX_UNUSED_VAR( pressure, temperature, phaseComposition, value, useMass );
  }

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
                bool useMass ) const
  {
    GEOSX_UNUSED_VAR( pressure, temperature,
                      phaseComposition, dPhaseComposition_dPressure, dPhaseComposition_dTemperature, dPhaseComposition_dTemperature, dPhaseComposition_dGlobalCompFraction,
                      value, dValue_dPressure, dValue_dTemperature, dValue_dGlobalCompFraction,
                      useMass );
  }

  virtual void move( LvArray::MemorySpace const space, bool const touch ) override
  {
    PVTFunctionBaseUpdate::move( space, touch );
  }

};


class NoOpPVTFunction : public PVTFunctionBase
{
public:

  NoOpPVTFunction( string const & name,
                   string_array const & inputPara,
                   string_array const & componentNames,
                   array1d< real64 > const & componentMolarWeight )
    : PVTFunctionBase( name,
                       componentNames,
                       componentMolarWeight )
  {
    GEOSX_UNUSED_VAR( inputPara );
  }

  virtual ~NoOpPVTFunction() override = default;

  static string catalogName() { return "NoOpPVTFunction"; }

  virtual string getCatalogName() const override final { return catalogName(); }

  virtual PVTFunctionType functionType() const override
  {
    return PVTFunctionType::UNKNOWN;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = NoOpPVTFunctionUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_componentMolarWeight );
  };

};

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_NOOPPVTFUNCTION_HPP_
