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
 * @file CubicEOSDensity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_CUBICEOSDENSITY_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_CUBICEOSDENSITY_HPP_

#include "FunctionBase.hpp"

#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

template< typename EOS_TYPE >
class CubicEOSDensityUpdate final : public FunctionBaseUpdate
{
public:
  CubicEOSDensityUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                         ComponentProperties const & componentProperties );

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
};

template< typename EOS_TYPE >
class CubicEOSDensity : public FunctionBase
{
public:

  CubicEOSDensity( string const & name,
                   array1d< string > const & componentNames,
                   array1d< real64 > const & componentMolarWeight,
                   ComponentProperties const & componentProperties );

  static string catalogName() { return "CubicEOSDensity"; }

  virtual string getCatalogName() const final { return catalogName(); }

  virtual FunctionType functionType() const override
  {
    return FunctionType::DENSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = CubicEOSDensityUpdate< EOS_TYPE >;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;
};

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

// Implementation
#include "CubicEOSDensityImpl.hpp"

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_CUBICEOSDENSITY_HPP_
