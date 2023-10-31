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
 * @file CompositionalDensity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_COMPOSITIONALDENSITY_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_COMPOSITIONALDENSITY_HPP_

#include "FunctionBase.hpp"

#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class CompositionalDensityUpdate final : public FunctionBaseUpdate
{
public:
  CompositionalDensityUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                              ComponentProperties const & componentProperties );

  template< int USD1 >
  GEOS_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                real64 & molarDensity,
                real64 & massDensity,
                bool useMass ) const;

  template< int USD1, int USD2, int USD3 >
  GEOS_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                arraySlice2d< real64 const, USD2 > const & dPhaseComposition,
                real64 & molarDensity,
                real64 & massDensity,
                arraySlice1d< real64, USD3 > const & dMolarDensity,
                arraySlice1d< real64, USD3 > const & dMassDensity,
                bool useMass ) const;
};

class CompositionalDensity : public FunctionBase
{
public:

  CompositionalDensity( string const & name,
                        array1d< string > const & componentNames,
                        array1d< real64 > const & componentMolarWeight,
                        ComponentProperties const & componentProperties );

  static string catalogName() { return "CompositionalDensity"; }

  virtual FunctionType functionType() const override
  {
    return FunctionType::DENSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = CompositionalDensityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;
};

template< int USD1 >
GEOS_HOST_DEVICE
void CompositionalDensityUpdate::compute( real64 const & pressure,
                                          real64 const & temperature,
                                          arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                          real64 & molarDensity,
                                          real64 & massDensity,
                                          bool useMass ) const
{
  GEOS_UNUSED_VAR( pressure, temperature, useMass );
  GEOS_UNUSED_VAR( phaseComposition );

  massDensity = 1000.0;
  molarDensity = massDensity/40.0;
}

template< int USD1, int USD2, int USD3 >
GEOS_HOST_DEVICE
void CompositionalDensityUpdate::compute( real64 const & pressure,
                                          real64 const & temperature,
                                          arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                          arraySlice2d< real64 const, USD2 > const & dPhaseComposition,
                                          real64 & molarDensity,
                                          real64 & massDensity,
                                          arraySlice1d< real64, USD3 > const & dMolarDensity,
                                          arraySlice1d< real64, USD3 > const & dMassDensity,
                                          bool useMass ) const
{
  GEOS_UNUSED_VAR( pressure, temperature, useMass );
  GEOS_UNUSED_VAR( phaseComposition, dPhaseComposition );

  massDensity = 1000.0;
  molarDensity = massDensity/40.0;

  LvArray::forValuesInSlice( dMolarDensity, setZero );
  LvArray::forValuesInSlice( dMassDensity, setZero );
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_COMPOSITIONALDENSITY_HPP_
