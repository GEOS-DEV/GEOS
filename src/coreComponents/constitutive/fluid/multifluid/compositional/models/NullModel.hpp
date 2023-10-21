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
 * @file NullModel.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NULLMODEL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NULLMODEL_HPP_

#include "FunctionBase.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class NullModelUpdate final : public FunctionBaseUpdate
{
public:

  NullModelUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                   ComponentProperties const & componentProperties ):
    FunctionBaseUpdate( componentMolarWeight,
                        componentProperties )
  {}

  template< int USD1 >
  GEOS_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                real64 & value,
                bool useMass ) const
  {
    GEOS_UNUSED_VAR( pressure, temperature, phaseComposition, value, useMass );
  }

  template< int USD1, int USD2, int USD3 >
  GEOS_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                arraySlice2d< real64 const, USD2 > const & dPhaseComposition,
                real64 & value,
                arraySlice1d< real64, USD3 > const & dValue,
                bool useMass ) const
  {
    GEOS_UNUSED_VAR( pressure, temperature,
                     phaseComposition, dPhaseComposition,
                     value, dValue,
                     useMass );
  }
};

class NullModel : public FunctionBase
{
public:

  NullModel( string const & name,
             string_array const & componentNames,
             array1d< real64 > const & componentMolarWeight,
             ComponentProperties const & componentProperties ):
    FunctionBase( name,
                  componentNames,
                  componentMolarWeight,
                  componentProperties )
  {}

  virtual ~NullModel() override = default;

  static string catalogName() { return "NullCompositionalPVTModel"; }

  virtual string getCatalogName() const override final { return catalogName(); }

  virtual FunctionType functionType() const override
  {
    return FunctionType::UNKNOWN;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = NullModelUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_componentMolarWeight,
                          m_componentProperties );
  };

};

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NULLMODEL_HPP_
