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
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/CompositionalProperties.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

template< typename EOS_TYPE >
class CompositionalDensityUpdate final : public FunctionBaseUpdate
{
public:
  explicit CompositionalDensityUpdate( arrayView1d< real64 const > const & volumeShift )
    : m_componentDimensionalVolumeShift( volumeShift )
  {}

  GEOS_HOST_DEVICE
  void compute( ComponentProperties::KernelWrapper const & componentProperties,
                real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const > const & phaseComposition,
                arraySlice2d< real64 const > const & dPhaseComposition,
                real64 & molarDensity,
                arraySlice1d< real64 > const & dMolarDensity,
                real64 & massDensity,
                arraySlice1d< real64 > const & dMassDensity,
                bool useMass ) const;

private:
  arrayView1d< real64 const > m_componentDimensionalVolumeShift;
};

template< typename EOS_TYPE >
class CompositionalDensity : public FunctionBase
{
public:
  CompositionalDensity( string const & name,
                        ComponentProperties const & componentProperties )
    : FunctionBase( name, componentProperties )
  {
    // Calculate the dimensional volume shift
    m_componentDimensionalVolumeShift.resize( componentProperties.getNumberOfComponents());
    EOS_TYPE::calculateDimensionalVolumeShift( componentProperties,
                                               m_componentDimensionalVolumeShift );
  }

  static string catalogName() { return "CompositionalDensity"; }

  virtual FunctionType functionType() const override
  {
    return FunctionType::DENSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = CompositionalDensityUpdate< EOS_TYPE >;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_componentDimensionalVolumeShift );
  }

private:
  array1d< real64 > m_componentDimensionalVolumeShift;
};

template< typename EOS_TYPE >
GEOS_HOST_DEVICE
void CompositionalDensityUpdate< EOS_TYPE >::
compute( ComponentProperties::KernelWrapper const & componentProperties,
         real64 const & pressure,
         real64 const & temperature,
         arraySlice1d< real64 const > const & phaseComposition,
         arraySlice2d< real64 const > const & dPhaseComposition,
         real64 & molarDensity,
         arraySlice1d< real64 > const & dMolarDensity,
         real64 & massDensity,
         arraySlice1d< real64 > const & dMassDensity,
         bool useMass ) const
{
  GEOS_UNUSED_VAR( useMass );
  GEOS_UNUSED_VAR( dPhaseComposition );

  integer const numComps = componentProperties.m_componentMolarWeight.size();
  integer const numDofs = 2 + numComps;

  real64 compressibilityFactor = 0.0;
  stackArray1d< real64, 2+MultiFluidConstants::MAX_NUM_COMPONENTS > compressibilityFactorDerivs( numDofs );

  EOS_TYPE::computeCompressibilityFactor( numComps,
                                          pressure,
                                          temperature,
                                          phaseComposition,
                                          componentProperties,
                                          compressibilityFactor,
                                          compressibilityFactorDerivs );

  CompositionalProperties::computeMolarDensity( numComps,
                                                pressure,
                                                temperature,
                                                phaseComposition,
                                                m_componentDimensionalVolumeShift,
                                                compressibilityFactor,
                                                compressibilityFactorDerivs,
                                                molarDensity,
                                                dMolarDensity );

  CompositionalProperties::computeMassDensity( numComps,
                                               phaseComposition,
                                               componentProperties.m_componentMolarWeight,
                                               molarDensity,
                                               dMolarDensity,
                                               massDensity,
                                               dMassDensity );
  std::cout << molarDensity << " " << massDensity << "\n";
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_COMPOSITIONALDENSITY_HPP_
