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

#include "constitutive/fluid/multifluid/compositional/functions/CompositionalProperties.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/CubicEOSPhaseModel.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class CompositionalDensityUpdate final : public FunctionBaseUpdate
{
public:
  CompositionalDensityUpdate( arrayView1d< real64 const > const & volumeShift,
                              integer const eosType );

  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  void compute( ComponentProperties::KernelWrapper const & componentProperties,
                EquationOfState::KernelWrapper const & equationOfState,
                real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                real64 & molarDensity,
                arraySlice1d< real64, USD2 > const & dMolarDensity,
                real64 & massDensity,
                arraySlice1d< real64, USD2 > const & dMassDensity,
                bool useMass ) const;

private:
  template< integer USD >
  GEOS_HOST_DEVICE
  void computeCompressibilityFactor( integer const numComps,
                                     real64 const & pressure,
                                     real64 const & temperature,
                                     arraySlice1d< real64 const, USD > const & composition,
                                     ComponentProperties::KernelWrapper const & componentProperties,
                                     real64 & compressibilityFactor,
                                     arraySlice1d< real64 > const & compressibilityFactorDerivs ) const;

private:
  arrayView1d< real64 const > m_componentDimensionalVolumeShift;
  integer const m_eosType;
};

class CompositionalDensity : public FunctionBase
{
public:
  CompositionalDensity( string const & name,
                        ComponentProperties const & componentProperties,
                        EquationOfState const & equationOfState,
                        integer const phaseIndex );

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

private:
  void calculateDimensionalVolumeShift( ComponentProperties const & componentProperties,
                                        EquationOfState const & equationOfState,
                                        arraySlice1d< real64 > componentDimensionalVolumeShift );
private:
  integer const m_phaseIndex;
  array1d< real64 > m_componentDimensionalVolumeShift;
};

template< integer USD1, integer USD2 >
GEOS_HOST_DEVICE
void CompositionalDensityUpdate::compute( ComponentProperties::KernelWrapper const & componentProperties,
                                          EquationOfState::KernelWrapper const & equationOfState,
                                          real64 const & pressure,
                                          real64 const & temperature,
                                          arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                          real64 & molarDensity,
                                          arraySlice1d< real64, USD2 > const & dMolarDensity,
                                          real64 & massDensity,
                                          arraySlice1d< real64, USD2 > const & dMassDensity,
                                          bool useMass ) const
{
  GEOS_UNUSED_VAR( useMass );
  GEOS_UNUSED_VAR( equationOfState );

  integer const numComps = componentProperties.m_componentMolarWeight.size();
  integer const numDofs = 2 + numComps;

  real64 compressibilityFactor = 0.0;
  stackArray1d< real64, 2+MultiFluidConstants::MAX_NUM_COMPONENTS > tempDerivs( numDofs );

  computeCompressibilityFactor( numComps,
                                pressure,
                                temperature,
                                phaseComposition,
                                componentProperties,
                                compressibilityFactor,
                                tempDerivs.toSlice() );

  CompositionalProperties::computeMolarDensity( numComps,
                                                pressure,
                                                temperature,
                                                phaseComposition,
                                                m_componentDimensionalVolumeShift.toSliceConst(),
                                                compressibilityFactor,
                                                tempDerivs.toSlice(),
                                                molarDensity,
                                                dMolarDensity );

  CompositionalProperties::computeMassDensity( numComps,
                                               phaseComposition,
                                               componentProperties.m_componentMolarWeight.toSliceConst(),
                                               molarDensity,
                                               dMolarDensity.toSliceConst(),
                                               massDensity,
                                               dMassDensity );
}

template< integer USD >
GEOS_HOST_DEVICE
void CompositionalDensityUpdate::computeCompressibilityFactor( integer const numComps,
                                                               real64 const & pressure,
                                                               real64 const & temperature,
                                                               arraySlice1d< real64 const, USD > const & composition,
                                                               ComponentProperties::KernelWrapper const & componentProperties,
                                                               real64 & compressibilityFactor,
                                                               arraySlice1d< real64 > const & compressibilityFactorDerivs ) const
{
  if( EquationOfState::equals( m_eosType, EquationOfStateType::PengRobinson ))
  {
    CubicEOSPhaseModel< PengRobinsonEOS >::
    computeCompressibilityFactor( numComps,
                                  pressure,
                                  temperature,
                                  composition,
                                  componentProperties,
                                  compressibilityFactor,
                                  compressibilityFactorDerivs );
  }
  else if( EquationOfState::equals( m_eosType, EquationOfStateType::SoaveRedlichKwong ))
  {
    CubicEOSPhaseModel< SoaveRedlichKwongEOS >::
    computeCompressibilityFactor( numComps,
                                  pressure,
                                  temperature,
                                  composition,
                                  componentProperties,
                                  compressibilityFactor,
                                  compressibilityFactorDerivs );
  }
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_COMPOSITIONALDENSITY_HPP_
