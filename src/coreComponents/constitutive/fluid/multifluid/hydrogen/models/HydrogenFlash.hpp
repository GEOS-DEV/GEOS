/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HydrogenFlash.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_HYDROGEN_MODELS_HYDROGENFLASH_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_HYDROGEN_MODELS_HYDROGENFLASH_HPP_

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"

namespace geos
{

namespace constitutive
{

class HydrogenFlashUpdate
{
public:
  using PhaseProp = MultiFluidBase::PhaseProp;
  using PhaseComp = MultiFluidBase::PhaseComp;

public:
  HydrogenFlashUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                       integer const h2ComponentIndex,
                       integer const h2oComponentIndex,
                       integer const gasPhaseIndex,
                       integer const watPhaseIndex );

  template< int USD1 >
  GEOS_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & compFraction,
                PhaseProp::SliceType const phaseFraction,
                PhaseComp::SliceType const phaseCompFraction ) const;

  void move( LvArray::MemorySpace const space, bool const touch );

protected:
  /// Number of components
  integer m_numComps{0};

  /// Index of the H2 component index
  integer m_h2Index{-1};

  /// Index of the water component index
  integer m_h2oIndex{-1};

  /// Index of the gas phase
  integer m_gasPhaseIndex{-1};

  /// Index of the water phase
  integer m_watPhaseIndex{-1};
};

class HydrogenFlash
{
public:
  HydrogenFlash( string const & name,
                 arrayView1d< real64 > const & componentMolarWeight,
                 integer const h2ComponentIndex,
                 integer const h2oComponentIndex,
                 integer const gasPhaseIndex,
                 integer const watPhaseIndex,
                 bool const printTable );

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = HydrogenFlashUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

private:
  /// Component molecular weights
  arrayView1d< real64 > const & m_componentMolarWeight;

  /// Index of the H2 component index
  integer m_h2Index{-1};

  /// Index of the water component index
  integer m_h2oIndex{-1};

  /// Index of the gas phase
  integer m_gasPhaseIndex{-1};

  /// Index of the water phase
  integer m_watPhaseIndex{-1};
};

template< int USD1 >
GEOS_HOST_DEVICE
void HydrogenFlashUpdate::compute( real64 const & pressure,
                                   real64 const & temperature,
                                   arraySlice1d< real64 const, USD1 > const & compFraction,
                                   PhaseProp::SliceType const phaseFraction,
                                   PhaseComp::SliceType const phaseCompFraction ) const
{
  GEOS_UNUSED_VAR( pressure );
  GEOS_UNUSED_VAR( temperature );

  // Zero out everything to start
  auto const setZero = []( real64 & val ){ val = 0.0; };
  LvArray::forValuesInSlice( phaseFraction.value, setZero );
  LvArray::forValuesInSlice( phaseCompFraction.value, setZero );
  LvArray::forValuesInSlice( phaseFraction.derivs, setZero );
  LvArray::forValuesInSlice( phaseCompFraction.derivs, setZero );

  // 1) Compute phase fractions

  phaseFraction.value[m_gasPhaseIndex] = 0.5;
  phaseFraction.value[m_watPhaseIndex] = 0.5;

  // 2) Compute phase component fractions
  for( integer ic = 0; ic < m_numComps; ++ic )
  {
    phaseCompFraction.value[m_gasPhaseIndex][ic] = compFraction[ic];
    phaseCompFraction.value[m_watPhaseIndex][ic] = compFraction[ic];
  }
}

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_HYDROGEN_MODELS_HYDROGENFLASH_HPP_
