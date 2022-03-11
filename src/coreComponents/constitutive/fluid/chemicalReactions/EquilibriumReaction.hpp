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
 * @file EquilibriumReaction.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_REACTIONBASE_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_REACTIONBASE_HPP_

#include "ReactionBase.hpp"

#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"
#include "constitutive/fluid/layouts.hpp"
#include "constitutive/fluid/MultiFluidUtils.hpp"
#include "functions/TableFunction.hpp"

namespace geosx
{

namespace constitutive
{

namespace chemicalReactions
{

constexpr real64 minForDivision = 1e-10;

class EquilibriumReactionUpdate final : public ReactionBaseUpdate
{
public:

  using PhaseProp = MultiFluidVar< real64, 3, multifluid::LAYOUT_PHASE, multifluid::LAYOUT_PHASE_DC >;
  using PhaseComp = MultiFluidVar< real64, 4, multifluid::LAYOUT_PHASE_COMP, multifluid::LAYOUT_PHASE_COMP_DC >;

  EquilibriumReactionUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                             TableFunction const & EquilibriumReactionTable,
                             integer const CO2Index,
                             integer const waterIndex,
                             integer const phaseGasIndex,
                             integer const phaseLiquidIndex )
    : ReactionBaseUpdate( componentMolarWeight ),
    m_EquilibriumReactionTable( EquilibriumReactionTable.createKernelWrapper() ),
    m_CO2Index( CO2Index ),
    m_waterIndex( waterIndex ),
    m_phaseGasIndex( phaseGasIndex ),
    m_phaseLiquidIndex( phaseLiquidIndex )
  {}

  template< int USD1, int USD2, int USD3 >
  GEOSX_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & compFraction,
                arraySlice1d< real64, USD2 > const & phaseFraction,
                arraySlice2d< real64, USD3 > const & phaseCompFraction ) const;

  template< int USD1 >
  GEOSX_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & compFraction,
                PhaseProp::SliceType const phaseFraction,
                PhaseComp::SliceType const phaseCompFraction ) const;

  virtual void move( LvArray::MemorySpace const space, bool const touch ) override
  {
    ReactionBaseUpdate::move( space, touch );
    m_EquilibriumReactionTable.move( space, touch );
  }

protected:

  /// Table with CO2 solubility tabulated as a function (P,T)
  TableFunction::KernelWrapper m_EquilibriumReactionTable;

  /// Index of the CO2 phase
  integer m_CO2Index;

  /// Index of the water phase
  integer m_waterIndex;

  /// Index of the gas phase
  integer m_phaseGasIndex;

  /// Index of the liquid phase
  integer m_phaseLiquidIndex;

};

class EquilibriumReaction : public ReactionBase
{
public:

  EquilibriumReaction( string const & name,
                       string_array const & inputParams,
                       string_array const & phaseNames,
                       string_array const & componentNames,
                       array1d< real64 > const & componentMolarWeight );

  static string catalogName() { return "EquilibriumReaction"; }

  virtual string getCatalogName() const final { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = EquilibriumReactionUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

private:

  m_array1d< real64 >

};

template< int USD1, int USD2, int USD3 >
GEOSX_HOST_DEVICE
inline void
EquilibriumReactionUpdate::compute( real64 const & pressure,
                                    real64 const & temperature,
                                    arraySlice1d< real64 const, USD1 > const & compFraction,
                                    arraySlice1d< real64, USD2 > const & phaseFraction,
                                    arraySlice2d< real64, USD3 > const & phaseCompFraction ) const
{}

template< int USD1 >
GEOSX_HOST_DEVICE
inline void
EquilibriumReactionUpdate::compute( real64 const & pressure,
                                    real64 const & temperature,
                                    arraySlice1d< real64 const, USD1 > const & compFraction,
                                    PhaseProp::SliceType const phaseFraction,
                                    PhaseComp::SliceType const phaseCompFraction ) const
{

}

} // end namespace chemicalReactions

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_REACTIONBASE_HPP_
