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

#include "ReactionsBase.hpp"

#include "constitutive/fluid/layouts.hpp"

namespace geosx
{

namespace constitutive
{

namespace chemicalReactions
{

class KineticReactions : public ReactionsBase
{
public:

  KineticReactions( string const & name, integer const numPrimarySpecies, integer const numSecSpecies );

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

private:

  m_array1d< real64 > m_reactionRateConstant;

  class KernelWrapper final : public ReactionsBase::KernelWrapper
{
public:

  static constexpr real64 RConst = 8.314;

  KernelWrapper( arrayView1d< real64 const > const & log10EqConst,
                 arrayView2d< real64 const > const &  stoichMatrix,
                 arrayView1d< integer const > const & chargePrimary,
                 arrayView1d< integer const > const & chargeSec, 
                 arrayView1d< real64 const > const & ionSizePrimary,  
                 arrayView1d< real64 const > const & ionSizeSec,
                 real64 const DebyeHuckelA,
                 real64 const DebyeHuckelB,
                 real64 const WATEQBDot,
                 arrayView1d< real64 > const & reactionRateConstant): 
  ReactionsBase::KernelWrapper( log10EqConst,
                                stoichMatrix,
                                chargePrimary,
                                chargeSec, 
                                ionSizePrimary,  
                                ionSizeSec,
                                DebyeHuckelA,
                                DebyeHuckelB,
                                WATEQBDot ),
  m_reactionRateConstant(reactionRateConstant)
  {}

  GEOSX_HOST_DEVICE
  void computeReactionRate( real64 const & temperature,
                            arraySlice1d< geosx::real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                            arraySlice1d< geosx::real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration, ) const;

private:

  arrayView1d< real64 > m_reactionRateConstant;

};

};

} // end namespace chemicalReactions

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_REACTIONBASE_HPP_
