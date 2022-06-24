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

#include "constitutive/fluid/layouts.hpp"

namespace geosx
{

namespace constitutive
{

namespace chemicalReactions
{

class EquilibriumReactionUpdate final : public ReactionBaseUpdate
{
public:

  EquilibriumReactionUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                             TableFunction const & EquilibriumReactionTable )
    : ReactionBaseUpdate( componentMolarWeight )
  {}

  computeReactionTerm( real64 const & temperature,
                       arraySlice1d< real64 > const & totalConc,
                       arraySlice2d< real64 > const & dTotalConc_dLog10PrimaryConc );

protected:

  array1Viewd< real64 >  m_log10EqConst;

};

class EquilibriumReaction : public ReactionBase
{
public:

  EquilibriumReaction( string const & name,
                       string_array const & inputParams,
                       string_array const & phaseNames,
                       string_array const & componentNames,
                       array1d< real64 > const & componentMolarWeight );

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = EquilibriumReactionUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

private:

  array1d< real64 >  m_log10EqConst;

};


// This computes the value of the function we are trying to solve (mass balance equation)
template< typename PHASE1, typename PHASE2, typename FLASH >
GEOSX_HOST_DEVICE inline void
EquilibriumReactionUpdate::computeReactionTerm( real64 const & temperature,
                                                arraySlice1d< real64 > const & totalConc,
                                                arraySlice2d< real64 > const & dTotalConc_dLog10PrimaryConc )

{ 
  stackArray1d< ReactionBase::maxNumPrimarySpecies > log10PrimaryConc(m_numPrimarySpecies), log10PrimaryActCoeff(m_numPrimarySpecies);
  stackArray1d< ReactionBase::maxNumSecondarySpecies > log10SecConc(m_numSecSpecies);
  stackArray2d< ReactionBase::maxNumSecondarySpecies > dLog10SecConc_dLog10PrimaryConc(m_numPrimarySpecies, m_numSecSpecies);

  computeLog10SecConcAndDerivative( temperature, log10PrimaryConc, log10SecConc, dLog10SecConc_dLog10PrimaryConc );
  
  computeTotalConcAndDerivative( temperature, log10PrimaryConc, log10SecConc, dLog10SecConc_dLog10PrimaryConc, totalConc, dTotalConc_dLog10PrimaryConc );
  
  //Matteo: I assume we want to solve this.
  stackArray1d<ReactionBase::maxNumPrimarySpecies> funValue (m_numPrimarySpecies);
  for(int i=0; i<funValue.size(); i++)
  {
    funValue[i] = totalConc[i] - inputTotalConc[i];
  }
}


} // end namespace chemicalReactions

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_REACTIONBASE_HPP_
