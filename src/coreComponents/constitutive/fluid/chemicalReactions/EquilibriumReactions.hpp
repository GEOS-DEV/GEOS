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
 * @file EquilibriumReactions.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_EQUILIBRIUMREACTIONS_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_EQUILIBRIUMREACTIONS_HPP_

#include "ReactionsBase.hpp"

#include "constitutive/fluid/layouts.hpp"

namespace geosx
{

namespace constitutive
{

namespace chemicalReactions
{

class EquilibriumReactions : public ReactionBase
{
public:

  EquilibriumReactions( string const & name );

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

private:

  array1d< real64 >  m_log10EqConst;

  class KernelWrapper final : public ReactionBase::KernelWrapper
{
public:

  KernelWrapper( arrayView1d< real64 const > const & componentMolarWeight ): 
    ReactionBase::KernelWrapper( componentMolarWeight )
  {}
  
  /**
   * @brief Construct a new update Concentrations object
   * 
   * @param temperature 
   * @param totalConc 
   * @param dLog10PrimaryConc_dTotalConc 
   */
  updateConcentrations( real64 const temperature,
                        arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                        arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesContentration, 
                        arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration );

protected:

  array1Viewd< real64 >  m_log10EqConst;

};

};

} // end namespace chemicalReactions

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_REQUILIBRIUMREACTIONS_HPP_
