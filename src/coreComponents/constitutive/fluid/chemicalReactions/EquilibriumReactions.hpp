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

  class KernelWrapper final : public ReactionBase::KernelWrapper
{
public:

  KernelWrapper( arrayView1d< real64 const > const & log10EqConst,
                 arrayView2d< real64 const > const &  stoichMatrix,
                 arrayView1d< integer const > const & chargePrimary,
                 arrayView1d< integer const > const & chargeSec, 
                 arrayView1d< real64 const > const & ionSizePrimary,  
                 arrayView1d< real64 const > const & ionSizeSec,
                 real64 const DebyeHuckelA,
                 real64 const DebyeHuckelB,
                 real64 const WATEQBDot ): 
  ReactionsBase::KernelWrapper( log10EqConst,
                                stoichMatrix,
                                chargePrimary,
                                chargeSec, 
                                ionSizePrimary,  
                                ionSizeSec,
                                DebyeHuckelA,
                                DebyeHuckelB,
                                WATEQBDot )
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
                        arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration ) const;

private:

  void assembleEquilibriumReactionSystem( real64 const & temperature,
                                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration
                                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                                          arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                                          DenseMatrix & matrix,
                                          DenseVector & rhs ) const;

  void computeSeondarySpeciesConcAndDerivative( real64 const temperature,
                                                arraySlice1d< real64 const > const & primarySpeciesConcentration,
                                                arraySlice1d< real64 > const & secondarySpeciesConcentration,
                                                arraySlice2d< real64 > const & dLog10SecConc_dLog10PrimaryConc ) const;

  void computeTotalConcAndDerivative( real64 const & temperature,
                                      arraySlice1d< real64 const > const & primarySpeciesConcentration,
                                      arraySlice1d< real64 const > const & secondarySpeciesConcentration,
                                      arraySlice2d< real64 const > const & dLog10SecConc_dLog10PrimaryConc,
                                      arraySlice1d< real64 > const & totalConc,
                                      arraySlice2d< real64 > const & dTotalConc_dLog10PrimaryConc) const;                                  


};

};

} // end namespace chemicalReactions

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_REQUILIBRIUMREACTIONS_HPP_
