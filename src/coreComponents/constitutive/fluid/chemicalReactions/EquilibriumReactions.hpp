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

class EquilibriumReactions : public ReactionsBase
{
public:

  EquilibriumReactions( string const & name, integer const numPrimarySpecies, integer const numSecSpecies );

  class KernelWrapper final : public ReactionsBase::KernelWrapper
  {
public:

    using DenseMatrix = stackArray2d< real64, ReactionsBase::maxNumPrimarySpecies * ReactionsBase::maxNumPrimarySpecies >;
    using DenseVector = stackArray1d< real64, ReactionsBase::maxNumPrimarySpecies >;

    KernelWrapper( integer const numPrimarySpecies,
                   integer const numSecondarySpecies,
                   arrayView1d< real64 > const & log10EqConst,
                   arrayView2d< real64 > const & stoichMatrix,
                   arrayView1d< integer > const & chargePrimary,
                   arrayView1d< integer > const & chargeSec,
                   arrayView1d< real64 > const & ionSizePrimary,
                   arrayView1d< real64 > const & ionSizeSec,
                   real64 const DebyeHuckelA,
                   real64 const DebyeHuckelB,
                   real64 const WATEQBDot ):
      ReactionsBase::KernelWrapper( numPrimarySpecies,
                                    numSecondarySpecies,
                                    log10EqConst,
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
    GEOSX_HOST_DEVICE
    void updateConcentrations( real64 const temperature,
                               arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                               arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesContentration,
                               arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration ) const;    
private:

    template<int USD>
    GEOSX_HOST_DEVICE
    void assembleEquilibriumReactionSystem( real64 const temperature,
                                            arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                                            arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                                            arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                                            arrayView2d< real64, USD > const matrix,
                                            arrayView1d< real64 > const rhs ) const;
    GEOSX_HOST_DEVICE
    void computeSeondarySpeciesConcAndDerivative( real64 const temperature,
                                                  arraySlice1d< real64 const > const & log10PrimaryActCoeff,
                                                  arraySlice1d< real64 const > const & dLog10PrimaryActCoeff_dIonicStrength,
                                                  arraySlice1d< real64 const > const & log10SecActCoeff,
                                                  arraySlice1d< real64 const > const & dLog10SecActCoeff_dIonicStrength,
                                                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                                                  arraySlice1d< real64, compflow::USD_COMP - 1  > const & secondarySpeciesConcentration,
                                                  arraySlice2d< real64 > const & dLog10SecConc_dLog10PrimaryConc ) const;
    GEOSX_HOST_DEVICE
    void computeTotalConcAndDerivative( real64 const temperature,
                                        arraySlice1d< real64 const, compflow::USD_COMP - 1  > const & primarySpeciesConcentration,
                                        arraySlice1d< real64 const, compflow::USD_COMP - 1  > const & secondarySpeciesConcentration,
                                        arraySlice2d< real64 const > const & dLog10SecConc_dLog10PrimaryConc,
                                        arraySlice1d< real64 > const & totalConc,
                                        arraySlice2d< real64 > const & dTotalConc_dLog10PrimaryConc ) const;

    GEOSX_HOST_DEVICE
    void updatePrimarySpeciesConcentrations( arrayView1d< real64 const > const solution,
                                             arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration ) const;

    GEOSX_HOST_DEVICE
    void setInitialGuess( arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                          arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration ) const;



    integer m_maxNumIterations = 30;
    real64 m_newtonTol = 1e-6;
  };

/**
 * @brief Create an update kernel wrapper.
 * @return the wrapper
 */
  KernelWrapper createKernelWrapper() const;

};


template< int USD >
GEOSX_HOST_DEVICE
void EquilibriumReactions::KernelWrapper::assembleEquilibriumReactionSystem( real64 const temperature,
                                                                             arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                                                                             arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                                                                             arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                                                                             arrayView2d< real64, USD > const matrix,
                                                                             arrayView1d< real64 > const rhs ) const
{

  stackArray1d< real64, ReactionsBase::maxNumPrimarySpecies > log10PrimaryActCoeff( m_numPrimarySpecies );
  stackArray1d< real64, ReactionsBase::maxNumSecondarySpecies > log10SecActCoeff( m_numSecondarySpecies );
  stackArray2d< real64, ReactionsBase::maxNumSecondarySpecies * ReactionsBase::maxNumPrimarySpecies > dLog10SecConc_dLog10PrimaryConc( m_numSecondarySpecies, m_numPrimarySpecies );
  stackArray1d< real64, ReactionsBase::maxNumPrimarySpecies > totalConcentration( m_numPrimarySpecies );
  stackArray2d< real64, ReactionsBase::maxNumPrimarySpecies * ReactionsBase::maxNumPrimarySpecies > dTotalConc_dLog10PrimaryConc( m_numPrimarySpecies, m_numPrimarySpecies );

  real64 ionicStrength = 0.0;
  stackArray1d< real64, ReactionsBase::maxNumPrimarySpecies > dIonicStrength_dPrimaryConcentration( m_numPrimarySpecies );
  stackArray1d< real64, ReactionsBase::maxNumPrimarySpecies > dLog10PrimaryActCoeff_dIonicStrength( m_numPrimarySpecies );
  stackArray1d< real64, ReactionsBase::maxNumSecondarySpecies > dLog10SecActCoeff_dIonicStrength( m_numSecondarySpecies );

  /// activity coefficients
  computeIonicStrength( primarySpeciesConcentration,
                        secondarySpeciesConcentration,
                        ionicStrength );

  computeLog10ActCoefBDotModel( temperature,
                                ionicStrength,
                                log10PrimaryActCoeff,
                                dLog10PrimaryActCoeff_dIonicStrength,
                                log10SecActCoeff,
                                dLog10SecActCoeff_dIonicStrength );

  computeSeondarySpeciesConcAndDerivative( temperature,
                                           log10PrimaryActCoeff,
                                           dLog10PrimaryActCoeff_dIonicStrength,
                                           log10SecActCoeff,
                                           dLog10SecActCoeff_dIonicStrength,
                                           primarySpeciesConcentration,
                                           secondarySpeciesConcentration,
                                           dLog10SecConc_dLog10PrimaryConc );

  computeTotalConcAndDerivative( temperature,
                                 primarySpeciesConcentration,
                                 secondarySpeciesConcentration,
                                 dLog10SecConc_dLog10PrimaryConc,
                                 totalConcentration,
                                 dTotalConc_dLog10PrimaryConc );

  for( int i=0; i<m_numPrimarySpecies; i++ )
  {
    rhs[i] = 1 - totalConcentration[i] / primarySpeciesTotalConcentration[i];
    rhs[i] = -rhs[i];
    for( int j=0; j<m_numPrimarySpecies; j++ )
    {
      matrix[i][j] = -dTotalConc_dLog10PrimaryConc[i][j] / primarySpeciesTotalConcentration[i];
    }
  }
}

} // end namespace chemicalReactions

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_REQUILIBRIUMREACTIONS_HPP_
