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
 * @file ReactionsBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_REACTIONSBASE_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_REACTIONSBASE_HPP_

#include "dataRepository/ObjectCatalog.hpp"

namespace geosx
{

namespace constitutive
{

namespace chemicalReactions
{

class ReactionsBase
{
public:

  ReactionsBase( string const & name,
                 string_array const & componentNames,
                 array1d< real64 > const & componentMolarWeight ):
    m_modelName( name ),
    m_componentNames( componentNames ),
    m_componentMolarWeight( componentMolarWeight )
  {}

  virtual ~ReactionsBase() = default;

  constexpr static integer maxNumPrimarySpecies = 7; 
  constexpr static integer maxNumSecondarySpecies = 11;

  string const & reactionName() const { return m_reactionName; }

protected:

  /// Name the solubility model
  string m_reactionName;

  /// Array storing the name of the components
  string_array m_componentNames;

  /// Array storing the component molar weights
  array1d< real64 > m_componentMolarWeight;

  array1d< real64 > m_log10EqConst;

  array2d< real64 >  m_stoichMatrix;

  array1d< integer > m_chargePrimary;
  array1d< integer > m_chargeSec;

  array1d< real64>  m_ionSizePrimary;  
  array1d< real64 >  m_ionSizeSec;

  real64 m_DebyeHuckelA;
  real64 m_DebyeHuckelB;
  real64 m_WATEQBDot; 

class KernelWrapper
{
public:

 /**
  * @brief Construct a new Kernel Wrapper object
  * 
  * @param log10EqConst 
  * @param stoichMatrix 
  * @param chargePrimary 
  * @param chargeSec 
  * @param m_ionSizePrimary 
  * @param ionSizeSec 
  * @param DebyeHuckelA 
  * @param DebyeHuckelB 
  * @param WATEQBDot 
  */
  KernelWrapper( arrayView1d< real64 const > const & log10EqConst,
                 arrayView2d< real64 const > const &  stoichMatrix,
                 arrayView1d< integer const > const & chargePrimary,
                 arrayView1d< integer const > const & chargeSec, 
                 arrayView1d< real64 const > const & ionSizePrimary,  
                 arrayView1d< real64 const > const & ionSizeSec,
                 real64 const DebyeHuckelA,
                 real64 const DebyeHuckelB,
                 real64 const WATEQBDot ):
  m_log10EqConst(log10EqConst),
  m_stoichMatrix(stoichMatrix),
  m_chargePrimary(chargePrimary),
  m_chargeSec(chargeSec),
  m_ionSizePrimary(m_ionSizePrimary),  
  m_ionSizeSec(ionSizeSec),
  m_DebyeHuckelA(DebyeHuckelA),
  m_DebyeHuckelB(DebyeHuckelB),
  m_WATEQBDot(WATEQBDot)
  {}

protected:

  /**
   * @brief 
   * 
   * @param temperature 
   * @param ionicStrength 
   * @param log10PrimaryActCoeff 
   * @param dLog10PrimaryActCoeff_dIonicStrength 
   * @param log10SecActCoeff 
   * @param dLog10SecActCoeff_dIonicStrength  
   */
  GEOSX_HOST_DEVICE
  void computeLog10ActCoefBDotModel( real64 const temperature,
                                     real64 const ionicStrength,
                                     arraySlice1d< real64 > & log10PrimaryActCoeff,
                                     arraySlice1d< real64 > & dLog10PrimaryActCoeff_dIonicStrength,
                                     arraySlice1d< real64 > & log10SecActCoeff,
                                     arraySlice1d< real64 > & dLog10SecActCoeff_dIonicStrength ) const;
  /**
   * @brief 
   * 
   * @return  
   */
  GEOSX_HOST_DEVICE
  real64 computeIonicStrength(  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                              arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration ) const;


  /// Hard coding the example case - eventually would have to be changed such that it is read from an input file
  integer m_numPrimarySpecies = 7;	// Currently not including H2O and O2gas
  integer m_numSecSpecies = 11;
  
  arrayView1d< real64 > m_log10EqConst;
  arrayView2d< real64 > m_stoichMatrix;

  arrayView1d< integer > m_chargePrimary;
  arrayView1d< integer > m_chargeSec;

  arrayView1d< real64 >  m_ionSizePrimary;  
  arrayView1d< real64 >  m_ionSizeSec; 

  real64 m_DebyeHuckelA;
  real64 m_DebyeHuckelB;
  real64 m_WATEQBDot; 
};

};

} // end namespace chemicalReactions

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_REACTIONSBASE_HPP_
