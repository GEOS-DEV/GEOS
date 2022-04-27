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
 * @file ReactionBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_REACTIONBASE_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_REACTIONBASE_HPP_

#include "dataRepository/ObjectCatalog.hpp"

namespace geosx
{

namespace constitutive
{

namespace chemicalReactions
{

class ReactionBaseUpdate
{
public:

  /**
   * @brief Constructor.
   * @param componentMolarWeight component molar weights
   */
  explicit ReactionBaseUpdate( arrayView1d< real64 const > const & componentMolarWeight )
    :
    m_componentMolarWeight( componentMolarWeight )
  {}

  /**
   * @brief Move the KernelWrapper to the given execution space, optionally touching it.
   * @param space the space to move the KernelWrapper to
   * @param touch whether the KernelWrapper should be touched in the new space or not
   * @note This function exists to enable holding KernelWrapper objects in an ArrayView
   *       and have their contents properly moved between memory spaces.
   */
  virtual void move( LvArray::MemorySpace const space, bool const touch )
  {
    m_componentMolarWeight.move( space, touch );
  }

  void computeLogActCoeff();


protected:

  /// Array storing the component molar weights
  arrayView1d< real64 const > m_componentMolarWeight;

  /// Hard coding the example case - eventually would have to be changed such that it is read from an input file
  localIndex m_numPrimarySpecies = 7;	// Currently not including H2O and O2gas
  localIndex m_numSecSpecies = 11;
  arrayView1d<real64>  m_log10EqConst;	
  m_log10EqConst.resize( m_numSecSpecies );	// Not sure if this is the correct way of allocating the size
  m_log10EqConst[0] = 13.99;
  m_log10EqConst[1] = -6.36;
  m_log10EqConst[2] = 10.33;
  m_log10EqConst[3] = -3.77;
  m_log10EqConst[4] = -1.09;	
  m_log10EqConst[5] = 7.07;
  m_log10EqConst[6] = -2.16;
  m_log10EqConst[7] = 0.67;
  m_log10EqConst[8] = 0.60;
  m_log10EqConst[9] = -2.43;
  m_log10EqConst[10] = -0.82;

  arrayView1d<real64>  m_log10SecActCoeff;	
  m_log10SecActCoeff.resize( m_numSecSpecies );	// Not sure if this is the correct way of allocating the size
// Not sure if this works
  m_log10SecActCoeff[0:10] = 0	//Assume dilute solution for first pass

  arrayView1d<real64>  m_log10PrimaryActCoeff;	
  m_log10PrimaryActCoeff.resize( m_numPrimarySpecies );	// Not sure if this is the correct way of allocating the size
// Not sure if this works
  m_log10PrimaryActCoeff[0:6] = 0	//Assume dilute solution for first pass

// First index: 0 = OH-, 1 = CO2, 2 = CO3-2, 3 = H2CO3, 4 = CaHCO3+, 5 = CaCO3, 6 = CaSO4, 7 = CaCl+, 8 = CaCl2, 9 = MgSO4, 10 = NaSO4-
// Second index: 0 = H+, 1 = HCO3-, 2 = Ca+2, 3 = SO4-2, 4 = Cl-, 5 = Mg+2, 6 = Na+1
  arrayView2d<real64>  m_stochMatrix;	
  m_stochMatrix.resize( m_numSecSpecies, m_numPrimarySpecies );	// Not sure if this is the correct way of allocating the size
// Not sure if this works
  m_stochMatrix[0:10][0:6] = 0;		
  m_stochMatrix[0][0] = -1;		
  m_stochMatrix[1][0] = 1;		
  m_stochMatrix[1][1] = 1;
  m_stochMatrix[2][0] = -1;
  m_stochMatrix[2][1] = 1;
  m_stochMatrix[3][0] = 1;
  m_stochMatrix[3][1] = 1;
  m_stochMatrix[4][1] = 1;
  m_stochMatrix[4][2] = 1;
  m_stochMatrix[5][0] = -1;
  m_stochMatrix[5][1] = 1;
  m_stochMatrix[5][2] = 1;
  m_stochMatrix[6][2] = 1;
  m_stochMatrix[6][3] = 1;
  m_stochMatrix[7][2] = 1;
  m_stochMatrix[7][4] = 1;
  m_stochMatrix[8][2] = 1;
  m_stochMatrix[8][4] = 2;
  m_stochMatrix[9][5] = 1;
  m_stochMatrix[9][3] = 1;
  m_stochMatrix[10][6] = 1;
  m_stochMatrix[10][3] = 1;


};

class ReactionBase
{
public:

  ReactionBase( string const & name,
                  string_array const & componentNames,
                  array1d< real64 > const & componentMolarWeight ):
    m_modelName( name ),
    m_componentNames( componentNames ),
    m_componentMolarWeight( componentMolarWeight )
  {}

  virtual ~ReactionBase() = default;

  using CatalogInterface = dataRepository::CatalogInterface< ReactionBase,
                                                             string const &,
                                                             string_array const &,
                                                             string_array const &,
                                                             string_array const &,
                                                             array1d< real64 > const & >;
  static typename CatalogInterface::CatalogType & getCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  virtual string getCatalogName() const = 0;

  string const & reactionName() const { return m_reactionName; }

protected:

  /// Name the solubility model
  string m_reactionName;

  /// Array storing the name of the components
  string_array m_componentNames;

  /// Array storing the component molar weights
  array1d< real64 > m_componentMolarWeight;

};

} // end namespace chemicalReactions

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_REACTIONBASE_HPP_
