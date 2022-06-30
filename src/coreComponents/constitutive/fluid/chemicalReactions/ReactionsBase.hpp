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

  using CatalogInterface = dataRepository::CatalogInterface< ReactionsBase,
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

  array2d< real64 >  m_stoichMatrix;

  array1d< integer > m_chargePrimary;
  array1d< integer > m_chargeSec;

  array1d< real64>  m_ionSizePrimary;  
  array1d< real64 >  m_ionSizeSec;

class KernelWrapper
{
public:

  /**
   * @brief Constructor.
   * @param componentMolarWeight component molar weights
   */
  explicit KernelWrapper( arrayView1d< real64 const > const & componentMolarWeight ):
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
  /// Hard coding the example case - eventually would have to be changed such that it is read from an input file
  integer m_numPrimarySpecies = 7;	// Currently not including H2O and O2gas
  integer m_numSecSpecies = 11;
  
  arrayView1d<real64> m_log10EqConst;
  arrayView2d<real64> m_stoichMatrix;

  arrayView1d< integer > m_chargePrimary;
  arrayView1d< integer > m_chargeSec;

  arrayView1d<real64>  m_ionSizePrimary;  
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
