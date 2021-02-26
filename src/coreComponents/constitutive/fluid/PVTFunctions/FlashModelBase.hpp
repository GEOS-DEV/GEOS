/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FlashModelBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_NEWFLASHMODELBASE_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_NEWFLASHMODELBASE_HPP_

#include "codingUtilities/StringUtilities.hpp"
#include "dataRepository/ObjectCatalog.hpp"

namespace geosx
{

namespace PVTProps
{

class FlashModelBaseUpdate
{
public:

  FlashModelBaseUpdate( arrayView1d< string const > const & componentNames,
                        arrayView1d< real64 const > const & componentMolarWeight )
    :
    m_componentNames( componentNames ),
    m_componentMolarWeight( componentMolarWeight )
  {}

  /// Default copy constructor
  FlashModelBaseUpdate( FlashModelBaseUpdate const & ) = default;

  /// Default move constructor
  FlashModelBaseUpdate( FlashModelBaseUpdate && ) = default;

  /// Deleted copy assignment operator
  FlashModelBaseUpdate & operator=( FlashModelBaseUpdate const & ) = delete;

  /// Deleted move assignment operator
  FlashModelBaseUpdate & operator=( FlashModelBaseUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  virtual void compute( real64 const & pressure,
                        real64 const & temperature,
                        arraySlice1d< real64 const > const & compFraction,
                        arraySlice1d< real64 > const & phaseFraction,
                        arraySlice1d< real64 > const & dPhaseFraction_dPressure,
                        arraySlice1d< real64 > const & dPhaseFraction_dTemperature,
                        arraySlice2d< real64 > const & dPhaseFraction_dCompFraction,
                        arraySlice2d< real64 > const & phaseCompFraction,
                        arraySlice2d< real64 > const & dPhaseCompFraction_dPressure,
                        arraySlice2d< real64 > const & dPhaseCompFraction_dTemperature,
                        arraySlice3d< real64 > const & dPhaseCompFraction_dCompFraction ) const = 0;

protected:

  /// Array storing the name of the components
  arrayView1d< string const > m_componentNames;

  /// Array storing the component molar weights
  arrayView1d< real64 const > m_componentMolarWeight;

};


class FlashModelBase
{
public:

  FlashModelBase( string const & name,
                  array1d< string > const & componentNames,
                  array1d< real64 > const & componentMolarWeight ):
    m_modelName( name ),
    m_componentNames( componentNames ),
    m_componentMolarWeight( componentMolarWeight )
  {}

  virtual ~FlashModelBase() = default;

  using CatalogInterface = dataRepository::CatalogInterface< FlashModelBase, array1d< string > const &,
                                                             array1d< string > const &,
                                                             array1d< string > const &,
                                                             array1d< real64 > const & >;
  static typename CatalogInterface::CatalogType & getCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  virtual string getCatalogName() const = 0;

  string const & flashModelName() const { return m_modelName; }

protected:

  /// Name the solubility model
  string m_modelName;

  /// Array storing the name of the components
  array1d< string > m_componentNames;

  /// Array storing the component molar weights
  array1d< real64 > m_componentMolarWeight;

};

} // end namespace PVTProps

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_NEWFLASHMODELBASE_HPP_
