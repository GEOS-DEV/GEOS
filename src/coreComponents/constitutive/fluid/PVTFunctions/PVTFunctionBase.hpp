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
 * @file PVTFunctionBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_PVTFUNCTIONBASE_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_PVTFUNCTIONBASE_HPP_

#include "codingUtilities/StringUtilities.hpp"
#include "dataRepository/ObjectCatalog.hpp"

namespace geosx
{

namespace constitutive
{

namespace PVTProps
{

enum class PVTFunctionType { UNKNOWN, DENSITY, VISCOSITY };

class PVTFunctionBaseUpdate
{
public:
  PVTFunctionBaseUpdate( arrayView1d< string const > const & componentNames,
                         arrayView1d< real64 const > const & componentMolarWeight )
    :
    m_componentNames( componentNames ),
    m_componentMolarWeight( componentMolarWeight )
  {}

  /// Default virtual destructor
  virtual ~PVTFunctionBaseUpdate() = default;

  /// Default copy constructor
  PVTFunctionBaseUpdate( PVTFunctionBaseUpdate const & ) = default;

  /// Default move constructor
  PVTFunctionBaseUpdate( PVTFunctionBaseUpdate && ) = default;

  /// Deleted copy assignment operator
  PVTFunctionBaseUpdate & operator=( PVTFunctionBaseUpdate const & ) = delete;

  /// Deleted move assignment operator
  PVTFunctionBaseUpdate & operator=( PVTFunctionBaseUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  virtual void compute( real64 const & pressure,
                        real64 const & temperature,
                        arraySlice1d< real64 const > const & phaseComposition,
                        real64 & value,
                        real64 & dValue_dPresure,
                        real64 & dValue_dTemperature,
                        arraySlice1d< real64 > const & dValue_dPhaseComposition,
                        bool useMass = 0 ) const = 0;

protected:

  /// Array storing the name of the components
  arrayView1d< string const > m_componentNames;

  /// Array storing the component molar weights
  arrayView1d< real64 const > m_componentMolarWeight;

};

class PVTFunctionBase
{

public:

  PVTFunctionBase( string const & name,
                   array1d< string > const & componentNames,
                   array1d< real64 > const & componentMolarWeight )
    :
    m_functionName( name ),
    m_componentNames( componentNames ),
    m_componentMolarWeight( componentMolarWeight )
  {}

  virtual ~PVTFunctionBase() = default;

  using CatalogInterface = dataRepository::CatalogInterface< PVTFunctionBase, array1d< string > const &,
                                                             array1d< string > const &,
                                                             array1d< real64 > const & >;
  static typename CatalogInterface::CatalogType & getCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  virtual string getCatalogName() const = 0;

  string const & functionName() const { return m_functionName; }

  virtual PVTFunctionType functionType() const = 0;

protected:

  /// Name of the PVT function
  string m_functionName;

  /// Array storing the name of the components
  array1d< string > m_componentNames;

  /// Array storing the component molar weights
  array1d< real64 > m_componentMolarWeight;

};

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_PVTFUNCTIONBASE_HPP_
