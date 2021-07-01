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

  explicit PVTFunctionBaseUpdate( arrayView1d< real64 const > const & componentMolarWeight )
    : m_componentMolarWeight( componentMolarWeight )
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

protected:

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
