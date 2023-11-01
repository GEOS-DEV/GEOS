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
 * @file FunctionBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_FUNCTIONBASE_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_FUNCTIONBASE_HPP_

#include "dataRepository/ObjectCatalog.hpp"

#include "ComponentProperties.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

enum class FunctionType : integer
{
  UNKNOWN,
  FLASH,
  DENSITY,
  VISCOSITY,
  ENTHALPY,
  INTERNAL_ENERGY
};

class FunctionBaseUpdate
{
public:

  explicit FunctionBaseUpdate( ComponentProperties const & componentProperties ):
    m_componentProperties ( componentProperties.createKernelWrapper() )
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
    m_componentProperties.move( space, touch );
  }

protected:
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void setZero( real64 & val ){ val = 0.0; }

protected:

  /// Array storing the component molar weights
  arrayView1d< real64 const > m_componentMolarWeight;

  // Standard compositional input
  ComponentProperties::KernelWrapper m_componentProperties;
};

class FunctionBase
{

public:

  FunctionBase( string const & name,
                ComponentProperties const & componentProperties ):
    m_functionName( name ),
    m_componentProperties ( componentProperties )
  {}

  virtual ~FunctionBase() = default;

  string const & functionName() const { return m_functionName; }

  virtual FunctionType functionType() const = 0;

protected:
  /// Name of the PVT function
  string m_functionName;

  // Standard compositional input
  ComponentProperties const & m_componentProperties;
};

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_FUNCTIONBASE_HPP_
