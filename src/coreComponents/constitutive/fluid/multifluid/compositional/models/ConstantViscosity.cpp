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
 * @file ConstantViscosity.cpp
 */

#include "ConstantViscosity.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

ConstantViscosity::ConstantViscosity( string const & name,
                                      array1d< string > const & componentNames,
                                      array1d< real64 > const & componentMolarWeight,
                                      ComponentProperties const & componentProperties ):
  FunctionBase( name,
                componentNames,
                componentMolarWeight,
                componentProperties )
{}

ConstantViscosity::KernelWrapper
ConstantViscosity::createKernelWrapper() const
{
  return KernelWrapper( m_componentMolarWeight,
                        m_componentProperties );
}

REGISTER_CATALOG_ENTRY( FunctionBase,
                        ConstantViscosity,
                        string const &,
                        string_array const &,
                        array1d< real64 > const &,
                        ComponentProperties const & )

} // end namespace compositional

} // namespace constitutive

} // end namespace geos
