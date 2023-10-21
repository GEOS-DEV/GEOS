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
 * @file CubicEOSDensity.cpp
 */

#include "CubicEOSDensity.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

template< typename EOS_TYPE >
CubicEOSDensityUpdate< EOS_TYPE >::CubicEOSDensityUpdate(
  arrayView1d< real64 const > const & componentMolarWeight,
  ComponentProperties const & componentProperties ):
  FunctionBaseUpdate( componentMolarWeight,
                      componentProperties )
{}

template< typename EOS_TYPE >
CubicEOSDensity< EOS_TYPE >::CubicEOSDensity( string const & name,
                                              array1d< string > const & componentNames,
                                              array1d< real64 > const & componentMolarWeight,
                                              ComponentProperties const & componentProperties ):
  FunctionBase( name,
                componentNames,
                componentMolarWeight,
                componentProperties )
{}

template< typename EOS_TYPE >
typename CubicEOSDensity< EOS_TYPE >::KernelWrapper
CubicEOSDensity< EOS_TYPE >::createKernelWrapper() const
{
  return KernelWrapper( m_componentMolarWeight,
                        m_componentProperties );
}

// Explicit instantiation of the model template.
template class CubicEOSDensity< PengRobinsonEOS >;
template class CubicEOSDensity< SoaveRedlichKwongEOS >;

REGISTER_CATALOG_ENTRY( FunctionBase, CubicEOSDensityPR,
                        string const &,
                        string_array const &,
                        array1d< real64 > const &,
                        ComponentProperties const & )

REGISTER_CATALOG_ENTRY( FunctionBase, CubicEOSDensitySRK,
                        string const &,
                        string_array const &,
                        array1d< real64 > const &,
                        ComponentProperties const & )

} // namespace PVTProps

} // namespace constitutive

} // namespace geos
