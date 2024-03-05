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
 * @file LohrenzBrayClarkViscosity.cpp
 */

#include "LohrenzBrayClarkViscosity.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

LohrenzBrayClarkViscosityUpdate::LohrenzBrayClarkViscosityUpdate( MixingType const mixing_type )
  : m_mixing_type( mixing_type )
{}

LohrenzBrayClarkViscosity::LohrenzBrayClarkViscosity( string const & name,
                                                      ComponentProperties const & componentProperties ):
  FunctionBase( name, componentProperties )
{}

LohrenzBrayClarkViscosity::KernelWrapper
LohrenzBrayClarkViscosity::createKernelWrapper() const
{
  return KernelWrapper( m_mixing_type );
}

} // end namespace compositional

} // namespace constitutive

} // end namespace geos
