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
 * @file WaterDensity.cpp
 */

#include "constitutive/fluid/multifluid/CO2Brine/functions/WaterDensity.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/PureWaterProperties.hpp"
#include "functions/FunctionManager.hpp"

namespace geos
{

using namespace stringutilities;

namespace constitutive
{

namespace PVTProps
{


WaterDensity::WaterDensity( string const & name,
                            integer const logLevel,
                            string_array const & inputParams,
                            string_array const & componentNames,
                            array1d< real64 > const & componentMolarWeight ):
  PVTFunctionBase( name,
                   logLevel,
                   componentNames,
                   componentMolarWeight )
{
  GEOS_UNUSED_VAR( inputParams );
  m_waterDensityTable = PureWaterProperties::makeSaturationDensityTable( m_logLevel, m_functionName, FunctionManager::getInstance() );
}

WaterDensity::KernelWrapper
WaterDensity::createKernelWrapper() const
{
  return KernelWrapper( m_componentMolarWeight,
                        *m_waterDensityTable );
}

REGISTER_CATALOG_ENTRY( PVTFunctionBase, WaterDensity, string const &, integer const, string_array const &, string_array const &, array1d< real64 > const & )

} // namespace PVTProps

} // namespace constitutive

} // namespace geos
