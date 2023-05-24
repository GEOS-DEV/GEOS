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
 * @file EzrokhiBrineDensity.cpp
 */

#include "constitutive/fluid/multifluid/CO2Brine/functions/EzrokhiBrineDensity.hpp"

#include "constitutive/fluid/multifluid/CO2Brine/functions/PureWaterProperties.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/PVTFunctionHelpers.hpp"
#include "functions/FunctionManager.hpp"

namespace geos
{

using namespace stringutilities;

namespace constitutive
{

namespace PVTProps
{

EzrokhiBrineDensity::EzrokhiBrineDensity( string const & name,
                                          string_array const & inputPara,
                                          string_array const & componentNames,
                                          array1d< real64 > const & componentMolarWeight ):
  PVTFunctionBase( name,
                   componentNames,
                   componentMolarWeight )
{
  string const expectedCO2ComponentNames[] = { "CO2", "co2" };
  m_CO2Index = PVTFunctionHelpers::findName( componentNames, expectedCO2ComponentNames, "componentNames" );

  string const expectedWaterComponentNames[] = { "Water", "water" };
  m_waterIndex = PVTFunctionHelpers::findName( componentNames, expectedWaterComponentNames, "componentNames" );

  makeCoefficients( inputPara );
  m_waterSatDensityTable = PureWaterProperties::makeSaturationDensityTable( m_functionName, FunctionManager::getInstance() );
  m_waterSatPressureTable = PureWaterProperties::makeSaturationPressureTable( m_functionName, FunctionManager::getInstance() );
}

void EzrokhiBrineDensity::makeCoefficients( string_array const & inputPara )
{
  // compute brine density following Ezrokhi`s method (referenced in Eclipse TD, Aqueous phase properties)
  // Reference : Zaytsev, I.D. and Aseyev, G.G. Properties of Aqueous Solutions of Electrolytes, Boca Raton, Florida, USA CRC Press (1993).

  m_waterCompressibility = 4.5e-10; // Pa-1
  GEOS_THROW_IF_LT_MSG( inputPara.size(), 5,
                        GEOS_FMT( "{}: insufficient number of model parameters", m_functionName ),
                        InputError );

  try
  {
    // assume CO2 is the only non-water component in the brine
    m_coef0 = stod( inputPara[2] );
    m_coef1 = stod( inputPara[3] );
    m_coef2 = stod( inputPara[4] );
  }
  catch( std::invalid_argument const & e )
  {
    GEOS_THROW( GEOS_FMT( "{}: invalid model parameter value '{}'", m_functionName, e.what() ), InputError );
  }
}

EzrokhiBrineDensity::KernelWrapper
EzrokhiBrineDensity::createKernelWrapper() const
{
  return KernelWrapper( m_componentMolarWeight,
                        *m_waterSatDensityTable,
                        *m_waterSatPressureTable,
                        m_CO2Index,
                        m_waterIndex,
                        m_waterCompressibility,
                        m_coef0,
                        m_coef1,
                        m_coef2 );
}

REGISTER_CATALOG_ENTRY( PVTFunctionBase, EzrokhiBrineDensity, string const &, string_array const &, string_array const &, array1d< real64 > const & )

} // end namespace PVTProps

} // namespace constitutive

} // end namespace geos
