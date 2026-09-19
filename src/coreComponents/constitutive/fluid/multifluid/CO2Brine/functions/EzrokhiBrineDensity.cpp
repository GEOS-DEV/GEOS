/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
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
                                          BrineFluidParameters const & brineFluidParameters,
                                          string_array const & componentNames,
                                          array1d< real64 > const & componentMolarWeight,
                                          TableFunction::OutputOptions const pvtOutputOpts ):
  PVTFunctionBase( name,
                   componentNames,
                   componentMolarWeight ),
  m_waterCompressibility( brineFluidParameters.m_waterCompressibility )
{
  string const expectedCO2ComponentNames[] = { "CO2", "co2" };
  m_CO2Index = PVTFunctionHelpers::findName( componentNames, expectedCO2ComponentNames, "componentNames" );

  string const expectedWaterComponentNames[] = { "Water", "water" };
  m_waterIndex = PVTFunctionHelpers::findName( componentNames, expectedWaterComponentNames, "componentNames" );

  makeCoefficients( brineFluidParameters.m_ezrokhiDensityCoefficients );
  m_waterSatDensityTable = PureWaterProperties::makeSaturationDensityTable( m_functionName, FunctionManager::getInstance() );
  m_waterSatPressureTable = PureWaterProperties::makeSaturationPressureTable( m_functionName, FunctionManager::getInstance() );

  m_waterSatPressureTable->outputTableData( pvtOutputOpts );
  m_waterSatDensityTable->outputTableData( pvtOutputOpts );
}

void EzrokhiBrineDensity::makeCoefficients( arrayView1d< real64 const > const & coefficients )
{
  // compute brine density following Ezrokhi`s method
  // Reference : Zaytsev, I.D. and Aseyev, G.G. Properties of Aqueous Solutions of Electrolytes, Boca Raton, Florida, USA CRC Press (1993).
  m_coef0 = coefficients[0];
  m_coef1 = coefficients[1];
  m_coef2 = coefficients[2];
}

void EzrokhiBrineDensity::checkTablesParameters( real64 const GEOS_UNUSED_PARAM( pressure ),
                                                 real64 const temperature ) const
{
  m_waterSatDensityTable->checkCoord( temperature, 0 );
  m_waterSatPressureTable->checkCoord( temperature, 0 );
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

} // end namespace PVTProps

} // namespace constitutive

} // end namespace geos
