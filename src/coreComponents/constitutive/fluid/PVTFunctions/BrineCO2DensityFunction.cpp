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
 * @file BrineCO2DensityFunction.cpp
 */

#include "constitutive/fluid/PVTFunctions/BrineCO2DensityFunction.hpp"

#include "managers/Functions/FunctionManager.hpp"
#include "managers/GeosxState.hpp"

namespace geosx
{

using namespace stringutilities;

namespace constitutive
{

namespace PVTProps
{

BrineCO2Density::BrineCO2Density( string_array const & inputPara,
                                  string_array const & componentNames,
                                  array1d< real64 > const & componentMolarWeight ):
  PVTFunctionBase( inputPara[1],
                   componentNames,
                   componentMolarWeight )
{
  char const * expectedCO2ComponentNames[] = { "CO2", "co2" };
  m_CO2Index = PVTFunctionHelpers::findName( componentNames, expectedCO2ComponentNames );
  GEOSX_ERROR_IF( m_CO2Index < 0 || m_CO2Index >= componentNames.size(), "Component CO2 is not found!" );

  char const * expectedWaterComponentNames[] = { "Water", "water" };
  m_waterIndex = PVTFunctionHelpers::findName( componentNames, expectedWaterComponentNames );
  GEOSX_ERROR_IF( m_waterIndex < 0 || m_waterIndex >= componentNames.size(), "Component Water/Brine is not found!" );

  makeTable( inputPara );
}

void BrineCO2Density::makeTable( string_array const & inputPara )
{
  PTTableCoordinates tableCoords;
  real64 const salinity = PVTFunctionHelpers::initializePropertyTableWithSalinity( inputPara, tableCoords );

  array1d< real64 > densities( tableCoords.nPressures() * tableCoords.nTemperatures() );
  calculateBrineDensity( tableCoords, salinity, densities );

  FunctionManager & functionManager = getGlobalState().getFunctionManager();
  m_brineDensityTable = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", "brineDensityTable" ) );
  m_brineDensityTable->setTableCoordinates( tableCoords.getCoords() );
  m_brineDensityTable->setTableValues( densities );
  m_brineDensityTable->reInitializeFunction();
  m_brineDensityTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
}

void BrineCO2Density::calculateBrineDensity( PTTableCoordinates const & tableCoords,
                                             real64 const & salinity,
                                             array1d< real64 > const & densities )
{
  // these coefficients come from Phillips et al. (1981), equations (4) and (5), pages 14 and 15
  constexpr real64 c1 = -9.9595;
  constexpr real64 c2 = 7.0845;
  constexpr real64 c3 = 3.9093;

  constexpr real64 a1 = -0.004539;
  constexpr real64 a2 = -0.0001638;
  constexpr real64 a3 = 0.00002551;

  constexpr real64 AA = -3.033405;
  constexpr real64 BB = 10.128163;
  constexpr real64 CC = -8.750567;
  constexpr real64 DD = 2.663107;

  for( localIndex i = 0; i < tableCoords.nPressures(); ++i )
  {
    real64 const P = tableCoords.getPressure( i ) / 1e5;

    for( localIndex j = 0; j < tableCoords.nTemperatures(); ++j )
    {
      // see Phillips et al. (1981), equations (4) and (5), pages 14 and 15
      real64 const x = c1 * exp( a1 * salinity )
                       + c2 * exp( a2 * tableCoords.getTemperature( j ) )
                       + c3 * exp( a3 * P );
      densities[j*tableCoords.nPressures()+i] = (AA + BB * x + CC * x * x + DD * x * x * x) * 1000.0;
    }
  }
}

BrineCO2Density::KernelWrapper BrineCO2Density::createKernelWrapper()
{
  return KernelWrapper( m_componentNames,
                        m_componentMolarWeight,
                        m_brineDensityTable,
                        m_CO2Index,
                        m_waterIndex );
}

REGISTER_CATALOG_ENTRY( PVTFunctionBase, BrineCO2Density, string_array const &, string_array const &, array1d< real64 > const & )

} // namespace PVTProps

} // namespace constitutive

} // namespace geosx
