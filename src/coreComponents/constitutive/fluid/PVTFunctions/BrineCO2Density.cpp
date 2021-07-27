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
 * @file BrineCO2Density.cpp
 */

#include "constitutive/fluid/PVTFunctions/BrineCO2Density.hpp"

#include "functions/FunctionManager.hpp"

namespace geosx
{

using namespace stringutilities;

namespace constitutive
{

namespace PVTProps
{

namespace
{

void calculateBrineDensity( PTTableCoordinates const & tableCoords,
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

  localIndex const nPressures = tableCoords.nPressures();
  localIndex const nTemperatures = tableCoords.nTemperatures();

  for( localIndex i = 0; i < nPressures; ++i )
  {
    real64 const P = tableCoords.getPressure( i ) / 1e5;

    for( localIndex j = 0; j < nTemperatures; ++j )
    {
      // see Phillips et al. (1981), equations (4) and (5), pages 14 and 15
      real64 const x = c1 * exp( a1 * salinity )
                       + c2 * exp( a2 * tableCoords.getTemperature( j ) )
                       + c3 * exp( a3 * P );
      densities[j*nPressures+i] = (AA + BB * x + CC * x * x + DD * x * x * x) * 1000.0;
    }
  }
}

TableFunction const * makeDensityTable( string_array const & inputParams,
                                        FunctionManager & functionManager )
{
  // initialize the (p,T) coordinates
  PTTableCoordinates tableCoords;
  PVTFunctionHelpers::initializePropertyTable( inputParams, tableCoords );

  // initialize salinity
  GEOSX_THROW_IF( inputParams.size() < 9,
                  "Invalid property input!",
                  InputError );
  real64 salinity = 0.0;
  try
  {
    salinity = stod( inputParams[8] );
  }
  catch( const std::invalid_argument & e )
  {
    GEOSX_THROW( "Invalid property argument:" + string( e.what() ), InputError );
  }

  array1d< real64 > densities( tableCoords.nPressures() * tableCoords.nTemperatures() );
  calculateBrineDensity( tableCoords, salinity, densities );

  TableFunction * const densityTable = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", "brineDensityTable" ) );
  densityTable->setTableCoordinates( tableCoords.getCoords() );
  densityTable->setTableValues( densities );
  densityTable->reInitializeFunction();
  densityTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );

  return densityTable;
}

} // namespace

BrineCO2Density::BrineCO2Density( string_array const & inputParams,
                                  string_array const & componentNames,
                                  array1d< real64 > const & componentMolarWeight ):
  PVTFunctionBase( inputParams[1],
                   componentNames,
                   componentMolarWeight )
{
  string const expectedCO2ComponentNames[] = { "CO2", "co2" };
  m_CO2Index = PVTFunctionHelpers::findName( componentNames, expectedCO2ComponentNames );

  string const expectedWaterComponentNames[] = { "Water", "water" };
  m_waterIndex = PVTFunctionHelpers::findName( componentNames, expectedWaterComponentNames );

  m_brineDensityTable = makeDensityTable( inputParams, FunctionManager::getInstance() );
}

BrineCO2Density::KernelWrapper BrineCO2Density::createKernelWrapper()
{
  return KernelWrapper( m_componentMolarWeight,
                        *m_brineDensityTable,
                        m_CO2Index,
                        m_waterIndex );
}

REGISTER_CATALOG_ENTRY( PVTFunctionBase, BrineCO2Density, string_array const &, string_array const &, array1d< real64 > const & )

} // namespace PVTProps

} // namespace constitutive

} // namespace geosx
