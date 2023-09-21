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
 * @file PhillipsBrineDensity.cpp
 */

#include "constitutive/fluid/multifluid/CO2Brine/functions/PhillipsBrineDensity.hpp"

#include "functions/FunctionManager.hpp"

namespace geos
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
                                        string const & functionName,
                                        FunctionManager & functionManager )
{
  // initialize the (p,T) coordinates
  PTTableCoordinates tableCoords;
  PVTFunctionHelpers::initializePropertyTable( inputParams, tableCoords );

  // initialize salinity
  GEOS_THROW_IF_LT_MSG( inputParams.size(), 9,
                        GEOS_FMT( "{}: insufficient number of model parameters", functionName ),
                        InputError );
  real64 salinity;
  try
  {
    salinity = stod( inputParams[8] );
  }
  catch( std::invalid_argument const & e )
  {
    GEOS_THROW( GEOS_FMT( "{}: invalid model parameter value: {}", functionName, e.what() ), InputError );
  }

  array1d< real64 > densities( tableCoords.nPressures() * tableCoords.nTemperatures() );
  calculateBrineDensity( tableCoords, salinity, densities );

  string const tableName = functionName + "_table";
  if( functionManager.hasGroup< TableFunction >( tableName ) )
  {
    return functionManager.getGroupPointer< TableFunction >( tableName );
  }
  else
  {
    TableFunction * const densityTable = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", tableName ) );
    densityTable->setTableCoordinates( tableCoords.getCoords(),
                                       { units::Pressure, units::TemperatureInC } );
    densityTable->setTableValues( densities, units::Density );
    densityTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    return densityTable;
  }
}

} // namespace

PhillipsBrineDensity::PhillipsBrineDensity( string const & name,
                                            string_array const & inputParams,
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

  m_brineDensityTable = makeDensityTable( inputParams, m_functionName, FunctionManager::getInstance() );
}

PhillipsBrineDensity::KernelWrapper
PhillipsBrineDensity::createKernelWrapper() const
{
  return KernelWrapper( m_componentMolarWeight,
                        *m_brineDensityTable,
                        m_CO2Index,
                        m_waterIndex );
}

void PhillipsBrineDensity::checkTablesParameters( real64 const pressure,
                                                  real64 const temperature ) const
{
  m_brineDensityTable->checkCoord( pressure, 0 );
  m_brineDensityTable->checkCoord( temperature, 1 );
}

REGISTER_CATALOG_ENTRY( PVTFunctionBase, PhillipsBrineDensity, string const &, string_array const &, string_array const &, array1d< real64 > const & )

} // namespace PVTProps

} // namespace constitutive

} // namespace geos
