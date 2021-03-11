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

#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"
#include "managers/Functions/FunctionManager.hpp"
#include "managers/GeosxState.hpp"

namespace geosx
{

using namespace stringutilities;

namespace constitutive
{

namespace PVTProps
{

BrineCO2Density::BrineCO2Density( array1d< string > const & inputPara,
                                  array1d< string > const & componentNames,
                                  array1d< real64 > const & componentMolarWeight ):
  PVTFunctionBase( inputPara[1],
                   componentNames,
                   componentMolarWeight )
{
  std::vector< string > const expectedCO2ComponentNames( { "CO2", "co2" } );
  m_CO2Index = PVTFunctionHelpers::findName( componentNames, expectedCO2ComponentNames );
  GEOSX_ERROR_IF( m_CO2Index < 0 || m_CO2Index >= componentNames.size(), "Component CO2 is not found!" );

  std::vector< string > const expectedWaterComponentNames( { "Water", "water" } );
  m_waterIndex = PVTFunctionHelpers::findName( componentNames, expectedWaterComponentNames );
  GEOSX_ERROR_IF( m_waterIndex < 0 || m_waterIndex >= componentNames.size(), "Component Water/Brine is not found!" );

  makeTable( inputPara );
}

void BrineCO2Density::makeTable( array1d< string > const & inputPara )
{
  real64 TStart = -1.0;
  real64 TEnd = -1.0;
  real64 dT = -1.0;
  real64 PStart = -1.0;
  real64 PEnd = -1.0;
  real64 dP = -1.0;
  real64 salinity = 0.0;

  GEOSX_ERROR_IF( inputPara.size() < 9, "Invalid BrineCO2Density input!" );

  try
  {
    PStart = stod( inputPara[2] );
    PEnd = stod( inputPara[3] );
    dP = stod( inputPara[4] );

    TStart = stod( inputPara[5] );
    TEnd = stod( inputPara[6] );
    dT = stod( inputPara[7] );

    salinity = stod( inputPara[8] );
  }
  catch( const std::invalid_argument & e )
  {
    GEOSX_ERROR( "Invalid BrineCO2Density argument:" + string( e.what()) );
  }

  PTTableCoordinates tableCoords;
  for( real64 P = PStart; P <= PEnd; P += dP )
  {
    tableCoords.appendPressure( P );
  }
  for( real64 T = TStart; T <= TEnd; T += dT )
  {
    tableCoords.appendTemperature( T );
  }

  array1d< real64 > values( tableCoords.nPressures() * tableCoords.nTemperatures() );
  calculateBrineDensity( tableCoords.get(), salinity, values );

  FunctionManager & functionManager = getGlobalState().getFunctionManager();
  m_brineDensityTable = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", "brineDensityTable" ) );
  m_brineDensityTable->setTableCoordinates( tableCoords.get() );
  m_brineDensityTable->setTableValues( values );
  m_brineDensityTable->reInitializeFunction();
  m_brineDensityTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
}

void BrineCO2Density::calculateBrineDensity( array1d< array1d< real64 > > const & coordinates,
                                             real64 const & salinity,
                                             array1d< real64 > const & values )
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

  localIndex const numPressures = coordinates[0].size();
  localIndex const numTemperatures = coordinates[1].size();

  for( localIndex i = 0; i < numPressures; ++i )
  {
    real64 const P = coordinates[0][i] / 1e5;

    for( localIndex j = 0; j < numTemperatures; ++j )
    {
      // see Phillips et al. (1981), equations (4) and (5), pages 14 and 15
      real64 const x = c1 * exp( a1 * salinity )
                       + c2 * exp( a2 * coordinates[1][j] )
                       + c3 * exp( a3 * P );
      values[j*numPressures+i] = (AA + BB * x + CC * x * x + DD * x * x * x) * 1000.0;
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

REGISTER_CATALOG_ENTRY( PVTFunctionBase, BrineCO2Density, array1d< string > const &, array1d< string > const &, array1d< real64 > const & )

} // namespace PVTProps

} // namespace constitutive

} // namespace geosx
