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

#include "constitutive/fluid/PVTFunctions/BrineCO2DensityFunction.hpp"

#include "managers/Functions/FunctionManager.hpp"
#include "managers/GeosxState.hpp"

namespace geosx
{

using namespace stringutilities;

namespace PVTProps
{

BrineCO2Density::BrineCO2Density( array1d< string > const & inputPara,
                                  array1d< string > const & componentNames,
                                  array1d< real64 > const & componentMolarWeight ):
  PVTFunctionBase( inputPara[1],
                   componentNames,
                   componentMolarWeight )
{
  bool notFound = true;
  for( localIndex i = 0; i < componentNames.size(); ++i )
  {
    if( componentNames[i] == "CO2" || componentNames[i] == "co2" )
    {
      m_CO2Index = i;
      notFound = false;
      break;
    }
  }
  GEOSX_ERROR_IF( notFound, "Component CO2 is not found!" );

  notFound = true;
  for( localIndex i = 0; i < componentNames.size(); ++i )
  {
    if( componentNames[i] =="Water" || componentNames[i] == "water" )
    {
      m_waterIndex = i;
      notFound = false;
      break;
    }
  }
  GEOSX_ERROR_IF( notFound, "Component Water/Brine is not found!" );

  makeTable( inputPara );
}

void BrineCO2Density::makeTable( array1d< string > const & inputPara )
{
  array1d< array1d< real64 > > coordinates;
  coordinates.resize( 2 );

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

  real64 P = PStart;
  while( P <= PEnd )
  {
    coordinates[0].emplace_back( P );
    P += dP;
  }

  real64 T = TStart;
  while( T <= TEnd )
  {
    coordinates[1].emplace_back( T );
    T += dT;
  }

  localIndex const nP = coordinates[0].size();
  localIndex const nT = coordinates[1].size();

  array1d< real64 > values( nP * nT );
  calculateBrineDensity( coordinates, salinity, values );

  FunctionManager & functionManager = getGlobalState().getFunctionManager();
  m_brineDensityTable = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", "brineDensityTable" ) );
  m_brineDensityTable->setTableCoordinates( coordinates );
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

} // namespace geosx
