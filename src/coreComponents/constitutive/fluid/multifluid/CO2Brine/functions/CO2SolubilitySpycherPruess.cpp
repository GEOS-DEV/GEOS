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
 * @file CO2SolubilitySpycherPruess.cpp
 */

#include "CO2SolubilitySpycherPruess.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/PVTFunctionHelpers.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/SpanWagnerCO2Density.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/CubicEOSPhaseModel.hpp"
#include "functions/TableFunction.hpp"
#include "functions/FunctionManager.hpp"
#include "common/Units.hpp"

namespace geos
{
namespace constitutive
{
namespace PVTProps
{

/**
 * Physical constants
 */
static constexpr real64 P_ref = 1.0;                        // Reference pressure is 1 bar
static constexpr real64 Pa_2_bar = 1.0e-5;                  // Conversion factor from Pa to Bar
static constexpr real64 R = 10.0*constants::gasConstant;    // Universal gas constant in [bar.cm3/mol.K]
static constexpr real64 molarMassCO2 = 44.01e-3;            // Molar mass of CO2 [kg/mol]
static constexpr real64 molarMassH2O = 18.01e-3;            // Molar mass of H2O [kg/mol]
static constexpr real64 v_av_H2O = 18.1;                    // Average partial molar volume of H2O [cm3/mol]
static constexpr real64 v_av_CO2 = 32.6;                    // Average partial molar volume of CO2 [cm3/mol]

TableFunction const * makeTable( string const & tableName,
                                 PTTableCoordinates const & tableCoords,
                                 array1d< real64 > && values,
                                 FunctionManager & functionManager )
{
  if( functionManager.hasGroup< TableFunction >( tableName ) )
  {
    return functionManager.getGroupPointer< TableFunction >( tableName );
  }
  else
  {
    TableFunction * table = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", tableName ) );
    table->setTableCoordinates( tableCoords.getCoords(),
                                { units::Pressure, units::TemperatureInC } );
    table->setTableValues( values, units::Solubility );
    table->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    return table;
  }
}

// Read the input parameters, populate tableCoords and return the salinity and tolerance
std::pair< real64 const, real64 const > readInputParameters( string_array const & inputParams,
                                                             string const & functionName,
                                                             PTTableCoordinates & tableCoords )
{
  PVTFunctionHelpers::initializePropertyTable( inputParams, tableCoords );

  // Initialize salinity and tolerance
  GEOS_THROW_IF_LT_MSG( inputParams.size(), 9,
                        GEOS_FMT( "{}: insufficient number of model parameters", functionName ),
                        InputError );

  real64 tolerance = 1e-9;
  real64 salinity = 0.0;
  try
  {
    salinity = stod( inputParams[8] );
    if( inputParams.size() >= 10 )
    {
      tolerance = stod( inputParams[9] );
    }
  }
  catch( const std::invalid_argument & e )
  {
    GEOS_THROW( GEOS_FMT( "{}: invalid model parameter value: {}", functionName, e.what() ), InputError );
  }
  return {salinity, tolerance};
}

/**
 * @brief Calculate the CO2 equilibrium constant from Spycher et al. (2003)
 * @details The correlation parameters are given in Table 2 of Spycher et al. (2003)
 * @param[in] T Temperature [C]
 */
real64 equilibriumConstantCO2( real64 const T )
{
  constexpr real64 c[] = {1.189, 1.304e-2, -5.446e-5};
  real64 const logk0_CO2 = c[0] + T * (c[1] + T * c[2]);
  return pow( 10.0, logk0_CO2 );
}

/**
 * @brief Calculate the H2O equilibrium constant from Spycher et al. (2003)
 * @details The correlation parameters are given in Table 2 of Spycher et al. (2003)
 * @param[in] T Temperature [C]
 */
real64 equilibriumConstantH2O( real64 const T )
{
  constexpr real64 c[] = {-2.209, 3.097e-2, -1.098e-4, 2.048e-7};
  real64 const logk0_H2O = c[0] + T * (c[1] + T * (c[2] + T* c[3]));
  return pow( 10.0, logk0_H2O );
}

/**
 * @brief Calculate the fugacity coefficient of CO2 from the RK-EOS
 * @param[in] P Pressure [bar]
 * @param[in] T Temperature [K]
 * @param[in] rhoCO2 Density of CO2 [kg/m3]
 * @param[in] salinity Salinity of water
 */
real64 fugacityCoefficientCO2( real64 const P, real64 const T, real64 const rhoCO2, real64 const salinity )
{
  GEOS_UNUSED_VAR( salinity );
  real64 const V = 1.0e6*molarMassCO2/rhoCO2;   // Molar volume [cm3/mol]

  // Mixture parameters from the tuned Redlich-Kwong EOS
  // These values are given in Table 1 of Spycher et al. (2003)
  real64 const a_CO2 = (7.54e7 - 4.13e4*T);     // [bar.cm3.K^(1/2)/mol^2]
  real64 constexpr b_CO2 = 27.8;                // [cm3/mol]

  real64 const lnPhiCO2 = log( V/(V - b_CO2)) + b_CO2/(V - b_CO2)
                          - 2*a_CO2/(R*pow( T, 1.5 )*b_CO2)*log((V + b_CO2)/V )
                          + a_CO2*b_CO2/(R*pow( T, 1.5 )*b_CO2*b_CO2) * (log((V + b_CO2)/V ) - b_CO2/(V + b_CO2))
                          - log( P*V/(R*T));
  return exp( lnPhiCO2 );
}

/**
 * @brief Calculate the fugacity coefficient of CO2 from the RK-EOS
 * @param[in] P Pressure [bar]
 * @param[in] T Temperature [K]
 * @param[in] rhoCO2 Density of CO2 [kg/m3]
 * @param[in] salinity Salinity of water
 */
real64 fugacityCoefficientH2O( real64 const P, real64 const T, real64 const rhoCO2, real64 const salinity )
{
  GEOS_UNUSED_VAR( salinity );
  real64 const V = 1.0e6*molarMassCO2/rhoCO2;   // Molar volume [cm3/mol]

  // Mixture parameters from the tuned Redlich-Kwong EOS
  // These values are given in Table 1 of Spycher et al. (2003)
  real64 const a_CO2 = (7.54e7 - 4.13e4*T);     // [bar.cm3.K^(1/2)/mol^2]
  real64 constexpr a_CO2_H2O = 7.89e7;          // [bar.cm3.K^(1/2)/mol^2]
  real64 constexpr b_CO2 = 27.8;                // [cm3/mol]
  real64 constexpr b_H2O = 18.18;               // [cm3/mol]

  real64 const lnPhiH2O = log( V/(V - b_CO2) ) + b_H2O/(V - b_CO2)
                          - 2.0 * a_CO2_H2O * log( (V + b_CO2)/V ) / ( R*pow( T, 1.5 )*b_CO2 )
                          + a_CO2 * b_H2O / ( R*pow( T, 1.5 )*b_CO2*b_CO2 )*( log( (V + b_CO2)/V ) - b_CO2/(V + b_CO2) )
                          - log( P*V/(R*T) );
  return exp( lnPhiH2O );
}

/**
 * @brief Calculate the parameter A from Eq (11) of Spycher et al. (2003)
 * @param[in] P Pressure [Pa]
 * @param[in] T Temperature [C]
 * @param[in] rhoCO2 CO2 density [kg/m3]
 * @param[in] salinity Salinity of the water
 */
real64 computeA( real64 const P, real64 const T, real64 const rhoCO2, real64 const salinity )
{
  real64 const P_in_bar = P * Pa_2_bar;
  real64 const deltaP = P_in_bar - P_ref;
  real64 const TinK = units::convertCToK( T );
  real64 const k0_H2O = equilibriumConstantH2O( T ); // K-value for H2O at 1 bar
  real64 const phi_H2O = fugacityCoefficientH2O( P_in_bar, TinK, rhoCO2, salinity ); // Fugacity coefficient of H2O for the water-CO2 system
  real64 const A = k0_H2O/(phi_H2O*P_in_bar) * exp( deltaP*v_av_H2O/(R*TinK));
  return A;
}

/**
 * @brief Calculate the parameter B from Eq (12) of Spycher et al. (2003)
 * @param[in] P Pressure [Pa]
 * @param[in] T Temperature [C]
 * @param[in] rhoCO2 CO2 density [kg/m3]
 * @param[in] salinity Salinity of the water
 */
real64 computeB( real64 const P, real64 const T, real64 const rhoCO2, real64 const salinity )
{
  real64 const P_in_bar = P * Pa_2_bar;
  real64 const deltaP = P_in_bar - P_ref;
  real64 const TinK = units::convertCToK( T );
  real64 const k0_CO2 = equilibriumConstantCO2( T ); // K-value for CO2 at 1 bar
  real64 const phi_CO2 = fugacityCoefficientCO2( P_in_bar, TinK, rhoCO2, salinity ); // Fugacity coefficient of CO2 for the water-CO2 system
  real64 const B = phi_CO2*P_in_bar/(55.508*k0_CO2) * exp( -(deltaP*v_av_CO2)/(R*TinK) );
  return B;
}

std::pair< TableFunction const *, TableFunction const * >
CO2SolubilitySpycherPruess::makeSolubilityTables( string_array const & inputParams,
                                                  string const & functionName,
                                                  FunctionManager & functionManager )
{
  // Initialize the (p,T) coordinates
  PTTableCoordinates tableCoords;
  auto [salinity, tolerance] = readInputParameters( inputParams, functionName, tableCoords );

  localIndex const nPressures = tableCoords.nPressures();
  localIndex const nTemperatures = tableCoords.nTemperatures();

  // Calculate the CO2 density
  array1d< real64 > densities( nPressures*nTemperatures );
  SpanWagnerCO2Density::calculateCO2Density( GEOS_FMT( "{}_co2_density", functionName ),
                                             tolerance,
                                             tableCoords,
                                             densities );

  array1d< real64 > co2Values( nPressures*nTemperatures );
  array1d< real64 > h2oValues( nPressures*nTemperatures );
  for( localIndex i = 0; i < nPressures; ++i )
  {
    real64 const P = tableCoords.getPressure( i );

    for( localIndex j = 0; j < nTemperatures; ++j )
    {
      real64 const T = tableCoords.getTemperature( j );

      // Get the CO2 density
      real64 const rhoCO2 = densities[j*nPressures+i];

      // Calculate A and B
      real64 const A = computeA( P, T, rhoCO2, salinity );
      real64 const B = computeB( P, T, rhoCO2, salinity );

      // Calculate the mole fractions
      // Eqns (13) and (14) from Spycher et al. (2003)
      real64 const y_H2O = (1.0 - B)/(1.0/A - B);
      real64 const x_CO2 = B*(1.0 - y_H2O);

      // Calculate the solubility
      co2Values[j*nPressures+i] = x_CO2/((1.0 - x_CO2)*molarMassH2O);
      h2oValues[j*nPressures+i] = y_H2O/((1.0 - y_H2O)*molarMassCO2);
    }
  }

  string const co2TableName = GEOS_FMT( "{}_co2Solubility_table", functionName );
  TableFunction const * co2SolubilityTable = makeTable( co2TableName, tableCoords, std::move( co2Values ), functionManager );
  string const h2oTableName = GEOS_FMT( "{}_waterVaporization_table", functionName );
  TableFunction const * h2oSolubilityTable = makeTable( h2oTableName, tableCoords, std::move( h2oValues ), functionManager );
  return {co2SolubilityTable, h2oSolubilityTable};
}

} // end namespace PVTProps
} // namespace constitutive
} // end namespace geos
