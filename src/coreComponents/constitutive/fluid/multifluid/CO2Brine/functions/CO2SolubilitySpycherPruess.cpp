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

namespace
{

// Some EOS parameters
struct PengRobinson
{
  static constexpr real64 omegaA = 0.457235529;
  static constexpr real64 omegaB = 0.077796074;

  GEOS_FORCE_INLINE
  static real64
  evaluate( real64 const & omega )
  {
    return ( omega < 0.49 )
      ? 0.37464 + 1.54226 * omega - 0.269920 * omega * omega
      : 0.37960 + 1.48500 * omega - 0.164423 * omega * omega + 0.016666 * omega * omega * omega;
  }

public:
  GEOS_FORCE_INLINE
  static real64
  aMixtureCoefficient( real64 const & ac, real64 const & pr, real64 const & tr )
  {
    real64 const m = evaluate( ac );
    return omegaA * pr / (tr*tr) * pow( 1.0 + m * ( 1.0 - sqrt( tr ) ), 2.0 );
  }

  GEOS_FORCE_INLINE
  static real64
  bMixtureCoefficient( real64 const & ac, real64 const & pr, real64 const & tr )
  {
    GEOS_UNUSED_VAR( ac );
    return omegaB * pr / tr;
  }
};

struct SoaveRedlichKwong
{
private:
  static constexpr real64 omegaA = 0.42748;
  static constexpr real64 omegaB = 0.08664;

  GEOS_FORCE_INLINE
  static real64
  evaluate( real64 const & omega )
  {
    return 0.480 + 1.574 * omega - 0.176 * omega * omega;
  }

public:
  GEOS_FORCE_INLINE
  static real64
  aMixtureCoefficient( real64 const & ac, real64 const & pr, real64 const & tr )
  {
    real64 const m = evaluate( ac );
    return omegaA * pr / (tr*tr) * pow( 1.0 + m * ( 1.0 - sqrt( tr ) ), 2.0 );
  }

  GEOS_FORCE_INLINE
  static real64
  bMixtureCoefficient( real64 const & ac, real64 const & pr, real64 const & tr )
  {
    GEOS_UNUSED_VAR( ac );
    return omegaB * pr / tr;
  }
};

struct RedlichKwong
{
  static constexpr real64 omegaA = 0.42748;
  static constexpr real64 omegaB = 0.08664;

public:
  GEOS_FORCE_INLINE
  static real64
  aMixtureCoefficient( real64 const & ac, real64 const & pr, real64 const & tr )
  {
    GEOS_UNUSED_VAR( ac, pr );
    return omegaA / sqrt( tr );
  }

  GEOS_FORCE_INLINE
  static real64
  bMixtureCoefficient( real64 const & ac, real64 const & pr, real64 const & tr )
  {
    GEOS_UNUSED_VAR( ac, pr, tr );
    return omegaB;
  }
};

}

/**
 * Physical constants
 */
constexpr real64 Pc_CO2 = 7377300.0;    // Critical pressure of CO2 in Pa
constexpr real64 Tc_CO2 = 304.1282;     // Critical temperature of CO2 in K
constexpr real64 Ac_CO2 = 0.22394;      // Accentric factor of CO2
constexpr real64 Pc_H2O = 221.2e5;      // Critical pressure of H2O in Pa
constexpr real64 Tc_H2O = 647.3;        // Critical temperature of H2O in K
constexpr real64 Ac_H2O = 0.344;        // Accentric factor of H2O

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
 * @param[in] T Temperature in C
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
 * @param[in] T Temperature in C
 */
real64 equilibriumConstantH2O( real64 const T )
{
  constexpr real64 c[] = {-2.209, 3.097e-2, -1.098e-4, 2.048e-7};
  real64 const logk0_H2O = c[0] + T * (c[1] + T * (c[2] + T* c[3]));
  return pow( 10.0, logk0_H2O );
}

real64 fugacityCoefficientCO2( real64 const P, real64 const T, real64 const rhoCO2, real64 const salinity )
{
  GEOS_UNUSED_VAR( salinity );
  real64 constexpr molarMassCO2 = 44.01e-3;   // Molar mass of CO2 [kg/mol]
  real64 const V = molarMassCO2/rhoCO2;   // Molar volume [m3/mol]

  // Mixture parameters from the Redlich-Kwong EOS
  real64 const a_CO2 = (7.54e7 - 4.13e4*T);
  real64 constexpr b_CO2 = 27.8;

  real64 constexpr R = 83.1446261815324e-6;   // Universal gas constant [bar.m3/mol.K]

  real64 const lnPhiCO2 = log( V/(V - b_CO2)) + b_CO2/(V - b_CO2)
                          - 2*a_CO2/(R*pow( T, 1.5 )*b_CO2)*log((V + b_CO2)/V )
                          + a_CO2*b_CO2/(R*pow( T, 1.5 )*b_CO2*b_CO2) * (log((V + b_CO2)/V ) - b_CO2/(V + b_CO2))
                          - log( P*V/(R*T));
  return exp( lnPhiCO2 );
}

real64 fugacityCoefficientH2O( real64 const P, real64 const T, real64 const rhoCO2, real64 const salinity )
{
  GEOS_UNUSED_VAR( salinity );
  real64 constexpr molarMassCO2 = 44.01e-3;   // Molar mass of CO2 [kg/mol]
  real64 const V = molarMassCO2/rhoCO2;   // Molar volume [m3/mol]

  // Mixture parameters from the Redlich-Kwong EOS
  real64 const a_CO2 = (7.54e7 - 4.13e4*T);
  real64 constexpr a_CO2_H2O = 7.89e7;
  real64 constexpr b_CO2 = 27.8;
  real64 constexpr b_H2O = 18.18;

  real64 constexpr R = 83.1446261815324e-6;   // Universal gas constant [bar.m3/mol.K]

  real64 const lnPhiH2O = log( V/(V - b_CO2)) + b_H2O/(V - b_CO2)
                          - 2*a_CO2_H2O/(R*pow( T, 1.5 )*b_CO2)*log((V + b_CO2)/V )
                          + a_CO2*b_H2O/(R*pow( T, 1.5 )*b_CO2*b_CO2)
                          *(log((V + b_CO2)/V ) - b_CO2/(V + b_CO2))
                          - log( P*V/(R*T));
  return exp( lnPhiH2O );
}

/**
 * @brief Calculate the parameter A from Eq (11) of Spycher et al. (2003)
 * @param[in] P Pressure in bar
 * @param[in] T Temperature in C
 * @param[in] rhoCO2 CO2 density in kg/m3
 * @param[in] salinity Salinity of the water
 */
real64 computeA( real64 const P, real64 const T, real64 const rhoCO2, real64 const salinity )
{
  real64 constexpr P0 = 1.0;      // Reference pressure is 1 bar
  real64 const deltaP = P - P0;
  real64 constexpr v_av_H2O = 18.1e-6;   // Average partial molar volume of H2O [m3/mol]
  real64 constexpr R = 83.1446261815324e-6;   // Universal gas constant [bar.m3/mol.C]
  real64 const k0_H2O = equilibriumConstantH2O( T ); // K-value for H2O at 1 bar
  real64 const phi_H2O = fugacityCoefficientH2O( P, T, rhoCO2, salinity ); // Fugacity coefficient of H2O for the water-CO2 system
  return k0_H2O/(phi_H2O*P) * exp( deltaP*v_av_H2O/(R*T));
}

/**
 * @brief Calculate the parameter B from Eq (12) of Spycher et al. (2003)
 * @param[in] P Pressure in bar
 * @param[in] T Temperature in C
 * @param[in] rhoCO2 CO2 density in kg/m3
 * @param[in] salinity Salinity of the water
 */
real64 computeB( real64 const P, real64 const T, real64 const rhoCO2, real64 const salinity )
{
  real64 constexpr P0 = 1.0;      // Reference pressure is 1 bar
  real64 const deltaP = P - P0;
  real64 constexpr v_av_CO2 = 32.6e-6;   // Average partial molar volume of CO2 [m3/mol]
  real64 constexpr R = 83.1446261815324e-6;   // Universal gas constant [bar.m3/mol.C]
  real64 const k0_CO2 = equilibriumConstantCO2( T ); // K-value for CO2 at 1 bar
  real64 const phi_CO2 = fugacityCoefficientCO2( P, T, rhoCO2, salinity ); // Fugacity coefficient of CO2 for the water-CO2 system
  return phi_CO2*P/(55.508*k0_CO2) * exp( -(deltaP*v_av_CO2)/(R*T));
}

template< typename LAMBDA >
void calculateSolubility( PTTableCoordinates const & tableCoords,
                          array1d< real64 > const & values,
                          real64 const salinity,
                          real64 const tolerance,
                          LAMBDA && solubilityFunction )
{
  GEOS_UNUSED_VAR( salinity );
  // Calculate the CO2 density
  string const functionName = "SpycherPruessSolubilityCO2Density";
  array1d< real64 > densities( tableCoords.nPressures() * tableCoords.nTemperatures() );
  SpanWagnerCO2Density::calculateCO2Density( functionName,
                                             tolerance,
                                             tableCoords,
                                             densities );

  localIndex const nPressures = tableCoords.nPressures();
  localIndex const nTemperatures = tableCoords.nTemperatures();
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

      // Calculate the solubility
      values[j*nPressures+i] = solubilityFunction( A, B );
    }
  }
}

TableFunction const *
CO2SolubilitySpycherPruess::makeSolubilityTable( string_array const & inputParams,
                                                 string const & functionName,
                                                 FunctionManager & functionManager )
{
  // Initialize the (p,T) coordinates
  PTTableCoordinates tableCoords;
  auto [salinity, tolerance] = readInputParameters( inputParams, functionName, tableCoords );

  array1d< real64 > values( tableCoords.nPressures() * tableCoords.nTemperatures() );
  calculateSolubility( tableCoords,
                       values,
                       salinity,
                       tolerance,
                       []( real64 const A, real64 const B ) -> real64 {
    real64 const x_CO2 = B*(1.0-A)/(1.0 - A*B);
    return x_CO2 / (1.0 - x_CO2);
  } );

  string const tableName = GEOS_FMT( "{}_co2Solubility_table", functionName );
  return makeTable( tableName, tableCoords, std::move( values ), functionManager );
}

TableFunction const *
CO2SolubilitySpycherPruess::makeVapourisationTable( string_array const & inputParams,
                                                    string const & functionName,
                                                    FunctionManager & functionManager )
{
  // Initialize the (p,T) coordinates
  PTTableCoordinates tableCoords;
  PVTFunctionHelpers::initializePropertyTable( inputParams, tableCoords );
  auto [salinity, tolerance] = readInputParameters( inputParams, functionName, tableCoords );

  array1d< real64 > values( tableCoords.nPressures() * tableCoords.nTemperatures() );
  calculateSolubility( tableCoords,
                       values,
                       salinity,
                       tolerance,
                       []( real64 const A, real64 const B ) -> real64 {
    real64 const y_H2O = (1.0-B)/(1.0/A - B);
    return y_H2O / (1.0 - y_H2O);
  } );

  string const tableName = GEOS_FMT( "{}_waterVaporization_table", functionName );
  return makeTable( tableName, tableCoords, std::move( values ), functionManager );
}

} // end namespace PVTProps
} // namespace constitutive
} // end namespace geos
