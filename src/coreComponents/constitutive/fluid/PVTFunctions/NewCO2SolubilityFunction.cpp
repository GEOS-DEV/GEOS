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
 * @file CO2Solubility.cpp
 */

#include "constitutive/fluid/PVTFunctions/NewCO2SolubilityFunction.hpp"

#include "managers/Functions/FunctionManager.hpp"
#include "managers/GeosxState.hpp"

namespace geosx
{

using namespace stringutilities;


namespace PVTProps
{

constexpr real64 T_K_f  = 273.15;
constexpr real64 P_Pa_f = 1e+5;
constexpr real64 P_c    = 73.773 * P_Pa_f;
constexpr real64 T_c    = 304.1282;
constexpr real64 Rgas   = 8.314467;
constexpr real64 V_c    = Rgas*T_c/P_c;

// these coefficients are in Table (A1) of Duan and Sun (2003)
constexpr real64 acoef[] =
{ 8.99288497e-2, -4.94783127e-1, 4.77922245e-2, 1.03808883e-2, -2.82516861e-2, 9.49887563e-2, 5.20600880e-4,
  -2.93540971e-4, -1.77265112e-3, -2.51101973e-5, 8.93353441e-5, 7.88998563e-5, -1.66727022e-2, 1.398, 2.96e-2 };

namespace detail
{

real64 ff( real64 const & T, real64 const & P, real64 const & V_r )
{
  // reduced pressure
  real64 const P_r = P*P_Pa_f/P_c;
  // reduced temperature
  real64 const T_r = (T_K_f+T)/T_c;

  // CO2 equation of state
  // see equation (A1) in Duan and Sun (2003)
  real64 const f_Z = 1.0
                     + ( acoef[0] + acoef[1]/(T_r * T_r) + acoef[2]/(T_r * T_r * T_r) )/V_r
                     + ( acoef[3] + acoef[4]/(T_r * T_r) + acoef[5]/(T_r * T_r * T_r) )/(V_r*V_r)
                     + ( acoef[6] + acoef[7]/(T_r * T_r) + acoef[8]/(T_r * T_r * T_r) )/(V_r*V_r*V_r*V_r)
                     + ( acoef[9] + acoef[10]/(T_r * T_r) + acoef[11]/(T_r * T_r * T_r) )/(V_r*V_r*V_r*V_r*V_r)
                     + acoef[12]/(T_r * T_r * T_r)/(V_r * V_r) * (acoef[13] + acoef[14]/(V_r * V_r)) * exp( -acoef[14]/(V_r * V_r)) - P_r * V_r / T_r;

  return f_Z;
}

real64 PWater( real64 const & T )
{
  // these coefficients are defined in Table (B1) of Duan and Sun (2003)
  constexpr real64 ccoef[] = { -38.640844, 5.8948420, 59.876516, 26.654627, 10.637097 };

  // H2O critical pressure (bars)
  real64 const P_c_w = 220.85;
  // H2O critical temperature (K)
  real64 const T_c_w = 647.29;
  real64 const tt = ( (T+T_K_f)-T_c_w )/T_c_w;
  // Empirical model for water pressure of equation (B1) of Duan and Sun (2003)
  real64 const x = (P_c_w*(T+T_K_f)/T_c_w)
                   * (1
                      + ccoef[0]*pow( -tt, 1.9 )
                      + ccoef[1]*tt
                      + ccoef[2]*tt*tt
                      + ccoef[3]*tt*tt*tt
                      + ccoef[4]*tt*tt*tt*tt);

  return x;
}

real64 logF( real64 const & T, real64 const & P, real64 const & V_r )
{
  // reduced pressure
  real64 const P_r = P*P_Pa_f/P_c;
  // reduced temperature
  real64 const T_r = (T_K_f+T)/T_c;
  real64 const Z   = P_r * V_r/T_r;

  // fugacity coefficient of CO2, equation (A6) of Duan and Sun (2003)
  real64 const log_f = Z - 1 - log( Z ) +
                       ( acoef[0] + acoef[1]/T_r/T_r + acoef[2]/T_r/T_r/T_r )/V_r
                       + ( acoef[3] + acoef[4]/T_r/T_r + acoef[5]/T_r/T_r/T_r )/2.0/V_r/V_r
                       + ( acoef[6] + acoef[7]/T_r/T_r + acoef[8]/T_r/T_r/T_r )/4.0/V_r/V_r/V_r/V_r
                       + ( acoef[9] + acoef[10]/T_r/T_r + acoef[11]/T_r/T_r/T_r )/5.0/V_r/V_r/V_r/V_r/V_r
                       + acoef[12]/2.0/T_r/T_r/T_r/acoef[14] * ( acoef[13] + 1.0 - (acoef[13] + 1.0 + acoef[14]/V_r/V_r) * exp( -acoef[14]/V_r/V_r ) );

  return log_f;
}

real64 Par( real64 const & T, real64 const & P, real64 const * cc )
{
  // "equation for the parameters", see equation (7) of Duan and Sun (2003)
  real64 x = cc[0]
             + cc[1]*T
             + cc[2]/T
             + cc[3]*T*T
             + cc[4]/(630.0-T)
             + cc[5]*P
             + cc[6]*P *log( T )
             + cc[7]*P/T
             + cc[8]*P/(630.0-T)
             + cc[9]*P*P/(630.0-T)/(630.0-T)
             + cc[10]*T *log( P );

  return x;
}

}

CO2Solubility::CO2Solubility( array1d< string > const & inputPara,
                              array1d< string > const & phaseNames,
                              array1d< string > const & componentNames,
                              array1d< real64 > const & componentMolarWeight ):
  FlashModelBase( inputPara[1],
                  componentNames,
                  componentMolarWeight )
{
  GEOSX_ERROR_IF( phaseNames.size() != 2, "The CO2Solubility model is a two-phase model" );
  GEOSX_ERROR_IF( componentNames.size() != 2, "The CO2Solubility model is a two-component model" );

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
    if( componentNames[i] == "Water" || componentNames[i] == "water" )
    {
      m_waterIndex = i;
      notFound = false;
      break;
    }
  }
  GEOSX_ERROR_IF( notFound, "Component Water/Brine is not found!" );

  notFound = true;
  for( localIndex i = 0; i < phaseNames.size(); ++i )
  {
    if( phaseNames[i] == "CO2" || phaseNames[i] == "co2" ||
        phaseNames[i] == "gas" || phaseNames[i] == "Gas" )
    {
      m_phaseGasIndex = i;
      notFound = false;
      break;
    }
  }
  GEOSX_ERROR_IF( notFound, "Phase co2/gas is not found!" );

  notFound = true;
  for( localIndex i = 0; i < phaseNames.size(); ++i )
  {
    if( streq( phaseNames[i], "Water" ) || streq( phaseNames[i], "water" ) ||
        streq( phaseNames[i], "Liquid" ) || streq( phaseNames[i], "liquid" ) )
    {
      m_phaseLiquidIndex = i;
      notFound = false;
      break;
    }
  }
  GEOSX_ERROR_IF( notFound, "Phase water/liquid is not found!" );

  makeTable( inputPara );
}

void CO2Solubility::makeTable( array1d< string > const & inputPara )
{
  array1d< array1d< real64 > > coordinates;
  coordinates.resize( 2 );

  real64 TStart = -1.0;
  real64 TEnd = -1.0;
  real64 dT = -1.0;
  real64 PStart = -1.0;
  real64 PEnd = -1.0;
  real64 dP = -1.0;
  real64 salinity = 0;

  GEOSX_ERROR_IF( inputPara.size() < 9, "Invalid CO2Solubility input!" );

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
    GEOSX_ERROR( "Invalid CO2Solubility argument:" + string( e.what()) );
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
  calculateCO2Solubility( coordinates, salinity, values );

  FunctionManager & functionManager = getGlobalState().getFunctionManager();
  m_CO2SolubilityTable = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", "CO2SolubilityTable" ) );
  m_CO2SolubilityTable->setTableCoordinates( coordinates );
  m_CO2SolubilityTable->setTableValues( values );
  m_CO2SolubilityTable->reInitializeFunction();
  m_CO2SolubilityTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
}

void CO2Solubility::calculateCO2Solubility( array1d< array1d< real64 > > const & coordinates,
                                            real64 const & salinity,
                                            array1d< real64 > const & values )
{
  // Interaction parameters, see Table 2 of Duan and Sun (2003)
  constexpr real64 mu[] =
  { 28.9447706, -0.0354581768, -4770.67077, 1.02782768e-5, 33.8126098, 9.04037140e-3,
    -1.14934031e-3, -0.307405726, -0.0907301486, 9.32713393e-4, 0 };
  constexpr real64 lambda[] = { -0.411370585, 6.07632013e-4, 97.5347708, 0, 0, 0, 0, -0.0237622469, 0.0170656236, 0, 1.41335834e-5 };
  constexpr real64 zeta[] = { 3.36389723e-4, -1.98298980e-5, 0, 0, 0, 0, 0, 2.12220830e-3, -5.24873303e-3, 0, 0 };

  localIndex const numPressures = coordinates[0].size();
  localIndex const numTemperatures = coordinates[1].size();

  for( localIndex i = 0; i < numPressures; ++i )
  {
    real64 const P = coordinates[0][i] / P_Pa_f;

    for( localIndex j = 0; j < numTemperatures; ++j )
    {
      real64 const T = coordinates[1][j];

      // compute reduced volume by solving the CO2 equation of state
      real64 V_r = 0.0;
      CO2SolubilityFunction( T, P, V_r, &detail::ff );

      // compute equation (6) of Duan and Sun (2003)
      real64 const logK = detail::Par( T+T_K_f, P, mu )
                          - detail::logF( T, P, V_r )
                          + 2*detail::Par( T+T_K_f, P, lambda ) * salinity
                          + detail::Par( T+T_K_f, P, zeta ) * salinity * salinity;

      // mole fraction of CO2 in vapor phase, equation (4) of Duan and Sun (2003)
      real64 const y_CO2 = (P - detail::PWater( T ))/P;
      values[j*numPressures+i] = y_CO2 * P / exp( logK );
    }
  }
}

void CO2Solubility::CO2SolubilityFunction( real64 const & T,
                                           real64 const & P,
                                           real64 & V_r,
                                           real64 (*f)( real64 const & x1, real64 const & x2, real64 const & x3 ) )
{

  constexpr real64 eps = 1e-9;
  constexpr real64 dx = 1e-10;
  int count = 0;
  real64 Vr_int = 0.05;

  V_r = 0.75*Rgas*(T_K_f+T)/(P*P_Pa_f)*(1/V_c);

  // iterate until the solution of the CO2 equation of state is found
  for(;; )
  {
    if( V_r < 0.0 )
    {
      V_r = Vr_int;
      Vr_int += 0.05;
    }

    real64 const v0 = (*f)( T, P, V_r );
    real64 const v1 = (*f)( T, P, V_r+dx );
    real64 const dre = -v0/((v1-v0)/dx);

    if( fabs( dre ) < eps )
    {
      break;
    }

    GEOSX_ERROR_IF( count > 50, "CO2Solubility NR convergence fails! " << "dre = " << dre << ", eps = " << eps );

    count++;
    V_r += dre;
  }
}


CO2Solubility::KernelWrapper CO2Solubility::createKernelWrapper()
{
  return KernelWrapper( m_componentNames,
                        m_componentMolarWeight,
                        m_CO2SolubilityTable,
                        m_CO2Index,
                        m_waterIndex,
                        m_phaseGasIndex,
                        m_phaseLiquidIndex );
}

REGISTER_CATALOG_ENTRY( FlashModelBase, CO2Solubility, array1d< string > const &, array1d< string > const &, array1d< string > const &, array1d< real64 > const & )

} // end namespace PVTProps

} // end namespace geosx
