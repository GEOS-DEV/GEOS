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
 * @file SpanWagnerCO2Density.cpp
 */

#include "constitutive/fluid/PVTFunctions/SpanWagnerCO2Density.hpp"

#include "functions/FunctionManager.hpp"

namespace geosx
{

using namespace stringutilities;

namespace constitutive
{

namespace PVTProps
{

namespace detail
{

real64 f( real64 const & T,
          real64 const & P,
          real64 const & rho )
{
  // all the coefficients below are defined in Table (31) of Span and Wagner (1996)
  // the variable names used below follow the notation of the Table
  constexpr real64 n[] =
  {0.38856823203161, 2.938547594274, -5.5867188534934, -0.76753199592477, 0.31729005580416, 0.54803315897767, 0.12279411220335, 2.165896154322,
   1.5841735109724, -0.23132705405503, 0.058116916431436, -0.55369137205382, 0.48946615909422, -0.024275739843501, 0.062494790501678, -0.12175860225246,
   -0.37055685270086, -0.016775879700426, -0.11960736637987, -0.045619362508778, 0.035612789270346, -0.0074427727132052, -0.0017395704902432,
   -0.021810121289527, 0.024332166559236, -0.037440133423463, 0.14338715756878, -0.13491969083286, -0.02315122505348, 0.012363125492901, 0.002105832197294,
   -0.00033958519026368, 0.0055993651771592, -0.00030335118055646, -213.65488688320, 26641.569149272, -24027.212204557, -283.41603423999, 212.47284400179,
   -0.66642276540751, 0.72608632349897, 0.055068668612842};

  constexpr real64 d[] = {1, 1, 1, 1, 2, 2, 3, 1, 2, 4, 5, 5, 5, 6, 6, 6, 1, 1, 4, 4, 4, 7, 8, 2, 3, 3, 5, 5, 6, 7, 8, 10, 4, 8, 2, 2, 2, 3, 3};

  constexpr real64 t[] =
  {0, 0.75, 1, 2, 0.75, 2, 0.75, 1.5, 1.5, 2.5, 0, 1.5, 2, 0, 1, 2, 3, 6, 3, 6, 8, 6, 0, 7, 12, 16, 22, 24, 16, 24, 8, 2, 28, 14, 1, 0, 1, 3, 3};

  constexpr real64 c[] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 6};

  constexpr real64 alpha[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25, 25, 25, 15, 20};

  constexpr real64 beta[] =
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 325, 300, 300, 275, 275, 0.3, 0.3, 0.3};

  constexpr real64 gamma0[] =
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.16, 1.19, 1.19, 1.25, 1.22};

  constexpr real64 epslon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1., 1., 1., 1., 1.};

  constexpr real64 a[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.5, 3.5, 3.};

  constexpr real64 b[] =
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.875, 0.925, 0.875};

  constexpr real64 A[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.7, 0.7, 0.7};

  constexpr real64 B[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3, 0.3, 1.0};

  constexpr real64 C[] =
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10.0, 10.0, 12.5};

  constexpr real64 D[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 275, 275, 275};

  constexpr real64 T_c = 304.1282;
  constexpr real64 rho_c = 467.6;

  constexpr real64 R = 188.9241;

  real64 const tau = T_c / T;
  real64 const delta = rho / rho_c;

  real64 theta, Delta, Phi, dPhi_delta, dDelta_delta, dDelta_delta_b;

  // residual part of the energy equation
  // this is defined in Table (32) of Span and Wagner (1996)
  real64 phi_r_delta = 0.0;

  for( int i=0; i<7; i++ )
  {
    phi_r_delta += n[i] * d[i] * pow( delta, d[i]-1.0 ) * pow( tau, t[i] );
  }

  for( int i=7; i<34; i++ )
  {
    phi_r_delta += n[i] * exp( -pow( delta, c[i] )) * pow( delta, d[i] - 1.0 ) * pow( tau, t[i] ) * (d[i] - c[i] * pow( delta, c[i] ));
  }

  for( int i=34; i<39; i++ )
  {
    phi_r_delta += n[i] *
                   pow( delta,
                        d[i] ) *
                   pow( tau,
                        t[i] ) *
                   exp( -alpha[i] *
                        pow( delta - epslon[i], 2.0 ) - beta[i] * pow( tau - gamma0[i], 2.0 )) * (d[i] / delta - 2.0 * alpha[i] * (delta - epslon[i]));
  }

  for( int i=39; i<42; i++ )
  {
    theta = 1.0 - tau + A[i] * pow( pow( delta - 1.0, 2.0 ), 1.0 / (2.0 * beta[i]));
    Delta = pow( theta, 2.0 ) + B[i] * pow( pow( delta - 1.0, 2.0 ), a[i] );

    Phi = exp( -C[i] * pow( delta - 1.0, 2.0 ) - D[i] * pow( tau - 1.0, 2.0 ));

    dPhi_delta = -2.0 * C[i] * (delta - 1.0) * Phi;

    dDelta_delta = (delta - 1.0) *
                   (A[i] * theta * 2.0 / beta[i] *
                    pow( pow( delta - 1.0, 2.0 ), 1.0 / (2.0 * beta[i]) - 1.0 ) + 2.0 * B[i] * a[i] * pow( pow( delta - 1.0, 2.0 ), a[i] - 1.0 ));

    dDelta_delta_b = b[i] * pow( Delta, b[i] - 1.0 ) * dDelta_delta;

    phi_r_delta += n[i] * (pow( Delta, b[i] ) * (Phi + delta *dPhi_delta) + dDelta_delta_b * delta * Phi);
  }

  // Helmholtz energy equation: see equation (2.2) from Span and Wagner (1996)
  return rho * (1.0 + delta * phi_r_delta) / (P / (R * T)) - 1.0;
}

}

SpanWagnerCO2Density::SpanWagnerCO2Density( string_array const & inputPara,
                                            string_array const & componentNames,
                                            array1d< real64 > const & componentMolarWeight ):
  PVTFunctionBase( inputPara[1],
                   componentNames,
                   componentMolarWeight )
{
  char const * expectedCO2ComponentNames[] = { "CO2", "co2" };
  m_CO2Index = PVTFunctionHelpers::findName( componentNames, expectedCO2ComponentNames );
  GEOSX_THROW_IF( m_CO2Index < 0 || m_CO2Index >= componentNames.size(),
                  "Component CO2 is not found!",
                  InputError );

  makeTable( inputPara );
}

void SpanWagnerCO2Density::makeTable( string_array const & inputPara )
{
  PTTableCoordinates tableCoords;
  PVTFunctionHelpers::initializePropertyTable( inputPara, tableCoords );

  real64 tolerance = 1e-10;
  try
  {
    if( inputPara.size() >= 9 )
    {
      tolerance = stod( inputPara[8] );
    }
  }
  catch( const std::invalid_argument & e )
  {
    GEOSX_THROW( "Invalid property argument:" + string( e.what()),
                 InputError );
  }

  array1d< real64 > densities( tableCoords.nPressures() * tableCoords.nTemperatures() );
  calculateCO2Density( tolerance, tableCoords, densities );

  FunctionManager & functionManager = FunctionManager::getInstance();
  m_CO2DensityTable = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", "CO2DensityTable" ) );
  m_CO2DensityTable->setTableCoordinates( tableCoords.getCoords() );
  m_CO2DensityTable->setTableValues( densities );
  m_CO2DensityTable->reInitializeFunction();
  m_CO2DensityTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
}


void SpanWagnerCO2Density::calculateCO2Density( real64 const & tolerance,
                                                PTTableCoordinates const & tableCoords,
                                                array1d< real64 > const & densities )
{

  constexpr real64 TK_f = 273.15;

  localIndex const nPressures = tableCoords.nPressures();
  localIndex const nTemperatures = tableCoords.nTemperatures();

  for( localIndex i = 0; i < nPressures; ++i )
  {
    real64 const PPa = tableCoords.getPressure( i );
    for( localIndex j = 0; j < nTemperatures; ++j )
    {
      real64 const TK = tableCoords.getTemperature( j ) + TK_f;
      spanWagnerCO2DensityFunction( tolerance, TK, PPa, densities[j*nPressures+i], &detail::f );
    }
  }
}

void SpanWagnerCO2Density::spanWagnerCO2DensityFunction( real64 const & tolerance,
                                                         real64 const & T,
                                                         real64 const & P,
                                                         real64 & rho,
                                                         real64 (*f)( real64 const & x1, real64 const & x2, real64 const & x3 ) )
{
  constexpr real64 P_Pa_f = 1e+5;
  constexpr real64 P_c = 73.773 * P_Pa_f;
  constexpr real64 T_c = 304.1282;
  constexpr real64 rho_c = 467.6;
  constexpr real64 R = 188.9241;

  constexpr real64 vpa[] = { -7.0602087, 1.9391218, -1.6463597, -3.2995634 };
  constexpr real64 vpt[] = { 1.0, 1.5, 2.0, 4.0 };

  constexpr real64 lda[] = { 1.9245108, -0.62385555, -0.32731127, 0.39245142 };
  constexpr real64 ldt[] = { 0.340, 0.5, 1.6666666667, 1.833333333 };

  real64 const dx = 1e-10;
  int count = 0;

  // initial guess
  if( T < T_c )
  {
    real64 sum = 0;
    for( int i=0; i < 4; i++ )
    {
      sum += vpa[i] * pow( 1.0 - T / T_c, vpt[i] );
    }
    real64 const psat = P_c * exp( T_c / T * sum );

    sum = 0;
    for( int i=0; i < 4; i++ )
    {
      sum += lda[i] * pow( 1.0 - T / T_c, ldt[i] );
    }
    real64 const denl = rho_c * exp( sum );

    if( P >= psat )
    {
      rho = denl;
    }
    else
    {
      rho = P / (R * T);
    }
  }
  else
  {
    rho = 1.0 / (30.0 + R * T / P);
  }

  // Newton loop
  for(;; )
  {
    real64 const v0 = (*f)( T, P, rho );
    real64 const v1 = (*f)( T, P, rho+dx );
    real64 const dre = -v0/((v1-v0)/dx);
    if( fabs( dre ) < tolerance )
    {
      break;
    }

    GEOSX_ERROR_IF( count > 50, "SpanWagnerCO2Density NR convergence fails! " << "dre = " << dre << ", tolerance = " << tolerance );

    count++;
    rho += dre;
  }
}

SpanWagnerCO2Density::KernelWrapper SpanWagnerCO2Density::createKernelWrapper()
{
  return KernelWrapper( m_componentNames,
                        m_componentMolarWeight,
                        m_CO2DensityTable,
                        m_CO2Index );
}

REGISTER_CATALOG_ENTRY( PVTFunctionBase, SpanWagnerCO2Density, string_array const &, string_array const &, array1d< real64 > const & )

} // namespace PVTProps

} // namespace constitutive

} // namespace geosx
