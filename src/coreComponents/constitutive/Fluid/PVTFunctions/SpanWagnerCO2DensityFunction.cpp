
/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file SpanWagnerCO2DensityFunction.cpp
 */


#include "constitutive/Fluid/PVTFunctions/SpanWagnerCO2DensityFunction.hpp"

using namespace std;

namespace geosx
{

using namespace stringutilities;

namespace PVTProps
{

  //calc density based on table 3
  real64 f(const real64 &T, const real64 &P, const real64 &rho)
  {
    constexpr real64 n[] = {0.38856823203161, 2.938547594274, -5.5867188534934, -0.76753199592477, 0.31729005580416, 0.54803315897767, 0.12279411220335, 2.165896154322, 1.5841735109724, -0.23132705405503, 0.058116916431436, -0.55369137205382, 0.48946615909422, -0.024275739843501, 0.062494790501678,-0.12175860225246, -0.37055685270086, -0.016775879700426, -0.11960736637987, -0.045619362508778, 0.035612789270346, -0.0074427727132052, -0.0017395704902432, -0.021810121289527,  0.024332166559236,-0.037440133423463, 0.14338715756878, -0.13491969083286, -0.02315122505348, 0.012363125492901, 0.002105832197294, -0.00033958519026368, 0.0055993651771592, -0.00030335118055646, -213.65488688320, 26641.569149272, -24027.212204557, -283.41603423999, 212.47284400179, -0.66642276540751, 0.72608632349897, 0.055068668612842};

    constexpr real64 d[] = {1, 1, 1, 1, 2, 2, 3, 1, 2, 4, 5, 5, 5, 6, 6, 6, 1, 1, 4, 4, 4, 7, 8, 2, 3, 3, 5, 5, 6, 7, 8, 10, 4, 8, 2, 2, 2, 3, 3};

    constexpr real64 t[] = {0, 0.75, 1, 2, 0.75, 2, 0.75, 1.5, 1.5, 2.5, 0, 1.5, 2, 0, 1, 2, 3, 6, 3, 6, 8, 6, 0, 7, 12, 16, 22, 24, 16, 24, 8, 2, 28, 14, 1, 0, 1, 3, 3};
    
    constexpr real64 c[] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2 ,3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 6};

    constexpr real64 alpha[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25, 25, 25, 15, 20};

    constexpr real64 beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 325, 300, 300, 275, 275, 0.3, 0.3, 0.3};

    constexpr real64 gamma0[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.16, 1.19, 1.19, 1.25, 1.22};

    constexpr real64 epslon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1., 1., 1., 1., 1.};

    constexpr real64 a[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.5, 3.5, 3.};

    constexpr real64 b[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.875, 0.925, 0.875};

    constexpr real64 A[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0.7, 0.7, 0.7};

    constexpr real64 B[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3, 0.3, 1.0};

    constexpr real64 C[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10.0, 10.0, 12.5};

    constexpr real64 D[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 275, 275, 275};

    constexpr real64 T_c = 304.1282;
    constexpr real64 rho_c = 467.6;

    constexpr real64 R = 188.9241;
    
    real64 tau = T_c / T;
    real64 delta = rho / rho_c;

    real64 theta, Delta, Phi, dPhi_delta, dDelta_delta, dDelta_delta_b;

    real64 phi_r_delta = 0.0;

    for(int i=0; i<7; i++)
    {
      phi_r_delta += n[i] * d[i] * pow(delta, d[i]-1.0) * pow(tau, t[i]);
    }

    for(int i=7; i<34; i++)
    {
      phi_r_delta += n[i] * exp(-pow(delta, c[i])) * pow(delta, d[i] - 1.0) * pow(tau, t[i]) * (d[i] - c[i] * pow(delta, c[i]));
    }

    for(int i=34; i<39; i++)
    {
      phi_r_delta += n[i] * pow(delta, d[i]) * pow(tau, t[i]) * exp(-alpha[i] * pow(delta - epslon[i], 2.0) - beta[i] * pow(tau - gamma0[i], 2.0)) * (d[i] / delta - 2.0 * alpha[i] * (delta - epslon[i]));
    }

    for(int i=39; i<42; i++)
    {
      theta = 1.0 - tau + A[i] * pow(pow(delta - 1.0, 2.0), 1.0 / (2.0 * beta[i]));
      Delta = pow(theta, 2.0) + B[i] * pow(pow(delta - 1.0, 2.0), a[i]);

      Phi = exp(-C[i] * pow(delta - 1.0, 2.0) - D[i] * pow(tau - 1.0, 2.0));

      dPhi_delta = -2.0 * C[i] * (delta - 1.0) * Phi;

      dDelta_delta = (delta - 1.0) * (A[i] * theta * 2.0 / beta[i] * pow(pow(delta - 1.0, 2.0), 1.0 / (2.0 * beta[i]) - 1.0) + 2.0 * B[i] * a[i] * pow(pow(delta - 1.0, 2.0), a[i] - 1.0));

      dDelta_delta_b = b[i] * pow(Delta, b[i] - 1.0) * dDelta_delta;

      phi_r_delta += n[i] * (pow(Delta, b[i]) * (Phi + delta *dPhi_delta) + dDelta_delta_b * delta * Phi);
    }

    return rho * (1.0 + delta * phi_r_delta) / (P / (R * T)) - 1.0;
  }

  
SpanWagnerCO2DensityFunction::SpanWagnerCO2DensityFunction( const string_array& inputPara, const string_array& componentNames, const real64_array& componentMolarWeight): PVTFunctionBase( inputPara[1], componentNames, componentMolarWeight)
{

  bool notFound = 1;

  for(localIndex i = 0; i < componentNames.size(); ++i)
  {

    if(streq(componentNames[i], "CO2") || streq(componentNames[i], "co2"))
    {
      m_CO2Index = i;
      notFound = 0;
      break;
    }

  }

  GEOS_ERROR_IF(notFound, "Component CO2 is not found!");

  MakeTable(inputPara);

}

void SpanWagnerCO2DensityFunction::MakeTable(const string_array& inputPara)
{

  real64_vector pressures;
  real64_vector temperatures;

  real64 PStart, PEnd, dP;
  real64 TStart, TEnd, dT;
  real64 P, T;

  PStart = stod(inputPara[2]);
  PEnd = stod(inputPara[3]);
  dP = stod(inputPara[4]);

  TStart = stod(inputPara[5]);
  TEnd = stod(inputPara[6]);
  dT = stod(inputPara[7]);


  P = PStart;

  while(P <= PEnd)
  {

    pressures.push_back(P);
    P += dP;

  }

  T = TStart;

  while(T <= TEnd)
  {

    temperatures.push_back(T);
    T += dT;

  }

  unsigned long nP = pressures.size();
  unsigned long nT = temperatures.size();

  array1dT<real64_vector> densities(nP);
  for(unsigned long i = 0; i < nP; ++i)
  {
    densities[i].resize(nT);
  }

  CalculateCO2Density(pressures, temperatures, densities);

  m_CO2DensityTable = make_shared<XYTable>("SpanWagnerCO2DensityTable", pressures, temperatures, densities);


}


void SpanWagnerCO2DensityFunction::Evaluation(const EvalVarArgs& pressure, const EvalVarArgs& temperature, const array1dT<EvalVarArgs>& phaseComposition, EvalVarArgs& value, bool useMass) const
{

  EvalArgs2D P, T, density;
  P.m_var = pressure.m_var;
  P.m_der[0] = 1.0;

  T.m_var = temperature.m_var;
  T.m_der[1] = 1.0;

  density = m_CO2DensityTable->Value(P, T);

  real64 CO2MW = m_componentMolarWeight[m_CO2Index];

  if(!useMass)
  {

    density /= CO2MW;

  }

  value.m_var = density.m_var;
  value.m_der[0] = density.m_der[0];

}

void SpanWagnerCO2DensityFunction::CalculateCO2Density(const real64_vector& pressure, const real64_vector& temperature, array1dT<real64_vector>& density)
{

  constexpr real64 T_K_f = 273.15;
  
  real64 PPa, TK;
  
  for(unsigned long i = 0; i < pressure.size(); ++i)
  {

    PPa = pressure[i];
      
    for(unsigned long j = 0; j < temperature.size(); ++j)    
    {

      TK = temperature[j] + T_K_f;

      SpanWagnerCO2Density(TK, PPa, density[i][j], &f);

    }

  }
}

void SpanWagnerCO2DensityFunction::SpanWagnerCO2Density(const real64 &T, const real64 &P, real64 &rho, real64 (*f)(const real64 &x1, const real64 &x2, const real64 &x3))
{
  constexpr real64 P_Pa_f = 1e+5;

  constexpr real64 P_c = 73.773 * P_Pa_f;
  constexpr real64 T_c = 304.1282;
  constexpr real64 rho_c = 467.6;

  constexpr real64 R = 188.9241;

  constexpr real64 Rgas = 8.314467;

  constexpr real64 V_c = Rgas*T_c/P_c;

  
  constexpr real64 vpa[] = {-7.0602087, 1.9391218, -1.6463597, -3.2995634};
  constexpr real64 vpt[] = {1.0, 1.5, 2.0, 4.0};

  constexpr real64 lda[] = {1.9245108, -0.62385555, -0.32731127, 0.39245142};
  constexpr real64 ldt[] = {0.340, 0.5, 1.6666666667, 1.833333333};

  constexpr real64 vda[] = {-1.7074879, -0.82274670, -4.6008549, -10.111178, -29.742252};
  constexpr real64 vdt[] = {0.340, 0.5, 1.0, 2.3333333333, 4.6666666667};


  real64 eps = 1e-10;
  int count = 0;
  real64 dx = 1e-10;

  real64 dre, sum;

  real64 psat, denl;

  //initial guess
  if(T < T_c)
  {
    sum = 0;
    for(int i=0; i<4; i++)
    {
      sum += vpa[i] * pow(1.0 - T / T_c, vpt[i]);
    }

    psat = P_c * exp(T_c / T * sum);

    sum = 0;
    for(int i=0; i<4; i++)
    {
      sum += lda[i] * pow(1.0 - T / T_c, ldt[i]);
    }

    denl = rho_c * exp(sum);

    /*
    sum = 0;
    for(int i=0; i<5; i++)
      sum += vda[i] * pow(1.0 - T / T_c, vdt[i]);

          denv = rho_c * exp(sum);
    */
    
    if(P >= psat)
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

  real64 v1, v0;

  for(;;)
  {
    v0 = (*f)(T, P, rho);
    v1 = (*f)(T, P, rho+dx);
    dre = -v0/((v1-v0)/dx);

    if(fabs(dre) < eps) break;

    GEOS_ERROR_IF(count > 50, "SpanWagnerCO2Density NR convergence fails! " << "dre = " << dre << ", eps = " << eps);

    count++;

    rho += dre;
  }
}



REGISTER_CATALOG_ENTRY( PVTFunctionBase,
                        SpanWagnerCO2DensityFunction,
                        string_array const &, string_array const &, real64_array const & )

}

}
