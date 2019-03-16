
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
  * @file PVTFunction.cpp
  */

#include "constitutive/Fluid/PVTFunctions/PVTFunction.hpp"

using namespace std;

namespace geosx
{

namespace PVTProps
{

  static real64 T_K_f = 273.15;
  static real64 P_Pa_f = 1e+5;

  static real64 P_c = 73.773 * P_Pa_f;
  static real64 T_c = 304.1282;
  static real64 rho_c = 467.6;

  static real64 R = 188.9241;

  static real64 Rgas = 8.314467;

  static real64 V_c = Rgas*T_c/P_c;
  
  //Table 31

  static real64 n[] = {0.38856823203161, 2.938547594274, -5.5867188534934, -0.76753199592477, 0.31729005580416, 0.54803315897767, 0.12279411220335, 2.165896154322, 1.5841735109724, -0.23132705405503, 0.058116916431436, -0.55369137205382, 0.48946615909422, -0.024275739843501, 0.062494790501678,-0.12175860225246, -0.37055685270086, -0.016775879700426, -0.11960736637987, -0.045619362508778, 0.035612789270346, -0.0074427727132052, -0.0017395704902432, -0.021810121289527,  0.024332166559236,-0.037440133423463, 0.14338715756878, -0.13491969083286, -0.02315122505348, 0.012363125492901, 0.002105832197294, -0.00033958519026368, 0.0055993651771592, -0.00030335118055646, -213.65488688320, 26641.569149272, -24027.212204557, -283.41603423999, 212.47284400179, -0.66642276540751, 0.72608632349897, 0.055068668612842};

  static real64 d[] = {1, 1, 1, 1, 2, 2, 3, 1, 2, 4, 5, 5, 5, 6, 6, 6, 1, 1, 4, 4, 4, 7, 8, 2, 3, 3, 5, 5, 6, 7, 8, 10, 4, 8, 2, 2, 2, 3, 3};

  static real64 t[] = {0, 0.75, 1, 2, 0.75, 2, 0.75, 1.5, 1.5, 2.5, 0, 1.5, 2, 0, 1, 2, 3, 6, 3, 6, 8, 6, 0, 7, 12, 16, 22, 24, 16, 24, 8, 2, 28, 14, 1, 0, 1, 3, 3};

  static real64 c[] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2 ,3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 6};

  static real64 alpha[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25, 25, 25, 15, 20};

  static real64 beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 325, 300, 300, 275, 275, 0.3, 0.3, 0.3};

  static real64 gamma0[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.16, 1.19, 1.19, 1.25, 1.22};

  static real64 epslon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1., 1., 1., 1., 1.};

  static real64 a[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.5, 3.5, 3.};

  static real64 b[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.875, 0.925, 0.875};

  static real64 A[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0.7, 0.7, 0.7};

  static real64 B[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3, 0.3, 1.0};

  static real64 C[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10.0, 10.0, 12.5};

  static real64 D[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 275, 275, 275};

  static real64 vpa[] = {-7.0602087, 1.9391218, -1.6463597, -3.2995634};
  static real64 vpt[] = {1.0, 1.5, 2.0, 4.0};

  static real64 lda[] = {1.9245108, -0.62385555, -0.32731127, 0.39245142};
  static real64 ldt[] = {0.340, 0.5, 1.6666666667, 1.833333333};

  static real64 vda[] = {-1.7074879, -0.82274670, -4.6008549, -10.111178, -29.742252};
  static real64 vdt[] = {0.340, 0.5, 1.0, 2.3333333333, 4.6666666667};

  static void SpanWagnerCO2Density(const real64 &T, const real64 &P, real64 &rho, real64 (*f)(const real64 &x1, const real64 &x2, const real64 &x3))
  {

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
	  sum += vpa[i] * pow(1.0 - T / T_c, vpt[i]);

	psat = P_c * exp(T_c / T * sum);

	sum = 0;
	for(int i=0; i<4; i++)
	  sum += lda[i] * pow(1.0 - T / T_c, ldt[i]);

	denl = rho_c * exp(sum);

	/*
	sum = 0;
	for(int i=0; i<5; i++)
	  sum += vda[i] * pow(1.0 - T / T_c, vdt[i]);

        denv = rho_c * exp(sum);
	*/
	
	if(P >= psat)
	  rho = denl;
	else  
	  rho = P / (R * T);

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

	if(fabs(dre) < eps) 
	  break;

	if(count > 50)
	  {

	    GEOS_ERROR( "SpanWagnerCO2Density NR convergence fails!");

	  }

	count++;

	rho += dre;

      }
    
  }

  //calc density based on table 3
  real64 f(const real64 &T, const real64 &P, const real64 &rho)
  {

    real64 tau = T_c / T;
    real64 delta = rho / rho_c;

    real64 theta, Delta, Phi, dPhi_delta, dDelta_delta, dDelta_delta_b;

    real64 phi_r_delta = 0.0;

    for(int i=0; i<7; i++) 
      phi_r_delta += n[i] * d[i] * pow(delta, d[i]-1.0) * pow(tau, t[i]);

    for(int i=7; i<34; i++) 
      phi_r_delta += n[i] * exp(-pow(delta, c[i])) * pow(delta, d[i] - 1.0) * pow(tau, t[i]) * (d[i] - c[i] * pow(delta, c[i]));

    for(int i=34; i<39; i++) 
      phi_r_delta += n[i] * pow(delta, d[i]) * pow(tau, t[i]) * exp(-alpha[i] * pow(delta - epslon[i], 2.0) - beta[i] * pow(tau - gamma0[i], 2.0)) * (d[i] / delta - 2.0 * alpha[i] * (delta - epslon[i]));

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

  //calc fugacity coef.  based on table 3

  real64 fugacity_c(const real64 &T, const real64 &P, const real64 &rho)
  {
    
    real64 tau = T_c / T;
    real64 delta = rho / rho_c;

    real64 theta, Delta, Phi, dPhi_delta, dDelta_delta, dDelta_delta_b;

    real64 phi_r = 0.0;
    real64 phi_r_delta = 0.0;

    for(int i=0; i<7; i++)
      {
	phi_r += n[i] * pow(delta, d[i]) * pow(tau, t[i]);
	phi_r_delta += n[i] * d[i] * pow(delta, d[i]-1.0) * pow(tau, t[i]);
      }

    for(int i=7; i<34; i++)
      {
	phi_r += n[i] * pow(delta, d[i]) * pow(tau, t[i]) * exp(-pow(delta, c[i]));
	phi_r_delta += n[i] * exp(-pow(delta, c[i])) * pow(delta, d[i] - 1.0) * pow(tau, t[i]) * (d[i] - c[i] * pow(delta, c[i]));
      }

    for(int i=34; i<39; i++)
      {
	phi_r += n[i] * pow(delta, d[i]) * pow(tau, t[i]) * exp(-alpha[i] * pow(delta - epslon[i], 2.0) - beta[i] * pow(tau - gamma0[i], 2.0));
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

	phi_r += n[i] * pow(Delta, b[i]) * delta * Phi;
	phi_r_delta += n[i] * (pow(Delta, b[i]) * (Phi + delta *dPhi_delta) + dDelta_delta_b * delta * Phi);

      }

    return exp(phi_r + delta * phi_r_delta - log(1.0 + delta * phi_r_delta));

  }

  static void FenghourCO2Viscosity(const real64 &Tcent,
				   const real64 &den, real64 &vis)
  {
    static const real64 espar = 251.196;
    static const real64 esparInv = 1.0 / espar;
    static const real64 aa[5] = { 0.235156,  -0.491266, 5.211155e-2, 5.347906e-2, -1.537102e-2 };
    static const real64 d11 =  0.4071119e-2;
    static const real64 d21 =  0.7198037e-4;
    static const real64 d64 =  0.2411697e-16;
    static const real64 d81 =  0.2971072e-22;
    static const real64 d82 = -0.1627888e-22;

    // temperature in Kelvin
    const real64 Tkelvin = Tcent + 273.15;
    // evaluate vlimit from eqns 3-5
    const real64 Tred   = Tkelvin * esparInv;
    const real64 x = log(Tred);
    const real64 lnGfun =
      aa[0] + x * (aa[1] + x * (aa[2] + x *(aa[3] + x * aa[4])));
    const real64 GfunInv = exp(-lnGfun);
    const real64 vlimit = 1.00697 * sqrt(Tkelvin) * GfunInv;

    const real64 d2 = den * den;
    const real64 vxcess =
      den * (d11 + den * (d21 + d2*d2*(d64 / (Tred*Tred*Tred) +
				   d2*(d81 + d82/Tred))));
    
    static const real64 vcrit = 0.0;

    vis = 1e-6 * (vlimit + vxcess + vcrit);

  }

  void CalculateCO2Density(const real64_vector& pressure, const real64_vector& temperature, array1dT<real64_vector>& density)
  {

    real64 PPa, TK, rho;
  
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


  void CalculateCO2Viscosity(const real64_vector& pressure, const real64_vector& temperature, array1dT<real64_vector>& viscosity)
  {

    real64 PPa, TK, rho;
  
    for(unsigned long i = 0; i < pressure.size(); ++i)
      {

	PPa = pressure[i];
      
	for(unsigned long j = 0; j < temperature.size(); ++j)    
	  {

	    TK = temperature[j] + T_K_f;
	    SpanWagnerCO2Density(TK, PPa, rho, &f);
	    FenghourCO2Viscosity(temperature[j], rho, viscosity[i][j]);
	    
	  }

      }

  }


  void CalculateBrineDensity(const real64_vector& pressure, const real64_vector& temperature, const real64& salinity, array1dT<real64_vector>& density)
  {
  
    static const real64 c1 = -9.9595;
    static const real64 c2 = 7.0845;  
    static const real64 c3 = 3.9093;

    static const real64 a1 = -0.004539;
    static const real64 a2 = -0.0001638;
    static const real64 a3 = 0.00002551;

    static const real64 AA = -3.033405;
    static const real64 BB = 10.128163;
    static const real64 CC = -8.750567;
    static const real64 DD = 2.663107;

    real64 P, x;
  
    for(unsigned long i = 0; i < pressure.size(); ++i)
      {

	P = pressure[i] / 1e5;

	for(unsigned long j = 0; j < temperature.size(); ++j)    
	  {
      
	    x = c1 * exp(a1 * salinity) + c2 * exp(a2 * temperature[j]) + c3 * exp(a3 * P);

	    density[i][j] = (AA + BB * x + CC * x * x + DD * x * x * x) * 1000.0;

	  }
      }
  
  }


  static const double acoef[] = {8.99288497e-2, -4.94783127e-1, 4.77922245e-2, 1.03808883e-2, -2.82516861e-2, 9.49887563e-2, 5.20600880e-4, -2.93540971e-4, -1.77265112e-3, -2.51101973e-5, 8.93353441e-5, 7.88998563e-5, -1.66727022e-2, 1.398, 2.96e-2};

  static const double ccoef[] = {-38.640844, 5.8948420, 59.876516, 26.654627, 10.637097};

  static void nr(const double &T, const double &P, double &V_r, double (*f)(const double &x1, const double &x2, const double &x3))
  {

    double eps = 1e-9;
    int count = 0;
    double dx = 1e-10;

    double dre;
    double Vr_int = 0.05;
  
    V_r = 0.75*Rgas*(T_K_f+T)/(P*P_Pa_f)*(1/V_c);

    double v1, v0;

    for(;;)
      {

	if(V_r < 0.0)
	  {
	    V_r = Vr_int;
	    Vr_int += 0.05;
	  }
    
	v0 = (*f)(T, P, V_r);
	v1 = (*f)(T, P, V_r+dx);
	dre = -v0/((v1-v0)/dx);

	if(fabs(dre) < eps) 
	  break;

	if(count > 50)
	  {
	    GEOS_ERROR( "CO2Solubiltiy NR convergence fails!");
	  }

	count++;

	V_r += dre;

      }
  }

  double ff(const double &T, const double &P, const double &V_r)
  {

    double P_r = P*P_Pa_f/P_c;
    double T_r = (T_K_f+T)/T_c;

    double f_Z = 1.0 + (acoef[0] + acoef[1]/T_r/T_r + acoef[2]/T_r/T_r/T_r)/V_r + (acoef[3] + acoef[4]/T_r/T_r + acoef[5]/T_r/T_r/T_r)/V_r/V_r + (acoef[6] + acoef[7]/T_r/T_r + acoef[8]/T_r/T_r/T_r)/V_r/V_r/V_r/V_r + (acoef[9] + acoef[10]/T_r/T_r + acoef[11]/T_r/T_r/T_r)/V_r/V_r/V_r/V_r/V_r + acoef[12]/T_r/T_r/T_r/V_r/V_r * (acoef[13] + acoef[14]/V_r/V_r) * exp(-acoef[14]/V_r/V_r) - P_r * V_r / T_r;

    return f_Z;  

  }

  double PWater(const double &T)
  {

    double P_c_w = 220.85;       // H2O critical pressure (bars)
    double T_c_w = 647.29;     // H2O critical temperature (K)
    double tt = ((T+T_K_f)-T_c_w)/T_c_w;
    double x = (P_c_w*(T+T_K_f)/T_c_w) * (1 + ccoef[0]*pow(-tt, 1.9) + ccoef[1]*tt + ccoef[2]*tt*tt + ccoef[3]*tt*tt*tt + ccoef[4]*tt*tt*tt*tt);

    return x;

  }

  double logF(const double &T, const double &P, const double &V_r)
  {

    double P_r = P*P_Pa_f/P_c;
    double T_r = (T_K_f+T)/T_c;

    double Z=P_r * V_r/T_r;

    double log_f = Z - 1 - log(Z) + (acoef[0] + acoef[1]/T_r/T_r + acoef[2]/T_r/T_r/T_r)/V_r + (acoef[3] + acoef[4]/T_r/T_r + acoef[5]/T_r/T_r/T_r)/2.0/V_r/V_r + (acoef[6] + acoef[7]/T_r/T_r + acoef[8]/T_r/T_r/T_r)/4.0/V_r/V_r/V_r/V_r + (acoef[9] + acoef[10]/T_r/T_r + acoef[11]/T_r/T_r/T_r)/5.0/V_r/V_r/V_r/V_r/V_r + acoef[12]/2.0/T_r/T_r/T_r/acoef[14] * (acoef[13] + 1.0 - (acoef[13] + 1.0 + acoef[14]/V_r/V_r) * exp(-acoef[14]/V_r/V_r));

    return log_f;

  }

  double Par(const double &T, const double &P, const double *cc)
  {

    double x = cc[0] + cc[1]*T +cc[2]/T + cc[3]*T*T + cc[4]/(630.0-T) + cc[5]*P + cc[6]*P*log(T) + cc[7]*P/T + cc[8]*P/(630.0-T) + cc[9]*P*P/(630.0-T)/(630.0-T) + cc[10]*T*log(P);

    return x;

  }

  void CalculateCO2Solubility(const real64_vector& pressure, const real64_vector& temperature, const real64& salinity, array1dT<real64_vector>& solubiltiy)
  {

    double T, P, V_r, m, logK, y_CO2, m_CO2;

    static double mu[] = {28.9447706, -0.0354581768, -4770.67077, 1.02782768e-5, 33.8126098, 9.04037140e-3, -1.14934031e-3, -0.307405726, -0.0907301486, 9.32713393e-4, 0};

    static double lambda[] = {-0.411370585, 6.07632013e-4, 97.5347708, 0, 0, 0, 0, -0.0237622469, 0.0170656236, 0, 1.41335834e-5};

    static double zeta[] = {3.36389723e-4, -1.98298980e-5, 0, 0, 0, 0, 0, 2.12220830e-3, -5.24873303e-3, 0, 0};

    m = salinity;

    for(unsigned long i = 0; i < pressure.size(); ++i)
      {

	P = pressure[i] / P_Pa_f;
      
	for(unsigned long j = 0; j < temperature.size(); ++j)    
	  {

	    T = temperature[j];
	  
	    nr(T, P, V_r, &ff);

	    logK = Par(T+T_K_f,P,mu) - logF(T, P, V_r) + 2*Par(T+T_K_f,P,lambda)*m + Par(T+T_K_f,P,zeta)*m*m;

	    y_CO2 = (P - PWater(T))/P;

	    solubiltiy[i][j] = y_CO2 * P / exp(logK);

	  }
  
      }

  }

}

}
