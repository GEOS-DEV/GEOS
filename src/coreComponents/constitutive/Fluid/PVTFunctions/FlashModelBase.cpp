
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
  * @file FlashModelBase.cpp
  */

#include "FlashModelBase.hpp"

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
