
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
 * @file FenghourCO2ViscosityFunction.cpp
 */


#include "constitutive/Fluid/PVTFunctions/FenghourCO2ViscosityFunction.hpp"
#include "constitutive/Fluid/PVTFunctions/SpanWagnerCO2DensityFunction.hpp"

using namespace std;

namespace geosx
{

namespace PVTProps
{

FenghourCO2ViscosityFunction::FenghourCO2ViscosityFunction( const string_array& inputPara,
                                                            const string_array& componentNames,
                                                            const real64_array& componentMolarWeight):
  PVTFunctionBase( inputPara[1], componentNames, componentMolarWeight)
{

  MakeTable(inputPara);

}

void FenghourCO2ViscosityFunction::MakeTable(const string_array& inputPara)
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

  array1dT<real64_vector> viscosities(nP);
  array1dT<real64_vector> densities(nP);  
  for(unsigned long i = 0; i < nP; ++i)
  {
    viscosities[i].resize(nT);
    densities[i].resize(nT);    
  }

  SpanWagnerCO2DensityFunction::CalculateCO2Density(pressures, temperatures, densities);
  
  CalculateCO2Viscosity(pressures, temperatures, densities, viscosities);

  m_CO2ViscosityTable = make_shared<XYTable>("FenghourCO2ViscosityTable", pressures, temperatures, viscosities);


}


void FenghourCO2ViscosityFunction::Evaluation(const EvalVarArgs& pressure, const EvalVarArgs& temperature, const array1dT<EvalVarArgs>& phaseComposition, EvalVarArgs& value, bool useMass) const
{

  EvalArgs2D P, T, viscosity;
  P.m_var = pressure.m_var;
  P.m_der[0] = 1.0;

  T.m_var = temperature.m_var;
  T.m_der[1] = 1.0;

  viscosity = m_CO2ViscosityTable->Value(P, T);

  value.m_var = viscosity.m_var;
  value.m_der[0] = viscosity.m_der[0];

}

void FenghourCO2ViscosityFunction::FenghourCO2Viscosity(const real64 &Tcent, const real64 &den, real64 &vis)
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
  const real64 lnGfun = aa[0] + x * (aa[1] + x * (aa[2] + x *(aa[3] + x * aa[4])));
  const real64 GfunInv = exp(-lnGfun);
  const real64 vlimit = 1.00697 * sqrt(Tkelvin) * GfunInv;

  const real64 d2 = den * den;
  const real64 vxcess = den * (d11 + den * (d21 + d2*d2*(d64 / (Tred*Tred*Tred) + d2*(d81 + d82/Tred))));
    
  static const real64 vcrit = 0.0;

  vis = 1e-6 * (vlimit + vxcess + vcrit);

}

  void FenghourCO2ViscosityFunction::CalculateCO2Viscosity(const real64_vector& pressure, const real64_vector& temperature, const array1dT<real64_vector>& density, array1dT<real64_vector>& viscosity)
{

  real64 PPa, TK, rho;
  
  for(unsigned long i = 0; i < pressure.size(); ++i)
    {

      PPa = pressure[i];
      
      for(unsigned long j = 0; j < temperature.size(); ++j)    
	{

	  FenghourCO2Viscosity(temperature[j], density[i][j], viscosity[i][j]);
	    
	}

    }

}
  
  
REGISTER_CATALOG_ENTRY( PVTFunctionBase,
                        FenghourCO2ViscosityFunction,
                        string_array const &, string_array const &, real64_array const & )

}

}
