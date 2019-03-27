
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

SpanWagnerCO2DensityFunction::SpanWagnerCO2DensityFunction( const string_array& inputPara,
                                                            const string_array& componentNames,
                                                            const real64_array& componentMolarWeight):
  PVTFunctionBase( inputPara[1], componentNames, componentMolarWeight)
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

REGISTER_CATALOG_ENTRY( PVTFunctionBase,
                        SpanWagnerCO2DensityFunction,
                        string_array const &, string_array const &, real64_array const & )

}

}
