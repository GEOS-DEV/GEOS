
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

using namespace std;

namespace geosx
{

namespace PVTProps
{

  FenghourCO2ViscosityFunction::FenghourCO2ViscosityFunction(const string_array& inputPara, const string_array& componentNames, const real64_array& componentMolarWeight) : PVTFunctionBase(componentNames, componentMolarWeight), m_functionName(inputPara[1])
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
    for(unsigned long i = 0; i < nP; ++i)
      {
	viscosities[i].resize(nT);
      }
    
    CalculateCO2Viscosity(pressures, temperatures, viscosities);

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

}

}
