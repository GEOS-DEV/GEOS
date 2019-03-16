
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
  * @file CO2SolubilityFunction.cpp
  */


#include "constitutive/Fluid/PVTFunctions/CO2SolubilityFunction.hpp"

using namespace std;

namespace geosx
{

  using namespace stringutilities;

  
namespace PVTProps
{

  CO2SolubilityFunction::CO2SolubilityFunction(const string_array& inputPara, const string_array& phaseNames, const string_array& componentNames, const real64_array& componentMolarWeight) : FlashModelBase(componentNames, componentMolarWeight), m_modelName(inputPara[1])
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

    notFound = 1;
    
    for(localIndex i = 0; i < componentNames.size(); ++i)
    {

      if(streq(componentNames[i], "Water") || streq(componentNames[i], "water"))
	{
	  m_waterIndex = i;
	  notFound = 0;
	  break;
	}

    }

    GEOS_ERROR_IF(notFound, "Component Water/Brine is not found!");   


    notFound = 1;
    
    for(localIndex i = 0; i < phaseNames.size(); ++i)
    {

      if(streq(phaseNames[i], "CO2") || streq(phaseNames[i], "co2") || streq(phaseNames[i], "gas") || streq(phaseNames[i], "Gas"))
	{
	  m_phaseGasIndex = i;
	  notFound = 0;
	  break;
	}

    }

    GEOS_ERROR_IF(notFound, "Phase co2/gas is not found!");   

    notFound = 1;
    
    for(localIndex i = 0; i < phaseNames.size(); ++i)
    {

      if(streq(phaseNames[i], "Water") || streq(phaseNames[i], "water") || streq(phaseNames[i], "Liquid") || streq(phaseNames[i], "liquid"))
	{
	  m_phaseLiquidIndex = i;
	  notFound = 0;
	  break;
	}

    }

    GEOS_ERROR_IF(notFound, "Phase water/liquid is not found!");   

    MakeTable(inputPara);
    
  }

  void CO2SolubilityFunction::MakeTable(const string_array& inputPara)
  {

    real64_vector pressures;
    real64_vector temperatures;    

    real64 PStart, PEnd, dP;
    real64 TStart, TEnd, dT;    
    real64 P, T, m;
    
    PStart = stod(inputPara[2]);
    PEnd = stod(inputPara[3]);    
    dP = stod(inputPara[4]);

    TStart = stod(inputPara[5]);
    TEnd = stod(inputPara[6]);    
    dT = stod(inputPara[7]);    

    m = stod(inputPara[8]);    

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

    array1dT<real64_vector> solubilities(nP);
    for(unsigned long i = 0; i < nP; ++i)
      {
	solubilities[i].resize(nT);
      }

    CalculateCO2Solubility(pressures, temperatures, m, solubilities);
  
    m_CO2SolubilityTable = make_shared<XYTable>("CO2SolubilityTable", pressures, temperatures, solubilities);

    
  }

  void CO2SolubilityFunction::Partition(const EvalVarArgs& pressure, const EvalVarArgs& temperature, const array1dT<EvalVarArgs>& compFraction, array1dT<EvalVarArgs>& phaseFraction, array1dT<array1dT<EvalVarArgs> >& phaseCompFraction) const      
  {

    EvalArgs2D P, T, solubility;
    P.m_var = pressure.m_var;
    P.m_der[0] = 1.0;

    T.m_var = temperature.m_var;
    T.m_der[1] = 1.0;    

    //solubiltiy mol/kg(water)  X = Csat/W
    solubility = m_CO2SolubilityTable->Value(P, T);

    static real64 waterMW = m_componentMolarWeight[m_waterIndex];
    static real64 CO2MW = m_componentMolarWeight[m_CO2Index];

    solubility *= waterMW;
    
    EvalVarArgs X, Y;

    X.m_var = solubility.m_var;
    X.m_der[0] = solubility.m_der[0];

    //Y = C/W = z/(1-z)

    Y = compFraction[m_CO2Index] / (1.0 - compFraction[m_CO2Index]);

    if(Y < X)
      {
	//liquid phase only
    
	phaseFraction[m_phaseLiquidIndex] = 1.0;
	phaseFraction[m_phaseGasIndex] = 0.0;    	

	phaseCompFraction[m_phaseLiquidIndex] = compFraction;

      }
    else
      {
	// two-phase
	// liquid phase fraction = (Csat + W) / (C + W) = (Csat/W + 1) / (C/W + 1)

	phaseFraction[m_phaseLiquidIndex] = (X + 1.0)/ (Y + 1.0);
	phaseFraction[m_phaseGasIndex] = 1.0 - phaseFraction[m_phaseLiquidIndex];      

	//liquid phase composition  CO2 = Csat / (Csat + W) = (Csat/W) / (Csat/W + 1)

	phaseCompFraction[m_phaseLiquidIndex][m_CO2Index] = X / (X + 1.0);	
	phaseCompFraction[m_phaseLiquidIndex][m_waterIndex] = 1.0 - phaseCompFraction[m_phaseLiquidIndex][m_CO2Index];

	//gas phase composition  CO2 = 1.0
	
	phaseCompFraction[m_phaseGasIndex][m_CO2Index] = 1.0;
	phaseCompFraction[m_phaseGasIndex][m_waterIndex] = 0.0;

      }
  }

}

}
