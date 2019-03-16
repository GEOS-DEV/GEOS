
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
  * @file BrineViscosityFunction.cpp
  */


#include "constitutive/Fluid/PVTFunctions/BrineViscosityFunction.hpp"

using namespace std;

namespace geosx
{

  using namespace stringutilities;  
  
namespace PVTProps
{

  BrineViscosityFunction::BrineViscosityFunction(const string_array& inputPara, const string_array& componentNames, const real64_array& componentMolarWeight) : PVTFunctionBase(componentNames, componentMolarWeight), m_functionName(inputPara[1])
  {

    MakeCoef(inputPara);
    
  }

  void BrineViscosityFunction::MakeCoef(const string_array& inputPara)
  {

    static const real64 a = 0.0816;
    static const real64 b = 0.0122;
    static const real64 c = 0.000128;
    static const real64 d = 0.000629;    
    static const real64 k = -0.7;    

    static const real64 waterVisc = 8.9e-4; //at 25C
    
    real64 m = stod(inputPara[2]);    

    m_coef0 = (1.0 + a * m + b * m * m + c * m * m * m) * waterVisc;

    m_coef1 =  d * (1.0 - exp(k * m)) * waterVisc;
    
  }
    
    
  void BrineViscosityFunction::Evaluation(const EvalVarArgs& pressure, const EvalVarArgs& temperature, const array1dT<EvalVarArgs>& phaseComposition, EvalVarArgs& value, bool useMass) const
  {

    value = m_coef0 + m_coef1 * temperature;
    
  }

}

}
