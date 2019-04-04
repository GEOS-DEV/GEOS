
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

BrineViscosityFunction::BrineViscosityFunction( string_array const & inputPara,
                                                string_array const & componentNames,
                                                real64_array const & componentMolarWeight):
  PVTFunctionBase( inputPara[1], componentNames, componentMolarWeight)
{

  MakeCoef(inputPara);

}

void BrineViscosityFunction::MakeCoef(string_array const & inputPara)
{

  constexpr real64 a = 0.0816;
  constexpr real64 b = 0.0122;
  constexpr real64 c = 0.000128;
  constexpr real64 d = 0.000629;
  constexpr real64 k = -0.7;

  constexpr real64 waterVisc = 8.9e-4; //at 25C

  real64 m = stod(inputPara[2]);

  m_coef0 = (1.0 + a * m + b * m * m + c * m * m * m) * waterVisc;

  m_coef1 =  d * (1.0 - exp(k * m)) * waterVisc;

}


void BrineViscosityFunction::Evaluation(const EvalVarArgs& pressure, const EvalVarArgs& temperature, const array1dT<EvalVarArgs>& phaseComposition, EvalVarArgs& value, bool useMass) const
{

  value = m_coef0 + m_coef1 * temperature;

}

REGISTER_CATALOG_ENTRY( PVTFunctionBase,
                        BrineViscosityFunction,
                        string_array const &, string_array const &, real64_array const & )

}

}
