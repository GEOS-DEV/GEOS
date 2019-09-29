
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
 * @file PVTFunctionBase.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PVTFUNCTIONBASE_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PVTFUNCTIONBASE_HPP

#include "constitutive/Fluid/PVTFunctions/UtilityFunctions.hpp"
#include "codingUtilities/StringUtilities.hpp"

namespace geosx
{

namespace PVTProps
{

enum class PVTFuncType {UNKNOWN, DENSITY, VISCOSITY};

class PVTFunction
{
public:

  PVTFunction( string const & name, string_array const & componentNames, real64_array const & componentMolarWeight ):
    m_functionName( name ),
    m_componentNames(componentNames),
    m_componentMolarWeight(componentMolarWeight)
  {}

  virtual ~PVTFunction(){}


  using CatalogInterface = cxx_utilities::CatalogInterface< PVTFunction, string_array const &,
							    string_array const &,
							    real64_array const & >;
  static typename CatalogInterface::CatalogType& GetCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }
  virtual string GetCatalogName() = 0;


  string const & FunctionName() const
  {
    return m_functionName;
  }

  virtual PVTFuncType FunctionType() const = 0;

  //phase density/viscosity
  //input: P, T, phaseCompFraction
  //output: phase density/viscoty

  virtual void Evaluation(EvalVarArgs const & pressure, EvalVarArgs const & temperature, arraySlice1d<EvalVarArgs const> const & phaseComposition, EvalVarArgs & value, bool useMass = 0) const = 0;

protected:

  string m_functionName;
  string_array m_componentNames;
  real64_array m_componentMolarWeight;


};

}

}

#endif
