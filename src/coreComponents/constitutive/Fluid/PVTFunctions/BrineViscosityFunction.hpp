
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
  * @file BrineViscosityFunction.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_BRINEVISCOSITYFUNCTION_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_BRINEVISCOSITYFUNCTION_HPP

#include "constitutive/Fluid/PVTFunctions/PVTFunction.hpp"

namespace geosx
{
  
namespace PVTProps
{

  class BrineViscosityFunction : public PVTFunctionBase
  {
  public:


    BrineViscosityFunction(const string_array& inputPara, const string_array& componentNames, const real64_array& componentMolarWeight);
    ~BrineViscosityFunction() {}

    virtual const string& FunctionName() const
    {
      return m_functionName;
    }
    
    virtual PVTFUNCTYPE FunctionType() const
    {
      return PVTFUNCTYPE::VISCOSITY;

    }  

    virtual void Evaluation(const EvalVarArgs& pressure, const EvalVarArgs& temperature, const array1dT<EvalVarArgs>& phaseComposition, EvalVarArgs& value, bool useMass = 0) const;

  private:

    void MakeCoef(const string_array& inputPara);
    
    string m_functionName;
    real64 m_coef0;
    real64 m_coef1;    
    
  };  

}

}  
#endif
