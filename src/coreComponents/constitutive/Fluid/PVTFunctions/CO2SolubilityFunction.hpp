
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
  * @file CO2SolubilityFunction.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_CO2SOLUBILITYFUNCTION_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_CO2SOLUBILITYFUNCTION_HPP

#include "constitutive/Fluid/PVTFunctions/PVTFunction.hpp"

namespace geosx
{

namespace PVTProps
{
  
  class CO2SolubilityFunction : public FlashModelBase
  {
  public:

    CO2SolubilityFunction(const string_array& inputPara, const string_array& phaseNames, const string_array& componentNames, const real64_array& componentMolarWeight);
    ~CO2SolubilityFunction() {}

    virtual const string& FlashModelName() const
    {
      return m_modelName;
    }
    
    virtual void Partition(const EvalVarArgs& pressure, const EvalVarArgs& temperature, const array1dT<EvalVarArgs>& compFraction, array1dT<EvalVarArgs>& phaseFraction, array1dT<array1dT<EvalVarArgs> >& phaseCompFraction) const;    

  private:

    void MakeTable(const string_array& inputPara);    
    
    TableFunctionPtr m_CO2SolubilityTable; 
    string m_modelName;
    localIndex m_CO2Index;
    localIndex m_waterIndex;    
    localIndex m_phaseGasIndex;
    localIndex m_phaseLiquidIndex;    
  };  

}

}
  
#endif
