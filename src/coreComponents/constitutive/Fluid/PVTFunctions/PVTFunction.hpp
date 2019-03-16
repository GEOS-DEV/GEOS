
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
  * @file PVTFunction.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PVTFUNCTION_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PVTFUNCTION_HPP

#include "constitutive/Fluid/PVTFunctions/UtilityFunctions.hpp"
#include "codingUtilities/StringUtilities.hpp"

namespace geosx
{

namespace PVTProps
{

  enum PVTFUNCTYPE {UNKNOWN, DENSITY, VISCOSITY};
  
  class PVTFunctionBase
  {
  public:

    PVTFunctionBase(const string_array& componentNames, const real64_array& componentMolarWeight) : m_componentNames(componentNames), m_componentMolarWeight(componentMolarWeight) {}

    virtual ~PVTFunctionBase(){}
    
    virtual const string& FunctionName() const = 0;
    virtual PVTFUNCTYPE FunctionType() const = 0;    

    //phase density/viscosity
    //input: P, T, phaseCompFraction
    //output: phase density/viscoty
    
    virtual void Evaluation(const EvalVarArgs& pressure, const EvalVarArgs& temperature, const array1dT<EvalVarArgs>& phaseComposition, EvalVarArgs& value, bool useMass = 0) const = 0;

  protected:

    string_array m_componentNames;
    real64_array m_componentMolarWeight;
    
  };  

  typedef std::shared_ptr<PVTFunctionBase> PVTFunction;


  class FlashModelBase
  {
  public:

    FlashModelBase(const string_array& componentNames, const real64_array& componentMolarWeight) : m_componentNames(componentNames), m_componentMolarWeight(componentMolarWeight) {}

    virtual ~FlashModelBase(){}
    
    virtual const string& FlashModelName() const = 0;

    //partition
    //input: P, T, totalCompFraction
    //output: phaseFraction, phaseCompFraction
    
    virtual void Partition(const EvalVarArgs& pressure, const EvalVarArgs& temperature, const array1dT<EvalVarArgs>& compFraction, array1dT<EvalVarArgs>& phaseFraction, array1dT<array1dT<EvalVarArgs> >& phaseCompFraction) const = 0;    

  protected:

    string_array m_componentNames;
    real64_array m_componentMolarWeight;
    
  };  

  typedef std::shared_ptr<FlashModelBase> FlashModel;


  void CalculateCO2Density(const real64_vector& pressure, const real64_vector& temperature, array1dT<real64_vector>& density);

  void CalculateCO2Viscosity(const real64_vector& pressure, const real64_vector& temperature, array1dT<real64_vector>& viscosity);  

  void CalculateBrineDensity(const real64_vector& pressure, const real64_vector& temperature, const real64& salinity, array1dT<real64_vector>& density);

  void CalculateCO2Solubility(const real64_vector& pressure, const real64_vector& temperature, const real64& salinity, array1dT<real64_vector>& solubiltiy);  
  
  
}

}
  
#endif
