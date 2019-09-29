
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
 * @file SpanWagnerCO2DensityFunction.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_SPANWAGNERCO2DENSITYFUNCTION_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_SPANWAGNERCO2DENSITYFUNCTION_HPP

#include "PVTFunctionBase.hpp"

namespace geosx
{

namespace PVTProps
{

class SpanWagnerCO2DensityFunction : public PVTFunction
{
public:


  SpanWagnerCO2DensityFunction( string_array const & inputPara,
                                string_array const & componentNames,
                                real64_array const & componentMolarWeight);

  ~SpanWagnerCO2DensityFunction() override {}


  static constexpr auto m_catalogName = "SpanWagnerCO2Density";
  static string CatalogName()                    { return m_catalogName; }
  virtual string GetCatalogName() override final { return CatalogName(); }


  virtual PVTFuncType FunctionType() const override
  {
    return PVTFuncType::DENSITY;
  }

  virtual void Evaluation(EvalVarArgs const & pressure, EvalVarArgs const & temperature, arraySlice1d<EvalVarArgs const> const & phaseComposition, EvalVarArgs & value, bool useMass = 0) const override;
  
  static void CalculateCO2Density(real64_const_array const & pressure, real64_const_array const & temperature, real64_array2d const & density);
  
private:

  void MakeTable(string_array const & inputPara);

  static void SpanWagnerCO2Density(real64 const & T, real64 const & P, real64 & rho, real64 (*f)(real64 const & x1, real64 const & x2, real64 const & x3));
  

  TableFunctionPtr m_CO2DensityTable;
  localIndex m_CO2Index;

};

}

}

#endif
