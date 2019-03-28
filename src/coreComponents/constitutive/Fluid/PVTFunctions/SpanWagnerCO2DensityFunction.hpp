
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

class SpanWagnerCO2DensityFunction : public PVTFunctionBase
{
public:


  SpanWagnerCO2DensityFunction( const string_array& inputPara,
                                const string_array& componentNames,
                                const real64_array& componentMolarWeight);

  ~SpanWagnerCO2DensityFunction() override {}


  static constexpr auto m_catalogName = "SpanWagnerCO2Density";
  static string CatalogName()                    { return m_catalogName; }
  virtual string GetCatalogName() override final { return CatalogName(); }


  virtual PVTFUNCTYPE FunctionType() const override
      {
    return PVTFUNCTYPE::DENSITY;

      }

  virtual void Evaluation( const EvalVarArgs& pressure,
                           const EvalVarArgs& temperature,
                           const array1dT<EvalVarArgs>& phaseComposition,
                           EvalVarArgs& value,
                           bool useMass = 0) const override;


  static void CalculateCO2Density(const real64_vector& pressure, const real64_vector& temperature, array1dT<real64_vector>& density);
  
private:

  void MakeTable(const string_array& inputPara);

  static void SpanWagnerCO2Density(const real64 &T, const real64 &P, real64 &rho, real64 (*f)(const real64 &x1, const real64 &x2, const real64 &x3));
  

  TableFunctionPtr m_CO2DensityTable;
  localIndex m_CO2Index;

};

}

}

#endif
