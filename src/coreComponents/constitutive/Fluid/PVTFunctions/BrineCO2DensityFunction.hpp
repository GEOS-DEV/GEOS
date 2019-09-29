
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
 * @file BrineCO2DensityFunction.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_BRINECO2DENSITYFUNCTION_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_BRINECO2DENSITYFUNCTION_HPP

#include "PVTFunctionBase.hpp"

namespace geosx
{

namespace PVTProps
{

class BrineCO2DensityFunction : public PVTFunction
{
public:

  BrineCO2DensityFunction( string_array const & inputPara,
                           string_array const & componentNames,
                           real64_array const & componentMolarWeight);
  ~BrineCO2DensityFunction() override {}


  static constexpr auto m_catalogName = "BrineCO2Density";
  static string CatalogName()                    { return m_catalogName; }
  virtual string GetCatalogName() override final { return CatalogName(); }


  virtual PVTFuncType FunctionType() const override
      {
    return PVTFuncType::DENSITY;

      }

  virtual void Evaluation( EvalVarArgs const & pressure,
                           EvalVarArgs const & temperature,
                           arraySlice1d<EvalVarArgs const> const & phaseComposition,
                           EvalVarArgs & value,
                           bool useMass = 0) const override;


private:

  void MakeTable(string_array const & inputPara);

  void CalculateBrineDensity(real64_const_array const & pressure, real64_const_array const & temperature, real64 const & salinity, real64_array2d const & density);
  
  TableFunctionPtr m_BrineDensityTable;
  localIndex m_CO2Index;
  localIndex m_waterIndex;

};

}

}

#endif
