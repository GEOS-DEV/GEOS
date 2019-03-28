
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

#include "FlashModelBase.hpp"

namespace geosx
{

namespace PVTProps
{

class CO2SolubilityFunction : public FlashModelBase
{
public:

  CO2SolubilityFunction( const string_array& inputPara,
                         const string_array& phaseNames,
                         const string_array& componentNames,
                         real64_array const & componentMolarWeight);

  ~CO2SolubilityFunction() override
  {}

  static constexpr auto m_catalogName = "CO2Solubility";
  static string CatalogName()                    { return m_catalogName; }
  virtual string GetCatalogName() override final { return CatalogName(); }

  virtual void Partition( const EvalVarArgs& pressure,
                          const EvalVarArgs& temperature,
                          const array1dT<EvalVarArgs>& compFraction,
                          array1dT<EvalVarArgs>& phaseFraction,
                          array1dT<array1dT<EvalVarArgs> >& phaseCompFraction) const override;

private:

  void MakeTable(const string_array& inputPara);

  void CO2Solubility(const double &T, const double &P, double &V_r, double (*f)(const double &x1, const double &x2, const double &x3));

  void CalculateCO2Solubility(const real64_vector& pressure, const real64_vector& temperature, const real64& salinity, array1dT<real64_vector>& solubiltiy);  
  
  TableFunctionPtr m_CO2SolubilityTable;
  localIndex m_CO2Index;
  localIndex m_waterIndex;
  localIndex m_phaseGasIndex;
  localIndex m_phaseLiquidIndex;
};

}

}

#endif
