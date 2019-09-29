
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

class CO2SolubilityFunction : public FlashModel
{
public:

  CO2SolubilityFunction( string_array const & inputPara,
                         string_array const & phaseNames,
                         string_array const & componentNames,
                         real64_array const & componentMolarWeight);

  ~CO2SolubilityFunction() override
  {}

  static constexpr auto m_catalogName = "CO2Solubility";
  static string CatalogName()                    { return m_catalogName; }
  virtual string GetCatalogName() override final { return CatalogName(); }

  virtual void Partition( EvalVarArgs const & pressure,
                          EvalVarArgs const & temperature,
                          arraySlice1d<EvalVarArgs const> const & compFraction,
                          arraySlice1d<EvalVarArgs> const & phaseFraction,
                          arraySlice2d<EvalVarArgs> const & phaseCompFraction) const override;

private:

  void MakeTable(const string_array& inputPara);

  TableFunctionPtr m_CO2SolubilityTable;
  localIndex m_CO2Index;
  localIndex m_waterIndex;
  localIndex m_phaseGasIndex;
  localIndex m_phaseLiquidIndex;
};

}

}

#endif
