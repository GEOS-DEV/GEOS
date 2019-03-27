
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

#include "PVTFunctionBase.hpp"

namespace geosx
{

namespace PVTProps
{

class BrineViscosityFunction : public PVTFunctionBase
{
public:


  BrineViscosityFunction( string_array const & inputPara,
                          string_array const & componentNames,
                          real64_array const & componentMolarWeight );
  ~BrineViscosityFunction() override
  {}

  static constexpr auto m_catalogName = "BrineViscosity";
  static string CatalogName()                    { return m_catalogName; }
  virtual string GetCatalogName() override final { return CatalogName(); }

  virtual PVTFUNCTYPE FunctionType() const override
  {
    return PVTFUNCTYPE::VISCOSITY;

  }

  virtual void Evaluation( const EvalVarArgs& pressure,
                           const EvalVarArgs& temperature,
                           const array1dT<EvalVarArgs>& phaseComposition,
                           EvalVarArgs& value, bool useMass = 0) const override;

private:

  void MakeCoef(string_array const & inputPara);

  real64 m_coef0;
  real64 m_coef1;

};

}

}  
#endif
