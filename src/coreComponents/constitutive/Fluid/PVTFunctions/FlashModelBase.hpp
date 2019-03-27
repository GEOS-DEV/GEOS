
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
 * @file FlashModelBase.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_FLASHMODELBASE_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_FLASHMODELBASE_HPP

#include "constitutive/Fluid/PVTFunctions/UtilityFunctions.hpp"
#include "codingUtilities/StringUtilities.hpp"

namespace geosx
{

namespace PVTProps
{


class FlashModelBase
{
public:

  FlashModelBase( string const & name,
                  const string_array& componentNames,
                  const real64_array& componentMolarWeight):
    m_modelName(name),
    m_componentNames(componentNames),
    m_componentMolarWeight(componentMolarWeight)
  {}

  virtual ~FlashModelBase(){}

  using CatalogInterface = cxx_utilities::CatalogInterface< FlashModelBase, string_array const &,
                                                                            string_array const &,
                                                                            string_array const &,
                                                                            real64_array const & >;
  static typename CatalogInterface::CatalogType& GetCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }
  virtual string GetCatalogName() = 0;


  string const & FlashModelName() const
  {
    return m_modelName;
  }

  //partition
  //input: P, T, totalCompFraction
  //output: phaseFraction, phaseCompFraction

  virtual void Partition( const EvalVarArgs& pressure,
                          const EvalVarArgs& temperature,
                          const array1dT<EvalVarArgs>& compFraction,
                          array1dT<EvalVarArgs>& phaseFraction,
                          array1dT<array1dT<EvalVarArgs> >& phaseCompFraction) const = 0;

protected:
  string m_modelName;
  string_array m_componentNames;
  real64_array m_componentMolarWeight;

};

typedef std::unique_ptr<FlashModelBase> FlashModel;

void CalculateCO2Solubility(const real64_vector& pressure, const real64_vector& temperature, const real64& salinity, array1dT<real64_vector>& solubiltiy);

}

}

#endif
