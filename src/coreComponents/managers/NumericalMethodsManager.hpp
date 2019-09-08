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

/*
 * NumericalMethodsManager.hpp
 *
 *  Created on: Apr 18, 2017
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MANAGERS_NUMERICALMETHODSMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_MANAGERS_NUMERICALMETHODSMANAGER_HPP_

#include "dataRepository/Group.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const numericalMethodsManager = "NumericalMethods";
string const basisFunctions = "BasisFunctions";
string const quadratureRules = "QuadratureRules";
string const finiteElementDiscretizations = "FiniteElements";
string const finiteVolumeManager = "FiniteVolume";
}
}


class NumericalMethodsManager : public dataRepository::Group
{
public:
  NumericalMethodsManager() = delete;
  NumericalMethodsManager(string const & name, Group * const parent);
  virtual ~NumericalMethodsManager() override;

  virtual Group * CreateChild( string const & childKey, string const & childName ) override;

  dataRepository::Group const * FindNumericalMethodByName(string const & name) const;


};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_NUMERICALMETHODSMANAGER_HPP_ */
