/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file NumericalMethodsManager.hpp
 */

#ifndef GEOSX_MANAGERS_NUMERICALMETHODSMANAGER_HPP_
#define GEOSX_MANAGERS_NUMERICALMETHODSMANAGER_HPP_

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
  NumericalMethodsManager( string const & name, Group * const parent );
  virtual ~NumericalMethodsManager() override;

  virtual Group * CreateChild( string const & childKey, string const & childName ) override;

  dataRepository::Group const * FindNumericalMethodByName( string const & name ) const;


};

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_NUMERICALMETHODSMANAGER_HPP_ */
