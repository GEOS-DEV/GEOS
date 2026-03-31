/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FieldSpecificationFactory.hpp
 */

#ifndef GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONFACTORY_HPP
#define GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONFACTORY_HPP


#include "common/DataTypes.hpp"
#include "dataRepository/Group.hpp"
#include "FieldSpecificationABC.hpp"

namespace geos
{

/**
 * @class FieldSpecificationFactory
 */
class FieldSpecificationFactory
{
public:

  /// @brief Generate FieldSpecifications based on the given "higher-level"
  ///        specification
  /// @param specification The higher-level specification used as a blueprint
  ///                      to create FieldSpecification
  /// @param manager The parent to store the created FieldSpecifications
  virtual void generate( FieldSpecificationABC const & specification,
                         dataRepository::Group & manager ) const = 0;

  /// @return The key that represents the element this factory is about.
  /// Purpose: link the factory to the specification it uses.
  virtual string const getKey() const = 0;

};

}


#endif //GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONFACTORY_HPP