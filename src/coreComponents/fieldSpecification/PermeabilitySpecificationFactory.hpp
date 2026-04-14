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
 * @file PermeabilitySpecificationFactory.hpp
 */

#ifndef GEOS_FIELDSPECIFICATION_PERMEABILITYSPECIFICATIONFACTORY_HPP
#define GEOS_FIELDSPECIFICATION_PERMEABILITYSPECIFICATIONFACTORY_HPP


#include "common/DataTypes.hpp"
#include "dataRepository/Group.hpp"
#include "PermeabilitySpecification.hpp"
#include "FieldSpecificationABC.hpp"
#include "FieldSpecificationBase.hpp"
#include "FieldSpecificationFactory.hpp"

namespace geos
{

/**
 * @class PermeabilitySpecificationFactory
 *
 * @copydoc geos::FieldSpecificationFactory
 * Field specification factory implementation for the PermeabilitySpecification
 */
class PermeabilitySpecificationFactory : public FieldSpecificationFactory
{
public:

  /// @copydoc geos::FieldSpecificationFactory::generate()
  void generate( FieldSpecificationABC const & specification,
                 dataRepository::Group & manager ) const;

  /// @copydoc geos::FieldSpecificationFactory::getKey()
  string const getKey() const
  {
    return PermeabilitySpecification::catalogName();
  }
};

}

#endif //GEOS_FIELDSPECIFICATION_PERMEABILITYSPECIFICATIONFACTORY_HPP
