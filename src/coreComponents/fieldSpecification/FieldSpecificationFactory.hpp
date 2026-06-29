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


#include "common/TypeDispatch.hpp"
#include "dataRepository/Group.hpp"

namespace geos
{

class PermeabilitySpecification;

/**
 * @brief List of high-level field specifications to expand into FieldSpecification objects.
 */
using ExpandableSpecTypes = types::TypeList< PermeabilitySpecification >;

/**
 * @brief Generate FieldSpecifications based on the given "higher-level" specification
 * @tparam SPEC_TYPE The type of the high-level specification
 * @param specification The high-level specification used as a blueprint
 *                      to create FieldSpecification
 * @param manager The parent to store the created FieldSpecifications
 */
template< typename SPEC_TYPE >
void generateFieldSpecifications( SPEC_TYPE const & specification, dataRepository::Group & manager ) = delete;

/**
 * @brief Expand high-level field specifications
 * @param manager The manager group
 * Expand all field specificiation listed in ExpandableSpecTypes
 */
void expandFieldSpecifications( dataRepository::Group & manager );

} // namespace geos


#endif //GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONFACTORY_HPP
