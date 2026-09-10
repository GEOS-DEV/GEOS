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
 * @file ReactiveSolid.cpp
 */

#include "ReactiveSolid.hpp"
#include "porosity/ReactivePorosityBase.hpp"
#include "constitutive/permeability/ConstantPermeability.hpp"
#include "constitutive/permeability/CarmanKozenyPermeability.hpp"
#include "constitutive/permeability/PressurePermeability.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

template< typename PORO_TYPE,
          typename PERM_TYPE >
ReactiveSolid< PORO_TYPE, PERM_TYPE >::ReactiveSolid( string const & name, Group * const parent ):
  CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >( name, parent )
{}

template< typename PORO_TYPE,
          typename PERM_TYPE >
ReactiveSolid< PORO_TYPE, PERM_TYPE >::~ReactiveSolid() = default;

// Register all ReactiveSolid model types.
typedef ReactiveSolid< ReactivePorosityBase, ConstantPermeability > ReactiveRockConstant;
typedef ReactiveSolid< ReactivePorosityBase, CarmanKozenyPermeability > ReactiveRockCK;
typedef ReactiveSolid< ReactivePorosityBase, PressurePermeability > ReactiveRockPressurePerm;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveRockConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveRockCK, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveRockPressurePerm, string const &, Group * const )

}
} /* namespace geos */
