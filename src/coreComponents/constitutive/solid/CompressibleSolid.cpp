/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file CompressibleSolid.cpp
 */

#include "CompressibleSolid.hpp"
#include "porosity/PressurePorosity.hpp"
#include "constitutive/permeability/ConstantPermeability.hpp"
#include "constitutive/permeability/ConstantPlusParallelPlatesPermeability.hpp"
#include "constitutive/permeability/CarmanKozenyPermeability.hpp"
#include "constitutive/permeability/ParallelPlatesPermeability.hpp"
#include "constitutive/permeability/SlipDependentPermeability.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

template< typename PORO_TYPE,
          typename PERM_TYPE >
CompressibleSolid< PORO_TYPE, PERM_TYPE >::CompressibleSolid( string const & name, Group * const parent ):
  CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >( name, parent )
{}

template< typename PORO_TYPE,
          typename PERM_TYPE >
CompressibleSolid< PORO_TYPE, PERM_TYPE >::~CompressibleSolid() = default;

// Register all CompressibleSolid model types.
typedef CompressibleSolid< PressurePorosity, ConstantPermeability > CompressibleRockConstant;
typedef CompressibleSolid< PressurePorosity, CarmanKozenyPermeability > CompressibleRockCK;
typedef CompressibleSolid< PressurePorosity, ParallelPlatesPermeability > FractureRock;
typedef CompressibleSolid< PressurePorosity, ConstantPlusParallelPlatesPermeability > FractureRockPlus;
typedef CompressibleSolid< PressurePorosity, SlipDependentPermeability > Fault;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompressibleRockConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompressibleRockCK, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, FractureRock, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, FractureRockPlus, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, Fault, string const &, Group * const )


}
} /* namespace geosx */
