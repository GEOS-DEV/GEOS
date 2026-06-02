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
 * @file PorousReactiveSolid.cpp
 */

#include "PorousReactiveSolid.hpp"
#include "ElasticIsotropic.hpp"
#include "constitutive/permeability/ConstantPermeability.hpp"
#include "constitutive/permeability/CarmanKozenyPermeability.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

template< typename SOLID_TYPE,
          typename PERM_TYPE >
PorousReactiveSolid< SOLID_TYPE, PERM_TYPE >::PorousReactiveSolid( string const & name, Group * const parent ):
  CoupledSolid< SOLID_TYPE, BiotReactivePorosity, PERM_TYPE >( name, parent )
{}

template< typename SOLID_TYPE,
          typename PERM_TYPE >
void PorousReactiveSolid< SOLID_TYPE, PERM_TYPE >::initializeState() const
{
  CoupledSolid< SOLID_TYPE, BiotReactivePorosity, PERM_TYPE >::initializeState();
}

// Register all PorousReactiveSolid model types.
typedef PorousReactiveSolid< ElasticIsotropic, ConstantPermeability > PorousReactiveElasticIsotropicConstant;
typedef PorousReactiveSolid< ElasticIsotropic, CarmanKozenyPermeability > PorousReactiveElasticIsotropicCK;


REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveElasticIsotropicConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveElasticIsotropicCK, string const &, Group * const )


}
} /* namespace geos */
