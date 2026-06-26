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
#include "Damage.hpp"
#include "DamageSpectral.hpp"
#include "DamageVolDev.hpp"
#include "constitutive/permeability/ConstantPermeability.hpp"
#include "constitutive/permeability/CarmanKozenyPermeability.hpp"
#include "constitutive/permeability/DamagePermeability.hpp"
#include "constitutive/diffusion/DamageDiffusion.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

template< typename SOLID_TYPE,
          typename PERM_TYPE,
          typename DIFF_TYPE >
PorousReactiveSolid< SOLID_TYPE, PERM_TYPE, DIFF_TYPE >::PorousReactiveSolid( string const & name, Group * const parent ):
  CoupledSolid< SOLID_TYPE, BiotReactivePorosity, PERM_TYPE >( name, parent )
{
  if constexpr( !std::is_same_v< DIFF_TYPE, NoDiffusion > )
  {
    this->registerWrapper( "diffusionModelName", &m_diffusionModelName ).
      setInputFlag( InputFlags::REQUIRED ).
      setDescription( "Name of the diffusion constitutive model" );
  }
}

template< typename SOLID_TYPE,
          typename PERM_TYPE,
          typename DIFF_TYPE >
void PorousReactiveSolid< SOLID_TYPE, PERM_TYPE, DIFF_TYPE >::initializeState() const
{
  CoupledSolid< SOLID_TYPE, BiotReactivePorosity, PERM_TYPE >::initializeState();
}

// Register all PorousReactiveSolid model types.
typedef PorousReactiveSolid< ElasticIsotropic, ConstantPermeability > PorousReactiveElasticIsotropicConstant;
typedef PorousReactiveSolid< Damage< ElasticIsotropic >, ConstantPermeability > PorousReactiveDamageConstant;
// typedef PorousReactiveSolid< DamageSpectral< ElasticIsotropic >, ConstantPermeability > PorousReactiveDamageSpectralConstant;
typedef PorousReactiveSolid< DamageVolDev< ElasticIsotropic >, ConstantPermeability > PorousReactiveDamageVolDevConstant;

typedef PorousReactiveSolid< ElasticIsotropic, CarmanKozenyPermeability > PorousReactiveElasticIsotropicCK;

typedef PorousReactiveSolid< Damage< ElasticIsotropic >, DamagePermeability > PorousReactiveDamageDamagePermeability;
// typedef PorousReactiveSolid< DamageSpectral< ElasticIsotropic >, DamagePermeability > PorousReactiveDamageSpectralDamagePermeability;
typedef PorousReactiveSolid< DamageVolDev< ElasticIsotropic >, DamagePermeability > PorousReactiveDamageVolDevDamagePermeability;

// Damage solid + damage permeability + damage diffusion
typedef PorousReactiveSolid< Damage< ElasticIsotropic >, DamagePermeability, DamageDiffusion > PorousReactiveDamageDamagePermeabilityDamageDiffusion;


REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveElasticIsotropicConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveDamageConstant, string const &, Group * const )
// REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveDamageSpectralConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveDamageVolDevConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveElasticIsotropicCK, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveDamageDamagePermeability, string const &, Group * const )
// REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveDamageSpectralDamagePermeability, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveDamageVolDevDamagePermeability, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveDamageDamagePermeabilityDamageDiffusion, string const &, Group * const )


}
} /* namespace geos */
