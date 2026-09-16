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
#include "constitutive/permeability/DamageCloggingPermeability.hpp"
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
  
  this->registerWrapper( "surfaceAreaDamageExponent", &m_surfaceAreaDamageExponent ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Exponent confining the reactive surface area to damaged cells: A = A0 (1-theta)^(2/3) d^s. "
                    "The default 0 leaves the area ungated, since d^0 = 1 even in an intact cell." );

  this->registerWrapper( "fluidModelName", &m_fluidModelName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the fluid constitutive model. Its (constant) compressibility "
                    "is handed to the porosity model for the pore-mineral-pressure / porosity coupling." );
}

template< typename SOLID_TYPE,
          typename PERM_TYPE,
          typename DIFF_TYPE >
void PorousReactiveSolid< SOLID_TYPE, PERM_TYPE, DIFF_TYPE >::initializeState() const
{
  CoupledSolid< SOLID_TYPE, BiotReactivePorosity, PERM_TYPE >::initializeState();
}

template< typename SOLID_TYPE,
          typename PERM_TYPE,
          typename DIFF_TYPE >
void PorousReactiveSolid< SOLID_TYPE, PERM_TYPE, DIFF_TYPE >::initializePreSubGroups()
{
  CoupledSolid< SOLID_TYPE, BiotReactivePorosity, PERM_TYPE >::initializePreSubGroups();

  // Read the fluid's (constant) compressibility and hand it to the porosity model, which needs the
  // fluid bulk modulus for the pore-mineral-pressure / porosity coupling but has no direct handle on
  // the fluid. The fluid is looked up through its CompressibleSinglePhaseFluid base, which also covers
  // the reactive and thermal variants.
  CompressibleSinglePhaseFluid const & fluid =
    this->getParent().template getGroup< CompressibleSinglePhaseFluid >( m_fluidModelName );
  GEOS_THROW_IF_LE_MSG( fluid.compressibility(), 0.0,
                        GEOS_FMT( "fluid model '{}' must have a positive compressibility, "
                                  "since the porosity model uses the fluid bulk modulus 1/compressibility",
                                  m_fluidModelName ),
                        InputError, this->getDataContext() );
  BiotReactivePorosity & porosity =
    dynamicCast< BiotReactivePorosity & >( this->getBasePorosityModel() );
  porosity.setFluidCompressibility( fluid.compressibility() );
}

// Register all PorousReactiveSolid model types.
typedef PorousReactiveSolid< ElasticIsotropic, ConstantPermeability > PorousReactiveElasticIsotropicConstant;
typedef PorousReactiveSolid< Damage< ElasticIsotropic >, ConstantPermeability > PorousReactiveDamageConstant;
typedef PorousReactiveSolid< DamageSpectral< ElasticIsotropic >, ConstantPermeability > PorousReactiveDamageSpectralConstant;
typedef PorousReactiveSolid< DamageVolDev< ElasticIsotropic >, ConstantPermeability > PorousReactiveDamageVolDevConstant;

typedef PorousReactiveSolid< ElasticIsotropic, CarmanKozenyPermeability > PorousReactiveElasticIsotropicCK;

typedef PorousReactiveSolid< Damage< ElasticIsotropic >, DamagePermeability > PorousReactiveDamageDamagePermeability;
typedef PorousReactiveSolid< DamageSpectral< ElasticIsotropic >, DamagePermeability > PorousReactiveDamageSpectralDamagePermeability;
typedef PorousReactiveSolid< DamageVolDev< ElasticIsotropic >, DamagePermeability > PorousReactiveDamageVolDevDamagePermeability;

// Damage solid + mineral-clogging damage permeability
typedef PorousReactiveSolid< Damage< ElasticIsotropic >, DamageCloggingPermeability > PorousReactiveDamageDamageCloggingPermeability;
typedef PorousReactiveSolid< DamageSpectral< ElasticIsotropic >, DamageCloggingPermeability > PorousReactiveDamageSpectralDamageCloggingPermeability;

// Damage solid + damage permeability + damage diffusion
typedef PorousReactiveSolid< Damage< ElasticIsotropic >, DamagePermeability, DamageDiffusion > PorousReactiveDamageDamagePermeabilityDamageDiffusion;


REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveElasticIsotropicConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveDamageConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveDamageSpectralConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveDamageVolDevConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveElasticIsotropicCK, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveDamageDamagePermeability, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveDamageSpectralDamagePermeability, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveDamageDamageCloggingPermeability, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveDamageSpectralDamageCloggingPermeability, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveDamageVolDevDamagePermeability, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveDamageDamagePermeabilityDamageDiffusion, string const &, Group * const )


}
} /* namespace geos */
