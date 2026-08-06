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
#include "constitutive/surfaceArea/ConstantSurfaceArea.hpp"
#include "constitutive/surfaceArea/PowerLawSurfaceArea.hpp"
#include "constitutive/surfaceArea/SubstrateCoverageSurfaceArea.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

template< typename SOLID_TYPE,
          typename PERM_TYPE,
          typename SURFACE_AREA_TYPE >
PorousReactiveSolid< SOLID_TYPE, PERM_TYPE, SURFACE_AREA_TYPE >::PorousReactiveSolid( string const & name, Group * const parent ):
  CoupledSolid< SOLID_TYPE, BiotReactivePorosity, PERM_TYPE >( name, parent )
{
  this->registerWrapper( "fluidModelName", &m_fluidModelName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the fluid constitutive model. When set, its (constant) compressibility "
                    "is handed to the porosity model for the pore-mineral-pressure / porosity coupling." );
}

template< typename SOLID_TYPE,
          typename PERM_TYPE,
          typename SURFACE_AREA_TYPE >
void PorousReactiveSolid< SOLID_TYPE, PERM_TYPE, SURFACE_AREA_TYPE >::initializeState() const
{
  CoupledSolid< SOLID_TYPE, BiotReactivePorosity, PERM_TYPE >::initializeState();
}

template< typename SOLID_TYPE,
          typename PERM_TYPE,
          typename SURFACE_AREA_TYPE >
void PorousReactiveSolid< SOLID_TYPE, PERM_TYPE, SURFACE_AREA_TYPE >::initializePreSubGroups()
{
  CoupledSolid< SOLID_TYPE, BiotReactivePorosity, PERM_TYPE >::initializePreSubGroups();

  // If a fluid model is specified, read its (constant) compressibility and hand it to the porosity
  // model. The porosity model needs the fluid bulk modulus for the pore-mineral-pressure / porosity
  // coupling but has no direct handle on the fluid. The fluid is looked up through its
  // CompressibleSinglePhaseFluid base, which also covers the reactive and thermal variants.
  if( !m_fluidModelName.empty() )
  {
    CompressibleSinglePhaseFluid const & fluid =
      this->getParent().template getGroup< CompressibleSinglePhaseFluid >( m_fluidModelName );
    BiotReactivePorosity & porosity =
      dynamicCast< BiotReactivePorosity & >( this->getBasePorosityModel() );
    porosity.setFluidCompressibility( fluid.compressibility() );
  }
}

// Register all PorousReactiveSolid model types.
typedef PorousReactiveSolid< ElasticIsotropic, ConstantPermeability, ConstantSurfaceArea > PorousReactiveElasticIsotropicConstant;
typedef PorousReactiveSolid< ElasticIsotropic, CarmanKozenyPermeability, ConstantSurfaceArea > PorousReactiveElasticIsotropicCK;
typedef PorousReactiveSolid< ElasticIsotropic, ConstantPermeability, PowerLawSurfaceArea > PorousReactiveElasticIsotropicConstantPowerLawSurfaceArea;
typedef PorousReactiveSolid< ElasticIsotropic, CarmanKozenyPermeability, PowerLawSurfaceArea > PorousReactiveElasticIsotropicCKPowerLawSurfaceArea;
typedef PorousReactiveSolid< ElasticIsotropic, ConstantPermeability, SubstrateCoverageSurfaceArea > PorousReactiveElasticIsotropicConstantSubstrateCoverageSurfaceArea;
typedef PorousReactiveSolid< ElasticIsotropic, CarmanKozenyPermeability, SubstrateCoverageSurfaceArea > PorousReactiveElasticIsotropicCKSubstrateCoverageSurfaceArea;


REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveElasticIsotropicConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveElasticIsotropicCK, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveElasticIsotropicConstantPowerLawSurfaceArea, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveElasticIsotropicCKPowerLawSurfaceArea, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveElasticIsotropicConstantSubstrateCoverageSurfaceArea, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousReactiveElasticIsotropicCKSubstrateCoverageSurfaceArea, string const &, Group * const )


}
} /* namespace geos */
