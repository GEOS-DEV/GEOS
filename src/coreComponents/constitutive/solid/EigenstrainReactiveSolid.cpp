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
 * @file EigenstrainReactiveSolid.cpp
 */

#include "EigenstrainReactiveSolid.hpp"
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
EigenstrainReactiveSolid< SOLID_TYPE, PERM_TYPE, SURFACE_AREA_TYPE >::EigenstrainReactiveSolid( string const & name, Group * const parent ):
  CoupledSolid< SOLID_TYPE, ReactivePorosityBase, PERM_TYPE >( name, parent )
{}

template< typename SOLID_TYPE,
          typename PERM_TYPE,
          typename SURFACE_AREA_TYPE >
void EigenstrainReactiveSolid< SOLID_TYPE, PERM_TYPE, SURFACE_AREA_TYPE >::initializeState() const
{
  CoupledSolid< SOLID_TYPE, ReactivePorosityBase, PERM_TYPE >::initializeState();
}

// Register all EigenstrainReactiveSolid model types.
typedef EigenstrainReactiveSolid< ElasticIsotropic, ConstantPermeability, ConstantSurfaceArea > EigenStrainReactiveElasticIsotropicConstant;
typedef EigenstrainReactiveSolid< ElasticIsotropic, CarmanKozenyPermeability, ConstantSurfaceArea > EigenStrainReactiveElasticIsotropicCK;
typedef EigenstrainReactiveSolid< ElasticIsotropic, ConstantPermeability, PowerLawSurfaceArea > EigenStrainReactiveElasticIsotropicConstantPowerLawSurfaceArea;
typedef EigenstrainReactiveSolid< ElasticIsotropic, CarmanKozenyPermeability, PowerLawSurfaceArea > EigenStrainReactiveElasticIsotropicCKPowerLawSurfaceArea;
typedef EigenstrainReactiveSolid< ElasticIsotropic, ConstantPermeability, SubstrateCoverageSurfaceArea > EigenStrainReactiveElasticIsotropicConstantSubstrateCoverageSurfaceArea;
typedef EigenstrainReactiveSolid< ElasticIsotropic, CarmanKozenyPermeability, SubstrateCoverageSurfaceArea > EigenStrainReactiveElasticIsotropicCKSubstrateCoverageSurfaceArea;


REGISTER_CATALOG_ENTRY( ConstitutiveBase, EigenStrainReactiveElasticIsotropicConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, EigenStrainReactiveElasticIsotropicCK, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, EigenStrainReactiveElasticIsotropicConstantPowerLawSurfaceArea, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, EigenStrainReactiveElasticIsotropicCKPowerLawSurfaceArea, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, EigenStrainReactiveElasticIsotropicConstantSubstrateCoverageSurfaceArea, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, EigenStrainReactiveElasticIsotropicCKSubstrateCoverageSurfaceArea, string const &, Group * const )


}
} /* namespace geos */
