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
#include "constitutive/permeability/PowerLawPermeability.hpp"
#include "constitutive/permeability/PressurePermeability.hpp"
#include "constitutive/surfaceArea/ConstantSurfaceArea.hpp"
#include "constitutive/surfaceArea/PowerLawSurfaceArea.hpp"
#include "constitutive/surfaceArea/SubstrateCoverageSurfaceArea.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

template< typename PORO_TYPE,
          typename PERM_TYPE,
          typename SURFACE_AREA_TYPE >
ReactiveSolid< PORO_TYPE, PERM_TYPE, SURFACE_AREA_TYPE >::ReactiveSolid( string const & name, Group * const parent ):
  CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >( name, parent )
{}

template< typename PORO_TYPE,
          typename PERM_TYPE,
          typename SURFACE_AREA_TYPE >
ReactiveSolid< PORO_TYPE, PERM_TYPE, SURFACE_AREA_TYPE >::~ReactiveSolid() = default;

template< typename PORO_TYPE,
          typename PERM_TYPE,
          typename SURFACE_AREA_TYPE >
void ReactiveSolid< PORO_TYPE, PERM_TYPE, SURFACE_AREA_TYPE >::initializePreSubGroups()
{
  CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >::initializePreSubGroups();

  GEOS_THROW_IF( !this->hasSurfaceAreaModel(),
                 GEOS_FMT( "{}: the reactive solid requires a surface area model, please set the attribute '{}'",
                           this->getFullName(),
                           CoupledSolidBase::viewKeyStruct::surfaceAreaModelNameString() ),
                 InputError );

  SURFACE_AREA_TYPE const & surfaceAreaModel = getSurfaceAreaModel();

  GEOS_ERROR_IF( SURFACE_AREA_TYPE::catalogName() != surfaceAreaModel.getCatalogName(),
                 GEOS_FMT( " The reactive solid expects a surface area model of type {} but the specified surface area model \"{}\" is of type {}",
                           SURFACE_AREA_TYPE::catalogName(),
                           this->m_surfaceAreaModelName,
                           surfaceAreaModel.getCatalogName() ),
                 this->getDataContext() );

  GEOS_THROW_IF_NE_MSG( surfaceAreaModel.numKineticReactions(), getPorosityModel().numKineticReactions(),
                        GEOS_FMT( "{}: mismatch in the number of kinetic reactions between the surface area model {} and the porosity model {}",
                                  this->getFullName(),
                                  surfaceAreaModel.getDataContext(),
                                  getPorosityModel().getDataContext() ),
                        InputError );
}

// Register all ReactiveSolid model types.
typedef ReactiveSolid< ReactivePorosityBase, ConstantPermeability, ConstantSurfaceArea > ReactiveRockConstant;
typedef ReactiveSolid< ReactivePorosityBase, CarmanKozenyPermeability, ConstantSurfaceArea > ReactiveRockCK;
typedef ReactiveSolid< ReactivePorosityBase, PressurePermeability, ConstantSurfaceArea > ReactiveRockPressurePerm;
typedef ReactiveSolid< ReactivePorosityBase, PowerLawPermeability, ConstantSurfaceArea > ReactiveRockPowerLawPerm;
typedef ReactiveSolid< ReactivePorosityBase, ConstantPermeability, PowerLawSurfaceArea > ReactiveRockConstantPowerLawSurfaceArea;
typedef ReactiveSolid< ReactivePorosityBase, CarmanKozenyPermeability, PowerLawSurfaceArea > ReactiveRockCKPowerLawSurfaceArea;
typedef ReactiveSolid< ReactivePorosityBase, PressurePermeability, PowerLawSurfaceArea > ReactiveRockPressurePermPowerLawSurfaceArea;
typedef ReactiveSolid< ReactivePorosityBase, PowerLawPermeability, PowerLawSurfaceArea > ReactiveRockPowerLawPermPowerLawSurfaceArea;
typedef ReactiveSolid< ReactivePorosityBase, ConstantPermeability, SubstrateCoverageSurfaceArea > ReactiveRockConstantSubstrateCoverageSurfaceArea;
typedef ReactiveSolid< ReactivePorosityBase, CarmanKozenyPermeability, SubstrateCoverageSurfaceArea > ReactiveRockCKSubstrateCoverageSurfaceArea;
typedef ReactiveSolid< ReactivePorosityBase, PressurePermeability, SubstrateCoverageSurfaceArea > ReactiveRockPressurePermSubstrateCoverageSurfaceArea;
typedef ReactiveSolid< ReactivePorosityBase, PowerLawPermeability, SubstrateCoverageSurfaceArea > ReactiveRockPowerLawPermSubstrateCoverageSurfaceArea;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveRockConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveRockCK, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveRockPressurePerm, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveRockPowerLawPerm, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveRockConstantPowerLawSurfaceArea, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveRockCKPowerLawSurfaceArea, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveRockPressurePermPowerLawSurfaceArea, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveRockPowerLawPermPowerLawSurfaceArea, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveRockConstantSubstrateCoverageSurfaceArea, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveRockCKSubstrateCoverageSurfaceArea, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveRockPressurePermSubstrateCoverageSurfaceArea, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveRockPowerLawPermSubstrateCoverageSurfaceArea, string const &, Group * const )

}
} /* namespace geos */
