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
 * @file PorousSolid.cpp
 */

#include "PorousSolid.hpp"
#include "ElasticIsotropic.hpp"
#include "ElasticIsotropicPressureDependent.hpp"
#include "ElasticTransverseIsotropic.hpp"
#include "ElasticOrthotropic.hpp"
#include "DelftEgg.hpp"
#include "DruckerPrager.hpp"
#include "DruckerPragerExtended.hpp"
#include "ModifiedCamClay.hpp"
#include "DuvautLionsSolid.hpp"
#include "constitutive/permeability/ConstantPermeability.hpp"
#include "constitutive/permeability/CarmanKozenyPermeability.hpp"
#include "constitutive/permeability/StrainDependentPermeability.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

template< typename SOLID_TYPE,
          typename PERM_TYPE >
PorousSolid< SOLID_TYPE, PERM_TYPE >::PorousSolid( string const & name, Group * const parent ):
  CoupledSolid< SOLID_TYPE, BiotPorosity, PERM_TYPE >( name, parent )
{}

template< typename SOLID_TYPE,
          typename PERM_TYPE >
void PorousSolid< SOLID_TYPE, PERM_TYPE >::initializeState() const
{
  CoupledSolid< SOLID_TYPE, BiotPorosity, PERM_TYPE >::initializeState();
}

// Register all PorousSolid model types.
typedef PorousSolid< ElasticIsotropic, ConstantPermeability > PorousElasticIsotropicConstant;
typedef PorousSolid< ElasticTransverseIsotropic, ConstantPermeability > PorousElasticTransverseIsotropicConstant;
typedef PorousSolid< ElasticOrthotropic, ConstantPermeability > PorousElasticOrthotropicConstant;
typedef PorousSolid< DelftEgg, ConstantPermeability > PorousDelftEggConstant;
typedef PorousSolid< DruckerPrager, ConstantPermeability > PorousDruckerPragerConstant;
typedef PorousSolid< DruckerPragerExtended, ConstantPermeability > PorousDruckerPragerExtendedConstant;
typedef PorousSolid< DuvautLionsSolid< DruckerPrager >, ConstantPermeability > PorousViscoDruckerPragerConstant;
typedef PorousSolid< DuvautLionsSolid< DruckerPragerExtended >, ConstantPermeability > PorousViscoDruckerPragerExtendedConstant;
typedef PorousSolid< DuvautLionsSolid< ModifiedCamClay >, ConstantPermeability > PorousViscoModifiedCamClayConstant;
typedef PorousSolid< ModifiedCamClay, ConstantPermeability > PorousModifiedCamClayConstant;
typedef PorousSolid< ElasticIsotropic, CarmanKozenyPermeability > PorousElasticIsotropicCK;
typedef PorousSolid< ElasticTransverseIsotropic, CarmanKozenyPermeability > PorousElasticTransverseIsotropicCK;
typedef PorousSolid< ElasticOrthotropic, CarmanKozenyPermeability > PorousElasticOrthotropicCK;
typedef PorousSolid< DelftEgg, CarmanKozenyPermeability > PorousDelftEggCK;
typedef PorousSolid< DruckerPrager, CarmanKozenyPermeability > PorousDruckerPragerCK;
typedef PorousSolid< DruckerPragerExtended, CarmanKozenyPermeability > PorousDruckerPragerExtendedCK;
typedef PorousSolid< DuvautLionsSolid< DruckerPrager >, CarmanKozenyPermeability > PorousViscoDruckerPragerCK;
typedef PorousSolid< DuvautLionsSolid< DruckerPragerExtended >, CarmanKozenyPermeability > PorousViscoDruckerPragerExtendedCK;
typedef PorousSolid< DuvautLionsSolid< ModifiedCamClay >, CarmanKozenyPermeability > PorousViscoModifiedCamClayCK;
typedef PorousSolid< ModifiedCamClay, CarmanKozenyPermeability > PorousModifiedCamClayCK;
typedef PorousSolid< ElasticIsotropic, StrainDependentPermeability > PorousElasticIsotropicSD;
typedef PorousSolid< ElasticIsotropicPressureDependent, StrainDependentPermeability > PorousElasticIsotropicPressureDependentSD;
typedef PorousSolid< ElasticTransverseIsotropic, StrainDependentPermeability > PorousElasticTransverseIsotropicSD;
typedef PorousSolid< ElasticOrthotropic, StrainDependentPermeability > PorousElasticOrthotropicSD;
typedef PorousSolid< DelftEgg, StrainDependentPermeability > PorousDelftEggSD;
typedef PorousSolid< DruckerPrager, StrainDependentPermeability > PorousDruckerPragerSD;
typedef PorousSolid< DruckerPragerExtended, StrainDependentPermeability > PorousDruckerPragerExtendedSD;
typedef PorousSolid< DuvautLionsSolid< DruckerPrager >, StrainDependentPermeability > PorousViscoDruckerPragerSD;
typedef PorousSolid< DuvautLionsSolid< DruckerPragerExtended >, StrainDependentPermeability > PorousViscoDruckerPragerExtendedSD;
typedef PorousSolid< DuvautLionsSolid< ModifiedCamClay >, StrainDependentPermeability > PorousViscoModifiedCamClaySD;
typedef PorousSolid< ModifiedCamClay, StrainDependentPermeability > PorousModifiedCamClaySD;


REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousElasticIsotropicConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousElasticTransverseIsotropicConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousElasticOrthotropicConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousDelftEggConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousDruckerPragerConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousDruckerPragerExtendedConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousModifiedCamClayConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousViscoDruckerPragerConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousViscoDruckerPragerExtendedConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousViscoModifiedCamClayConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousElasticIsotropicCK, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousElasticTransverseIsotropicCK, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousElasticOrthotropicCK, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousDelftEggCK, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousDruckerPragerCK, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousDruckerPragerExtendedCK, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousModifiedCamClayCK, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousViscoDruckerPragerCK, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousViscoDruckerPragerExtendedCK, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousViscoModifiedCamClayCK, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousElasticIsotropicSD, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousElasticIsotropicPressureDependentSD, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousElasticTransverseIsotropicSD, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousElasticOrthotropicSD, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousDelftEggSD, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousDruckerPragerSD, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousDruckerPragerExtendedSD, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousModifiedCamClaySD, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousViscoDruckerPragerSD, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousViscoDruckerPragerExtendedSD, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousViscoModifiedCamClaySD, string const &, Group * const )


}
} /* namespace geos */
