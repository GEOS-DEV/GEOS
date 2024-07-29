/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file ThermoPoroMechanics.cpp
 */

#include "ThermoPoroMechanics.hpp"
#include "constitutive/solid/ElasticIsotropic.hpp"
#include "constitutive/solid/ElasticTransverseIsotropic.hpp"
#include "constitutive/solid/ElasticOrthotropic.hpp"
#include "constitutive/solid/DelftEgg.hpp"
#include "constitutive/solid/DruckerPrager.hpp"
#include "constitutive/solid/DruckerPragerExtended.hpp"
#include "constitutive/solid/Damage.hpp"
#include "constitutive/solid/DamageSpectral.hpp"
#include "constitutive/solid/DamageVolDev.hpp"
#include "constitutive/solid/ModifiedCamClay.hpp"
#include "constitutive/solid/DuvautLionsSolid.hpp"
#include "constitutive/thermalConductivity/SinglePhaseThermalConductivity.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

template< typename SOLID_TYPE >
ThermoPoroMechanics< SOLID_TYPE >::ThermoPoroMechanics( string const & name, Group * const parent ):
  CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >( name, parent )
{}

template< typename SOLID_TYPE >
ThermoPoroMechanics< SOLID_TYPE >::~ThermoPoroMechanics() = default;

template< typename SOLID_TYPE >
void ThermoPoroMechanics< SOLID_TYPE >::initializeState() const
{
  CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >::initializeState();
}

// Register all ThermoPoroMechanics model types.
typedef ThermoPoroMechanics< ElasticIsotropic > ThermoPoroElasticIsotropic;
typedef ThermoPoroMechanics< ElasticTransverseIsotropic > ThermoPoroElasticTransverseIsotropic;
typedef ThermoPoroMechanics< ElasticOrthotropic > ThermoPoroElasticOrthotropic;
typedef ThermoPoroMechanics< DelftEgg > ThermoPoroDelftEgg;
typedef ThermoPoroMechanics< DruckerPrager > ThermoPoroDruckerPrager;
typedef ThermoPoroMechanics< DruckerPragerExtended > ThermoPoroDruckerPragerExtended;
typedef ThermoPoroMechanics< Damage< ElasticIsotropic > > ThermoPoroDamageElasticIsotropic;
typedef ThermoPoroMechanics< DamageSpectral< ElasticIsotropic > > ThermoPoroDamageSpectralElasticIsotropic;
typedef ThermoPoroMechanics< DamageVolDev< ElasticIsotropic > > ThermoPoroDamageVolDevElasticIsotropic;
typedef ThermoPoroMechanics< DuvautLionsSolid< DruckerPrager > > ThermoPoroViscoDruckerPrager;
typedef ThermoPoroMechanics< DuvautLionsSolid< DruckerPragerExtended > > ThermoPoroViscoDruckerPragerExtended;
typedef ThermoPoroMechanics< DuvautLionsSolid< ModifiedCamClay > > ThermoPoroViscoModifiedCamClay;
typedef ThermoPoroMechanics< ModifiedCamClay > ThermoPoroModifiedCamClay;


REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermoPoroElasticIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermoPoroElasticTransverseIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermoPoroElasticOrthotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermoPoroDelftEgg, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermoPoroDruckerPrager, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermoPoroDruckerPragerExtended, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermoPoroDamageElasticIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermoPoroDamageSpectralElasticIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermoPoroDamageVolDevElasticIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermoPoroModifiedCamClay, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermoPoroViscoDruckerPrager, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermoPoroViscoDruckerPragerExtended, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermoPoroViscoModifiedCamClay, string const &, Group * const )


}
} /* namespace geos */
