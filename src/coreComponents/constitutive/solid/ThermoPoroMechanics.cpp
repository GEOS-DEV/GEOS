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

template< typename SOLID_TYPE,
          typename THERMAL_COND_TYPE >
ThermoPoroMechanics< SOLID_TYPE, THERMAL_COND_TYPE >::ThermoPoroMechanics( string const & name, Group * const parent ):
  CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >( name, parent )
{}

template< typename SOLID_TYPE,
          typename THERMAL_COND_TYPE >
ThermoPoroMechanics< SOLID_TYPE, THERMAL_COND_TYPE >::~ThermoPoroMechanics() = default;

template< typename SOLID_TYPE,
          typename THERMAL_COND_TYPE >
void ThermoPoroMechanics< SOLID_TYPE, THERMAL_COND_TYPE >::initializeState() const
{
  CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >::initializeState();
}

// Register all ThermoPoroMechanics model types.
typedef ThermoPoroMechanics< ElasticIsotropic, SinglePhaseThermalConductivity > ThermoPoroElasticIsotropic;
typedef ThermoPoroMechanics< ElasticTransverseIsotropic, SinglePhaseThermalConductivity > ThermoPoroElasticTransverseIsotropic;
typedef ThermoPoroMechanics< ElasticOrthotropic, SinglePhaseThermalConductivity > ThermoPoroElasticOrthotropic;
typedef ThermoPoroMechanics< DelftEgg, SinglePhaseThermalConductivity > ThermoPoroDelftEgg;
typedef ThermoPoroMechanics< DruckerPrager, SinglePhaseThermalConductivity > ThermoPoroDruckerPrager;
typedef ThermoPoroMechanics< DruckerPragerExtended, SinglePhaseThermalConductivity > ThermoPoroDruckerPragerExtended;
typedef ThermoPoroMechanics< Damage< ElasticIsotropic >, SinglePhaseThermalConductivity > ThermoPoroDamageElasticIsotropic;
typedef ThermoPoroMechanics< DamageSpectral< ElasticIsotropic >, SinglePhaseThermalConductivity > ThermoPoroDamageSpectralElasticIsotropic;
typedef ThermoPoroMechanics< DamageVolDev< ElasticIsotropic >, SinglePhaseThermalConductivity > ThermoPoroDamageVolDevElasticIsotropic;
typedef ThermoPoroMechanics< DuvautLionsSolid< DruckerPrager >, SinglePhaseThermalConductivity > ThermoPoroViscoDruckerPrager;
typedef ThermoPoroMechanics< DuvautLionsSolid< DruckerPragerExtended >, SinglePhaseThermalConductivity > ThermoPoroViscoDruckerPragerExtended;
typedef ThermoPoroMechanics< DuvautLionsSolid< ModifiedCamClay >, SinglePhaseThermalConductivity > ThermoPoroViscoModifiedCamClay;
typedef ThermoPoroMechanics< ModifiedCamClay, SinglePhaseThermalConductivity > ThermoPoroModifiedCamClay;


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
