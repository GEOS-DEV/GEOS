/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
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
#include "ElasticTransverseIsotropic.hpp"
#include "DruckerPrager.hpp"
#include "DruckerPragerExtended.hpp"
#include "Damage.hpp"
#include "DamageSpectral.hpp"
#include "DamageVolDev.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

template< typename SOLID_TYPE >
PorousSolid< SOLID_TYPE >::PorousSolid( string const & name, Group * const parent ):
  CoupledSolid< SOLID_TYPE, BiotPorosity, StrainDependentPermeability >( name, parent )
{}

template< typename SOLID_TYPE >
PorousSolid< SOLID_TYPE >::~PorousSolid() = default;

// Register all PorousSolid model types.
typedef PorousSolid< ElasticIsotropic > PoroElasticIsotropic;
typedef PorousSolid< ElasticTransverseIsotropic > PoroElasticTransverseIsotropic;
typedef PorousSolid< DruckerPrager > PoroDruckerPrager;
typedef PorousSolid< DruckerPragerExtended > PoroDruckerPragerExtended;
typedef PorousSolid< Damage< ElasticIsotropic > > PoroDamageElasticIsotropic;
typedef PorousSolid< DamageSpectral< ElasticIsotropic > > PoroDamageSpectralElasticIsotropic;
typedef PorousSolid< DamageVolDev< ElasticIsotropic > > PoroDamageVolDevElasticIsotropic;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoroElasticIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoroElasticTransverseIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoroDruckerPrager, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoroDruckerPragerExtended, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoroDamageElasticIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoroDamageSpectralElasticIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoroDamageVolDevElasticIsotropic, string const &, Group * const )

}
} /* namespace geosx */
