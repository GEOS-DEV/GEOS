/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DamageSpectral.cpp
 */

#include "Damage.hpp"
#include "DamageSpectral.hpp"

#include "ElasticIsotropic.hpp"

namespace geos
{

using namespace dataRepository;
namespace constitutive
{

template< typename BASE >
DamageSpectral< BASE >::DamageSpectral( string const & name, Group * const parent ):
  Damage< BASE >( name, parent )
{
  this->registerWrapper( viewKeyStruct::dissipationFuncitonTypeString(), &m_dissipationFuctionType ).
    setApplyDefaultValue( damageSpectral::LocalDissipation::Linear ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Dissipation function type: "
                    "``" + EnumStrings< damageSpectral::LocalDissipation >::concat( "|" ) + "``" );
}

template< typename BASE >
DamageSpectral< BASE >::~DamageSpectral()
{}

typedef DamageSpectral< ElasticIsotropic > DamageSpectralElasticIsotropic;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DamageSpectralElasticIsotropic, string const &, Group * const )

}
} /* namespace geos */
