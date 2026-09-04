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
  // Preserve historical behavior: DamageSpectral defaults to the Cohesive degradation function.
  this->template getWrapper< FractureModelType >( Damage< BASE >::viewKeyStruct::fractureModelTypeString() ).
    setApplyDefaultValue( FractureModelType::Cohesive );
}

template< typename BASE >
void DamageSpectral< BASE >::postInputInitialization()
{
  Damage< BASE >::postInputInitialization();

  GEOS_ERROR_IF( this->getFractureModelType() == FractureModelType::Nucleation,
                 "the Nucleation crack model is not supported with the Spectral split",
                 this->getDataContext() );
}

typedef DamageSpectral< ElasticIsotropic > DamageSpectralElasticIsotropic;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DamageSpectralElasticIsotropic, string const &, Group * const )

}
} /* namespace geos */
