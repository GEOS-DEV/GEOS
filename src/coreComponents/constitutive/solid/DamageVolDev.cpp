
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
 * @file Damage.cpp
 */

#include "Damage.hpp"
#include "DamageVolDev.hpp"

#include "LinearElasticAnisotropic.hpp"
#include "LinearElasticIsotropic.hpp"
#include "LinearElasticTransverseIsotropic.hpp"

namespace geosx
{

using namespace dataRepository;
namespace constitutive
{

template< typename BASE >
DamageVolDev< BASE >::DamageVolDev( string const & name, Group * const parent ):
  Damage<BASE>( name, parent )
{}

template< typename BASE >
DamageVolDev< BASE >::~DamageVolDev()
{}

typedef DamageVolDev< LinearElasticIsotropic > DamageVolDevLinearElasticIsotropic;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DamageVolDevLinearElasticIsotropic, string const &, Group * const )

}
} /* namespace geosx */
