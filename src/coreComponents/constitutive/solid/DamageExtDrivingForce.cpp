
/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DamageExtDrivingForce.cpp
 */

#include "Damage.hpp"
#include "DamageExtDrivingForce.hpp"

#include "ElasticIsotropic.hpp"

namespace geosx
{

using namespace dataRepository;
namespace constitutive
{

template< typename BASE >
DamageExtDrivingForce< BASE >::DamageExtDrivingForce( string const & name, Group * const parent ):
  Damage< BASE >( name, parent ),
  m_tensileStrength(), 
  m_compressStrength(), 
{
  this->registerWrapper( viewKeyStruct::tensileStrengthString(), &m_tensileStrength ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Tensile strength from the uniaxial tension test" );

  this->registerWrapper( viewKeyStruct::compressStrengthString(), &m_compressStrength ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Compressive strength from the uniaxial compression test" );

  this->registerWrapper( viewKeyStruct::deltaCoefficientString(), &m_deltaCoefficient ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Coefficient in the calculation of the external driving force" );

}

typedef DamageExtDrivingForce< ElasticIsotropic > DamageExtDrivingForceElasticIsotropic;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DamageExtDrivingForceElasticIsotropic, string const &, Group * const )

}
} /* namespace geosx */
