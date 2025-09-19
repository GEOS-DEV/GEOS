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
 *  @file CeramicDamage.cpp
 */

#include "CeramicDamage.hpp"
#include "SolidFields.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

CeramicDamage::CeramicDamage( string const & name, Group * const parent ):
  ElasticIsotropic( name, parent )
{
  // register default values
  registerWrapper( viewKeyStruct::tensileStrengthString(), &m_tensileStrength ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Tensile strength" );

  registerWrapper( viewKeyStruct::compressiveStrengthString(), &m_compressiveStrength ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Compressive strength" );

  registerWrapper( viewKeyStruct::maximumStrengthString(), &m_maximumStrength ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum theoretical strength" );

  registerWrapper( viewKeyStruct::crackSpeedString(), &m_crackSpeed ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Crack speed" );

  // register fields
  registerField< fields::solid::damage >( &m_damage );

  registerField< fields::solid::jacobian >( &m_jacobian );

  registerField< fields::solid::lengthScale >( &m_lengthScale );
}


void CeramicDamage::allocateConstitutiveData( Group & parent, localIndex const numPts )
{
  m_damage.resize( 0, numPts );
  m_jacobian.resize( 0, numPts );

  ElasticIsotropic::allocateConstitutiveData( parent, numPts );
}


void CeramicDamage::postInputInitialization()
{
  ElasticIsotropic::postInputInitialization();

  GEOS_THROW_IF( m_tensileStrength < 0.0, "Tensile strength must be a positive number.", InputError );
  GEOS_THROW_IF( m_compressiveStrength < m_tensileStrength, "Compressive strength must be greater than tensile strength.", InputError );
  GEOS_THROW_IF( m_maximumStrength < m_compressiveStrength, "Maximum theoretical strength must be greater than compressive strength.", InputError );
  GEOS_THROW_IF( m_crackSpeed < 0.0, "Crack speed must be a positive number.", InputError );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, CeramicDamage, std::string const &, Group * const )
}
} /* namespace geos */
