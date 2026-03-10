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
 *  @file ElasticIsotropic.cpp
 */

#include "ElasticIsotropic.hpp"
#include "SolidFields.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

ElasticIsotropic::ElasticIsotropic( string const & name, Group * const parent ):
  SolidBase( name, parent )
{
  registerWrapper( viewKeyStruct::defaultBulkModulusString(), &m_defaultBulkModulus ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Bulk Modulus Parameter" );

  registerWrapper( viewKeyStruct::defaultShearModulusString(), &m_defaultShearModulus ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Shear Modulus Parameter" );

  registerWrapper< real64 >( viewKeyStruct::defaultYoungModulusString() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Young's Modulus" );

  registerWrapper< real64 >( viewKeyStruct::defaultPoissonRatioString() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Poisson's Ratio" );

  // register fields

  registerField< fields::solid::bulkModulus >( &m_bulkModulus );

  registerField< fields::solid::shearModulus >( &m_shearModulus );
}

void ElasticIsotropic::postInputInitialization()
{
  // check what constants the user actually input, and do conversions as needed

  SolidBase::postInputInitialization();

  real64 & nu = getReference< real64 >( viewKeyStruct::defaultPoissonRatioString() );
  real64 & E  = getReference< real64 >( viewKeyStruct::defaultYoungModulusString() );
  real64 & K  = m_defaultBulkModulus;
  real64 & G  = m_defaultShearModulus;

  // Poisson ratio range is: -0.5 < nu < 0.5
  // Zero bulk, shear or Young modulus is not accepted to avoid devision to zero

  string errorCheck( "( " );
  int numConstantsSpecified = 0;
  if( nu > -0.5 && nu < 0.5 )
  {
    ++numConstantsSpecified;
    errorCheck += "nu, ";
  }
  if( E > 0.0 )
  {
    ++numConstantsSpecified;
    errorCheck += "E, ";
  }
  if( K > 0.0 )
  {
    ++numConstantsSpecified;
    errorCheck += "K, ";
  }
  if( G > 0.0 )
  {
    ++numConstantsSpecified;
    errorCheck += "G, ";
  }
  errorCheck += ")";

  GEOS_ERROR_IF( numConstantsSpecified != 2,
                 getFullName() << ": A specific pair of elastic constants is required. " <<
                 "Either (K,G), (K,E), (G,E), (K,nu), (G,nu) or (E,nu). " <<
                 "You have specified " << errorCheck );

  if( nu > -0.5 && nu < 0.5 && E > 0.0 )
  {
    K = conversions::youngModAndPoissonRatio::toBulkMod( E, nu );
    G = conversions::youngModAndPoissonRatio::toShearMod( E, nu );
  }
  else if( nu > -0.5 && nu < 0.5 && G > 0.0 )
  {
    E = conversions::shearModAndPoissonRatio::toYoungMod( G, nu );
    K = conversions::shearModAndPoissonRatio::toBulkMod( G, nu );
  }
  else if( nu > -0.5 && nu < 0.5 && K > 0.0 )
  {
    E = conversions::bulkModAndPoissonRatio::toYoungMod( K, nu );
    G = conversions::bulkModAndPoissonRatio::toShearMod( K, nu );
  }
  else if( E > 0.0 && K > 0.0 )
  {
    nu = conversions::bulkModAndYoungMod::toPoissonRatio( K, E );
    G  = conversions::bulkModAndYoungMod::toShearMod( K, E );
  }
  else if( E > 0.0 && G > 0.0 )
  {
    nu = conversions::shearModAndYoungMod::toPoissonRatio( G, E );
    K  = conversions::shearModAndYoungMod::toBulkMod( G, E );
  }
  else if( K > 0.0 && G > 0.0 )
  {
    E  = conversions::bulkModAndShearMod::toYoungMod( K, G );
    nu = conversions::bulkModAndShearMod::toPoissonRatio( K, G );
  }
  else
  {
    GEOS_ERROR( getFullName() << ": Invalid specification for default elastic constants. " <<
                errorCheck << " has been specified." );
  }

  // set results as array default values

  getField< fields::solid::bulkModulus >().
    setApplyDefaultValue( m_defaultBulkModulus );

  getField< fields::solid::shearModulus >().
    setApplyDefaultValue( m_defaultShearModulus );

  this->getWrapper< array2d< real64 > >( viewKeyStruct::wavespeedString() ).
    setApplyDefaultValue( sqrt( ( m_defaultBulkModulus + (4.0/3.0) * m_defaultShearModulus ) / m_defaultDensity ) );
}

void ElasticIsotropic::allocateConstitutiveData( dataRepository::Group & parent,
                                                 localIndex const numConstitutivePointsPerParentIndex )
{
  ContinuumBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_bulkModulus.resize( 0 );
  m_shearModulus.resize( 0 );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ElasticIsotropic, string const &, Group * const )
}
} /* namespace geos */
