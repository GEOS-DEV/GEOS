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

  registerField< fields::solid::youngModulus >( &m_youngModulus );

  registerField< fields::solid::poissonRatio >( &m_poissonRatio );
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
                 GEOS_FMT( "A specific pair of elastic constants is required. "
                           "Either (K,G), (K,E), (G,E), (K,nu), (G,nu) or (E,nu). "
                           "You have specified {}",
                           errorCheck ),
                 getDataContext() );

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
    GEOS_ERROR( GEOS_FMT( "Invalid specification for default elastic constants. {} has been specified.",
                          errorCheck ),
                getDataContext() );
  }

  // set results as array default values

  getField< fields::solid::bulkModulus >().
    setApplyDefaultValue( m_defaultBulkModulus );

  getField< fields::solid::shearModulus >().
    setApplyDefaultValue( m_defaultShearModulus );
}

void ElasticIsotropic::initializePostInitialConditionsPreSubGroups()
{
  SolidBase::initializePostInitialConditionsPreSubGroups();

  // If per-cell Young's modulus and Poisson's ratio were imported from an external mesh
  // (indicated by youngModulus[k] > 0 and poissonRatio[k] in (-0.5, 0.5)),
  // convert them to bulk and shear modulus on a cell-by-cell basis.
  arrayView1d< real64 const > const youngMod = m_youngModulus;
  arrayView1d< real64 const > const nu       = m_poissonRatio;
  arrayView1d< real64 > const bulkMod        = m_bulkModulus;
  arrayView1d< real64 > const shearMod       = m_shearModulus;

  localIndex numConverted = 0;
  localIndex numInvalidE  = 0;

  for( localIndex k = 0; k < m_youngModulus.size(); ++k )
  {
    // youngModulus default is 0: negative values are invalid, zero means not imported
    if( youngMod[k] < 0.0 )
    {
      GEOS_WARNING( GEOS_FMT( "ElasticIsotropic '{}': element {} has negative Young's modulus ({:.6e}). "
                              "Skipping per-cell conversion; default bulk/shear modulus will be used.",
                              getName(), k, youngMod[k] ) );
      ++numInvalidE;
      continue;
    }

    // youngModulus default is 0: not imported — skip silently
    if( !( youngMod[k] > 0.0 ) )
      continue;

    // E was imported and is positive: nu must also be valid
    GEOS_ERROR_IF( nu[k] <= -0.5 || nu[k] >= 0.5,
                   GEOS_FMT( "ElasticIsotropic '{}': element {} has Young's modulus imported ({:.6e}) "
                             "but Poisson's ratio ({:.6f}) is outside the valid range (-0.5, 0.5). "
                             "Both youngModulus and poissonRatio must be provided together.",
                             getName(), k, youngMod[k], nu[k] ) );

    bulkMod[k]  = conversions::youngModAndPoissonRatio::toBulkMod( youngMod[k], nu[k] );
    shearMod[k] = conversions::youngModAndPoissonRatio::toShearMod( youngMod[k], nu[k] );
    ++numConverted;
  }

  GEOS_LOG_RANK_0_IF( numConverted > 0,
                      GEOS_FMT( "ElasticIsotropic '{}': converted per-cell Young's modulus / Poisson's ratio "
                                "to bulk / shear modulus for {} element(s).",
                                getName(), numConverted ) );
  GEOS_WARNING_IF( numInvalidE > 0,
                   GEOS_FMT( "ElasticIsotropic '{}': {} element(s) had non-positive Young's modulus and were skipped.",
                             getName(), numInvalidE ) );

  // Back-compute E and nu for all cells from the final K/G so that output fields are meaningful
  // for both imported and default cells.
  for( localIndex k = 0; k < m_bulkModulus.size(); ++k )
  {
    m_youngModulus[k] = conversions::bulkModAndShearMod::toYoungMod( bulkMod[k], shearMod[k] );
    m_poissonRatio[k] = conversions::bulkModAndShearMod::toPoissonRatio( bulkMod[k], shearMod[k] );
  }
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ElasticIsotropic, string const &, Group * const )
}
} /* namespace geos */
