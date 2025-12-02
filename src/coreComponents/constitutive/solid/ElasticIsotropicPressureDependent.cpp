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
 *  @file ElasticIsotropicPressureDependent.cpp
 */

#include "ElasticIsotropicPressureDependent.hpp"
#include "SolidFields.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

ElasticIsotropicPressureDependent::ElasticIsotropicPressureDependent( string const & name, Group * const parent ):
  SolidBase( name, parent )
{
  registerWrapper( viewKeyStruct::defaultRefPressureString(), &m_defaultRefPressure ).
    setApplyDefaultValue( -1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference Pressure" );

  registerWrapper( viewKeyStruct::defaultRefStrainVolString(), &m_defaultRefStrainVol ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference Volumetric Strain" );

  registerWrapper( viewKeyStruct::defaultRecompressionIndexString(), &m_defaultRecompressionIndex ).
    setApplyDefaultValue( 2e-3 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Recompresion Index" );

  registerWrapper( viewKeyStruct::defaultShearModulusString(), &m_defaultShearModulus ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Elastic Shear Modulus Parameter" );

  registerWrapper( viewKeyStruct::refPressureString(), &m_refPressure ).
    setApplyDefaultValue( -1 ).
    setDescription( "Reference Pressure Field" );

  registerWrapper( viewKeyStruct::refStrainVolString(), &m_refStrainVol ).
    setApplyDefaultValue( -1 ).
    setDescription( "Reference Volumetric Strain" );

  // register fields

  registerField< fields::solid::recompressionIndex >( &m_recompressionIndex );

  registerField< fields::solid::shearModulus >( &m_shearModulus );
}


void ElasticIsotropicPressureDependent::postInputInitialization()
{
  // check what constants the user actually input, and do conversions as needed

  SolidBase::postInputInitialization();

  real64 & G  = m_defaultShearModulus;
  real64 & Cr  = m_defaultRecompressionIndex;

  string errorCheck( "( " );
  int numConstantsSpecified = 0;

  if( Cr >= 0.0 )
  {
    ++numConstantsSpecified;
    errorCheck += "Cr ";
  }
  if( G >= 0.0 )
  {
    ++numConstantsSpecified;
    errorCheck += "G ";
  }
  errorCheck += ")";

  GEOS_ERROR_IF( numConstantsSpecified != 2,
                 GEOS_FMT( "{}: A specific pair of elastic constants is required: ( Cr G ), specified: {}",
                           getFullName(), errorCheck ) );
  GEOS_THROW_IF( m_defaultRefPressure >= 0,
                 GEOS_FMT( "{}: Reference pressure must be negative", getFullName() ),
                 InputError );
  GEOS_THROW_IF( m_defaultRecompressionIndex <= 0,
                 GEOS_FMT( "{}: Non-positive recompression index detected {}", getFullName(), m_defaultRecompressionIndex ),
                 InputError, getDataContext() );
  real64 poisson =
    conversions::bulkModAndShearMod::toPoissonRatio( -1 * m_defaultRefPressure / m_defaultRecompressionIndex, m_defaultShearModulus );
  GEOS_THROW_IF( poisson < 0,
                 GEOS_FMT( "{}: Elastic parameters lead to negative Poisson ratio at reference pressure", getFullName() ),
                 InputError, getDataContext() );

  // set results as array default values

  getWrapper< real64 >( viewKeyStruct::refPressureString() ).
    setApplyDefaultValue( m_defaultRefPressure );

  getWrapper< real64 >( viewKeyStruct::refStrainVolString() ).
    setApplyDefaultValue( m_defaultRefStrainVol );

  getField< fields::solid::recompressionIndex >().
    setApplyDefaultValue( m_defaultRecompressionIndex );

  getField< fields::solid::shearModulus >().
    setApplyDefaultValue( m_defaultShearModulus );

}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ElasticIsotropicPressureDependent, string const &, Group * const )
}
} /* namespace geos */
