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
 *  @file ElasticIsotropicPressureDependent.cpp
 */

#include "ElasticIsotropicPressureDependent.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{

ElasticIsotropicPressureDependent::ElasticIsotropicPressureDependent( string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_defaultRefPressure(),
  m_defaultRefStrainVol(),
  m_defaultRecompressionIndex(),
  m_defaultShearModulus(),
  m_refPressure(),
  m_refStrainVol(),
  m_recompressionIndex(),
  m_shearModulus()
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

  registerWrapper( viewKeyStruct::recompressionIndexString(), &m_recompressionIndex ).
    setApplyDefaultValue( -1 ).
    setDescription( "Recompression Index Field" );

  registerWrapper( viewKeyStruct::shearModulusString(), &m_shearModulus ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Shear Modulus" );
}


ElasticIsotropicPressureDependent::~ElasticIsotropicPressureDependent()
{}


void ElasticIsotropicPressureDependent::postProcessInput()
{
  // check what constants the user actually input, and do conversions as needed

  SolidBase::postProcessInput();

  real64 & G  = m_defaultShearModulus;
  real64 & Cr  = m_defaultRecompressionIndex;

  string errorCheck( "( " );
  int numConstantsSpecified = 0;

  if( Cr >= 0.0 )
  {
    ++numConstantsSpecified;
    errorCheck += "Cr, ";
  }
  if( G >= 0.0 )
  {
    ++numConstantsSpecified;
    errorCheck += "G, ";
  }
  errorCheck += ")";

  GEOSX_ERROR_IF( numConstantsSpecified != 2,
                  "A specific pair of elastic constants is required: (Cr, G). " );
  GEOSX_THROW_IF( m_defaultRecompressionIndex <= 0, "Non-positive recompression index detected " << m_defaultRecompressionIndex, InputError );
  real64 poisson = conversions::BulkModAndShearMod::toPoissonRatio( -1*m_defaultRefPressure/m_defaultRecompressionIndex, m_defaultShearModulus );
  GEOSX_THROW_IF( poisson < 0,
                  "Elastic parameters lead to negative Poisson ratio at reference pressure ", InputError );


  // set results as array default values
  this->getWrapper< real64 >( viewKeyStruct::refPressureString() ).
    setApplyDefaultValue( m_defaultRefPressure );

  this->getWrapper< real64 >( viewKeyStruct::refStrainVolString() ).
    setApplyDefaultValue( m_defaultRefStrainVol );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::recompressionIndexString() ).
    setApplyDefaultValue( m_defaultRecompressionIndex );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::shearModulusString() ).
    setApplyDefaultValue( m_defaultShearModulus );

}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ElasticIsotropicPressureDependent, string const &, Group * const )
}
} /* namespace geosx */
