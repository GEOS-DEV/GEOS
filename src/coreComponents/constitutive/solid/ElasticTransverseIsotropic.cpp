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
 *  @file ElasticTransverseIsotropic.cpp
 */

#include "ElasticTransverseIsotropic.hpp"
#include "SolidFields.hpp"

namespace geos
{
using namespace dataRepository;

namespace constitutive
{

ElasticTransverseIsotropic::ElasticTransverseIsotropic( string const & name, Group * const parent ):
  SolidBase( name, parent )
{
  registerWrapper( viewKeyStruct::defaultYoungModulusTransverseString(), &m_defaultYoungModulusTransverse ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Transverse Young's Modulus" );

  registerWrapper( viewKeyStruct::defaultYoungModulusAxialString(), &m_defaultYoungModulusAxial ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Axial Young's Modulus" );

  registerWrapper( viewKeyStruct::defaultPoissonRatioTransverseString(), &m_defaultPoissonRatioTransverse ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Transverse Poisson's Ratio" );

  registerWrapper( viewKeyStruct::defaultPoissonRatioAxialTransverseString(), &m_defaultPoissonRatioAxialTransverse ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Axial-Transverse Poisson's Ratio" );

  registerWrapper( viewKeyStruct::defaultShearModulusAxialTransverseString(), &m_defaultShearModulusAxialTransverse ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Axial-Transverse Shear Modulus" );

  registerWrapper< real64 >( viewKeyStruct::defaultC11String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Stiffness Parameter C11" );

  registerWrapper< real64 >( viewKeyStruct::defaultC13String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Stiffness Parameter C13" );

  registerWrapper< real64 >( viewKeyStruct::defaultC33String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Stiffness Parameter C33" );

  registerWrapper< real64 >( viewKeyStruct::defaultC44String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Stiffness Parameter C44" );

  registerWrapper< real64 >( viewKeyStruct::defaultC66String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Stiffness Parameter C66" );

  // register fields - TODO merge in one

  registerField< fields::solid::c11 >( &m_c11 );
  registerField< fields::solid::c13 >( &m_c13 );
  registerField< fields::solid::c33 >( &m_c33 );
  registerField< fields::solid::c44 >( &m_c44 );
  registerField< fields::solid::c66 >( &m_c66 );
}

void ElasticTransverseIsotropic::postInputInitialization()
{
  SolidBase::postInputInitialization();

  real64 & c11  = getReference< real64 >( viewKeyStruct::defaultC11String() );
  real64 & c13  = getReference< real64 >( viewKeyStruct::defaultC13String() );
  real64 & c33  = getReference< real64 >( viewKeyStruct::defaultC33String() );
  real64 & c44  = getReference< real64 >( viewKeyStruct::defaultC44String() );
  real64 & c66  = getReference< real64 >( viewKeyStruct::defaultC66String() );

  real64 & Et = m_defaultYoungModulusTransverse;
  real64 & Ea = m_defaultYoungModulusAxial;
  real64 & Nut = m_defaultPoissonRatioTransverse;
  real64 & Nuat = m_defaultPoissonRatioAxialTransverse;
  real64 & Gat = m_defaultShearModulusAxialTransverse;

  if( Et > 0.0 && Ea > 0.0 && Gat > 0.0 && Nut > -0.5 && Nut < 0.5 )
  {
    real64 const Nuta = Nuat * ( Et / Ea );
    real64 const delta = ( 1.0 + Nut ) * ( 1.0 - Nut - 2.0 * Nuta * Nuat );

    if( delta > 0.0 && Nuta * Nuat < 1.0 )
    {
      c11 = ( 1.0 - Nuta * Nuat ) * Et / delta;
      c13 = Nuat * ( 1.0 + Nut ) * Et / delta;
      c33 = ( 1.0 - Nut * Nut ) * Ea / delta;
      c44 = Gat;
      c66 = 0.5 * Et / ( 1.0 + Nut );
    }
  }

  getField< fields::solid::c11 >().
    setApplyDefaultValue( c11 );
  getField< fields::solid::c13 >().
    setApplyDefaultValue( c13 );
  getField< fields::solid::c33 >().
    setApplyDefaultValue( c33 );
  getField< fields::solid::c44 >().
    setApplyDefaultValue( c44 );
  getField< fields::solid::c66 >().
    setApplyDefaultValue( c66 );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ElasticTransverseIsotropic, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
