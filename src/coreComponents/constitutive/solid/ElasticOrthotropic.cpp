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
 *  @file ElasticOrthotropic.cpp
 */

#include "ElasticOrthotropic.hpp"
#include "SolidFields.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

ElasticOrthotropic::ElasticOrthotropic( string const & name, Group * const parent ):
  SolidBase( name, parent )
{
  registerWrapper( viewKeyStruct::defaultE1String(), &m_defaultE1 ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Young's Modulus E1" );

  registerWrapper( viewKeyStruct::defaultE2String(), &m_defaultE2 ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Young's Modulus E2" );

  registerWrapper( viewKeyStruct::defaultE3String(), &m_defaultE3 ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Young's Modulus E3" );

  registerWrapper( viewKeyStruct::defaultNu12String(), &m_defaultNu12 ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Poission's Ratio Nu12" );

  registerWrapper( viewKeyStruct::defaultNu13String(), &m_defaultNu13 ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Poission's Ratio Nu13" );

  registerWrapper( viewKeyStruct::defaultNu23String(), &m_defaultNu23 ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Poission's Ratio Nu23" );

  registerWrapper( viewKeyStruct::defaultG12String(), &m_defaultG12 ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Shear Modulus G12" );

  registerWrapper( viewKeyStruct::defaultG13String(), &m_defaultG13 ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Shear Modulus G13" );

  registerWrapper( viewKeyStruct::defaultG23String(), &m_defaultG23 ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Shear Modulus G23" );

  registerWrapper< real64 >( viewKeyStruct::defaultC11String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default C11 Component of Voigt Stiffness Tensor" );

  registerWrapper< real64 >( viewKeyStruct::defaultC12String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default C12 Component of Voigt Stiffness Tensor" );

  registerWrapper< real64 >( viewKeyStruct::defaultC13String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default C13 Component of Voigt Stiffness Tensor" );

  registerWrapper< real64 >( viewKeyStruct::defaultC22String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default C22 Component of Voigt Stiffness Tensor" );

  registerWrapper< real64 >( viewKeyStruct::defaultC23String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default C23 Component of Voigt Stiffness Tensor" );

  registerWrapper< real64 >( viewKeyStruct::defaultC33String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default C33 Component of Voigt Stiffness Tensor" );

  registerWrapper< real64 >( viewKeyStruct::defaultC44String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default C44 Component of Voigt Stiffness Tensor" );

  registerWrapper< real64 >( viewKeyStruct::defaultC55String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default C55 Component of Voigt Stiffness Tensor" );

  registerWrapper< real64 >( viewKeyStruct::defaultC66String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default C66 Component of Voigt Stiffness Tensor" );

  // register fields - TODO merge in one

  registerField< fields::solid::c11 >( &m_c11 );
  registerField< fields::solid::c12 >( &m_c12 );
  registerField< fields::solid::c13 >( &m_c13 );
  registerField< fields::solid::c22 >( &m_c22 );
  registerField< fields::solid::c23 >( &m_c23 );
  registerField< fields::solid::c33 >( &m_c33 );
  registerField< fields::solid::c44 >( &m_c44 );
  registerField< fields::solid::c55 >( &m_c55 );
  registerField< fields::solid::c66 >( &m_c66 );
}

void ElasticOrthotropic::postInputInitialization()
{
  SolidBase::postInputInitialization();

  real64 & c11 = getReference< real64 >( viewKeyStruct::defaultC11String() );
  real64 & c12 = getReference< real64 >( viewKeyStruct::defaultC12String() );
  real64 & c13 = getReference< real64 >( viewKeyStruct::defaultC13String() );
  real64 & c22 = getReference< real64 >( viewKeyStruct::defaultC22String() );
  real64 & c23 = getReference< real64 >( viewKeyStruct::defaultC23String() );
  real64 & c33 = getReference< real64 >( viewKeyStruct::defaultC33String() );
  real64 & c44 = getReference< real64 >( viewKeyStruct::defaultC44String() );
  real64 & c55 = getReference< real64 >( viewKeyStruct::defaultC55String() );
  real64 & c66 = getReference< real64 >( viewKeyStruct::defaultC66String() );

  real64 const E1   = m_defaultE1;
  real64 const E2   = m_defaultE2;
  real64 const E3   = m_defaultE3;
  real64 const Nu12 = m_defaultNu12;
  real64 const Nu13 = m_defaultNu13;
  real64 const Nu23 = m_defaultNu23;
  real64 const G12  = m_defaultG12;
  real64 const G13  = m_defaultG13;
  real64 const G23  = m_defaultG23;

  if( E1 > 0.0 && E2 > 0.0 && E3 > 0.0 && G12 > 0.0 && G13 > 0.0 && G23 > 0.0 )
  {
    real64 const Nu21 = Nu12 * E2 / E1;
    real64 const Nu31 = Nu13 * E3 / E1;
    real64 const Nu32 = Nu23 * E3 / E2;

    real64 const delta = 1.0 - Nu12 * Nu21 - Nu13 * Nu31 - Nu23 * Nu32 - 2.0 * Nu12 * Nu23 * Nu31;

    if( delta > 0.0 && Nu23 * Nu32 < 1.0 && Nu31 * Nu13 < 1.0 && Nu12 * Nu21 < 1.0 )
    {
      c11 = ( 1.0 - Nu23 * Nu32 )  * E1 / delta;
      c12 = ( Nu21 + Nu31 * Nu23 ) * E1 / delta;
      c13 = ( Nu31 + Nu21 * Nu32 ) * E1 / delta;

      c22 = ( 1.0 - Nu31 * Nu13 )  * E2 / delta;
      c23 = ( Nu32 + Nu31 * Nu12 ) * E2 / delta;

      c33 = ( 1.0 - Nu12 * Nu21 )  * E3 / delta;

      c44 = G23;
      c55 = G13;
      c66 = G12;
    }
    else
    {
      GEOS_ERROR( "Invalid specification for default elastic constants.", getDataContext());
    }
  }
  else if( c11 <= 0.0 || c22 <= 0.0 || c33 <= 0.0 || c44 <= 0.0 || c55 <= 0.0 || c66 <= 0.0 )
  {
    GEOS_ERROR( "Invalid specification for default elastic stiffnesses.", getDataContext());
  }

  // TODO merge in one
  getField< fields::solid::c11 >().
    setApplyDefaultValue( c11 );
  getField< fields::solid::c12 >().
    setApplyDefaultValue( c12 );
  getField< fields::solid::c13 >().
    setApplyDefaultValue( c13 );
  getField< fields::solid::c22 >().
    setApplyDefaultValue( c22 );
  getField< fields::solid::c23 >().
    setApplyDefaultValue( c23 );
  getField< fields::solid::c33 >().
    setApplyDefaultValue( c33 );
  getField< fields::solid::c44 >().
    setApplyDefaultValue( c44 );
  getField< fields::solid::c55 >().
    setApplyDefaultValue( c55 );
  getField< fields::solid::c66 >().
    setApplyDefaultValue( c66 );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ElasticOrthotropic, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
