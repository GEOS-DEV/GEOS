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
 *  @file LinearElasticTransverseIsotropic.cpp
 */

#include "LinearElasticTransverseIsotropic.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{

LinearElasticTransverseIsotropic::LinearElasticTransverseIsotropic( string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_defaultYoungsModulusTransverse(),
  m_defaultYoungsModulusAxial(),
  m_defaultPoissonsRatioTransverse(),
  m_defaultPoissonsRatioAxialTransverse(),
  m_defaultShearModulusAxialTransverse(),
  m_c11(),
  m_c13(),
  m_c33(),
  m_c44(),
  m_c66()
{
  registerWrapper( viewKeyStruct::defaultYoungsModulusTransverse, &m_defaultYoungsModulusTransverse )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Elastic Transverse Young's Modulus Parameter" );

  registerWrapper( viewKeyStruct::defaultYoungsModulusAxial, &m_defaultYoungsModulusAxial )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Elastic Axial Young's Modulus Parameter" );

  registerWrapper( viewKeyStruct::defaultPoissonsRatioTransverse, &m_defaultPoissonsRatioTransverse )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Elastic Transverse Poisson's ratio Parameter" );

  registerWrapper( viewKeyStruct::defaultPoissonsRatioAxialTransverse, &m_defaultPoissonsRatioAxialTransverse )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Elastic Axial-Transverse Poisson's ratio Parameter" );

  registerWrapper( viewKeyStruct::defaultShearModulusAxialTransverse, &m_defaultShearModulusAxialTransverse )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Elastic Axial-Transverse (out-plane) Shear Modulus Parameter" );

  registerWrapper< real64 >( viewKeyStruct::defaultC11 )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Elastic Stiffness Parameter C11" );

  registerWrapper< real64 >( viewKeyStruct::defaultC13 )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Elastic Stiffness Parameter C13" );

  registerWrapper< real64 >( viewKeyStruct::defaultC33 )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Elastic Stiffness Parameter C33" );

  registerWrapper< real64 >( viewKeyStruct::defaultC44 )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Elastic Stiffness Parameter C44" );

  registerWrapper< real64 >( viewKeyStruct::defaultC66 )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Elastic Stiffness Parameter C66" );

  registerWrapper( viewKeyStruct::c11, &m_c11 )->
    setApplyDefaultValue( -1 )->
    setDescription( "Elastic Stiffness Field C11" );

  registerWrapper( viewKeyStruct::c13, &m_c13 )->
    setApplyDefaultValue( -1 )->
    setDescription( "Elastic Stiffness Field C13" );

  registerWrapper( viewKeyStruct::c33, &m_c33 )->
    setApplyDefaultValue( -1 )->
    setDescription( "Elastic Stiffness Field C33" );

  registerWrapper( viewKeyStruct::c44, &m_c44 )->
    setApplyDefaultValue( -1 )->
    setDescription( "Elastic Stiffness Field C44" );

  registerWrapper( viewKeyStruct::c66, &m_c66 )->
    setApplyDefaultValue( -1 )->
    setDescription( "Elastic Stiffness Field C66" );

}

LinearElasticTransverseIsotropic::~LinearElasticTransverseIsotropic()
{}

void LinearElasticTransverseIsotropic::postProcessInput()
{

  SolidBase::postProcessInput();

  real64 & c11  = getReference< real64 >( viewKeyStruct::defaultC11 );
  real64 & c13  = getReference< real64 >( viewKeyStruct::defaultC13 );
  real64 & c33  = getReference< real64 >( viewKeyStruct::defaultC33 );
  real64 & c44  = getReference< real64 >( viewKeyStruct::defaultC44 );
  real64 & c66  = getReference< real64 >( viewKeyStruct::defaultC66 );

  real64 & Et = m_defaultYoungsModulusTransverse;
  real64 & Ea = m_defaultYoungsModulusAxial;
  real64 & Nut = m_defaultPoissonsRatioTransverse;
  real64 & Nuat = m_defaultPoissonsRatioAxialTransverse;
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

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c11 )->
    setApplyDefaultValue( c11 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c13 )->
    setApplyDefaultValue( c13 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c33 )->
    setApplyDefaultValue( c33 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c44 )->
    setApplyDefaultValue( c44 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c66 )->
    setApplyDefaultValue( c66 );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, LinearElasticTransverseIsotropic, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geosx */
