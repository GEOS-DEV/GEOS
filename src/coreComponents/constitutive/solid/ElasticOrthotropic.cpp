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
 *  @file ElasticOrthotropic.cpp
 */

#include "ElasticOrthotropic.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{

ElasticOrthotropic::ElasticOrthotropic( string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_defaultYoungsModulusTransverse(),
  m_defaultYoungsModulusAxial(),
  m_defaultPoissonTransverse(),
  m_defaultPoissonAxialTransverse(),
  m_defaultShearModulusAxialTransverse(),
  m_defaultC11(),
  m_defaultC12(),
  m_defaultC13(),
  m_defaultC22(),
  m_defaultC23(),
  m_defaultC33(),
  m_defaultC44(),
  m_defaultC55(),
  m_defaultC66(),
  m_c11(),
  m_c12(),
  m_c13(),
  m_c22(),
  m_c23(),
  m_c33(),
  m_c44(),
  m_c55(),
  m_c66()
{
  registerWrapper( viewKeyStruct::defaultYoungsModulusTransverse, &m_defaultYoungsModulusTransverse )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Elastic Bulk Modulus Parameter" );

  registerWrapper( viewKeyStruct::defaultYoungsModulusAxial, &m_defaultYoungsModulusAxial )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Elastic Shear Modulus Parameter" );

  registerWrapper( viewKeyStruct::defaultPoissonRatioTransverse, &m_defaultPoissonTransverse )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Elastic Shear Modulus Parameter" );

  registerWrapper( viewKeyStruct::defaultPoissonRatioAxialTransverse, &m_defaultPoissonAxialTransverse )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Elastic Shear Modulus Parameter" );

  registerWrapper( viewKeyStruct::defaultShearModulusAxialTransverse, &m_defaultShearModulusAxialTransverse )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Elastic Shear Modulus Parameter" );

  registerWrapper( viewKeyStruct::defaultC11, &m_defaultC11 )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Default C11 Component of Voigt Stiffness Tensor" );

  registerWrapper( viewKeyStruct::defaultC12, &m_defaultC12 )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Default C12 Component of Voigt Stiffness Tensor" );

  registerWrapper( viewKeyStruct::defaultC13, &m_defaultC13 )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Default C13 Component of Voigt Stiffness Tensor" );

  registerWrapper( viewKeyStruct::defaultC22, &m_defaultC22 )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Default C22 Component of Voigt Stiffness Tensor" );

  registerWrapper( viewKeyStruct::defaultC23, &m_defaultC23 )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Default C23 Component of Voigt Stiffness Tensor" );

  registerWrapper( viewKeyStruct::defaultC33, &m_defaultC33 )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Default C33 Component of Voigt Stiffness Tensor" );

  registerWrapper( viewKeyStruct::defaultC44, &m_defaultC44 )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Default C44 Component of Voigt Stiffness Tensor" );

  registerWrapper( viewKeyStruct::defaultC55, &m_defaultC55 )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Default C55 Component of Voigt Stiffness Tensor" );

  registerWrapper( viewKeyStruct::defaultC66, &m_defaultC66 )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Default C66 Component of Voigt Stiffness Tensor" );

  registerWrapper( viewKeyStruct::c11, &m_c11 )->
    setApplyDefaultValue( -1 )->
    setDescription( "Elastic Stiffness Field C11" );

  registerWrapper( viewKeyStruct::c12, &m_c12 )->
    setApplyDefaultValue( -1 )->
    setDescription( "Elastic Stiffness Field C12" );

  registerWrapper( viewKeyStruct::c13, &m_c13 )->
    setApplyDefaultValue( -1 )->
    setDescription( "Elastic Stiffness Field C13" );

  registerWrapper( viewKeyStruct::c22, &m_c22 )->
    setApplyDefaultValue( -1 )->
    setDescription( "Elastic Stiffness Field C22" );

  registerWrapper( viewKeyStruct::c23, &m_c23 )->
    setApplyDefaultValue( -1 )->
    setDescription( "Elastic Stiffness Field C23" );

  registerWrapper( viewKeyStruct::c33, &m_c33 )->
    setApplyDefaultValue( -1 )->
    setDescription( "Elastic Stiffness Field C33" );

  registerWrapper( viewKeyStruct::c44, &m_c44 )->
    setApplyDefaultValue( -1 )->
    setDescription( "Elastic Stiffness Field C44" );

  registerWrapper( viewKeyStruct::c55, &m_c55 )->
    setApplyDefaultValue( -1 )->
    setDescription( "Elastic Stiffness Field C55" );

  registerWrapper( viewKeyStruct::c66, &m_c66 )->
    setApplyDefaultValue( -1 )->
    setDescription( "Elastic Stiffness Field C66" );
}

ElasticOrthotropic::~ElasticOrthotropic()
{}

void ElasticOrthotropic::postProcessInput()
{

  real64 & c11 = getReference< real64 >( viewKeyStruct::defaultC11 );
  real64 & c12 = getReference< real64 >( viewKeyStruct::defaultC12 );
  real64 & c13 = getReference< real64 >( viewKeyStruct::defaultC13 );
  real64 & c22 = getReference< real64 >( viewKeyStruct::defaultC22 );
  real64 & c23 = getReference< real64 >( viewKeyStruct::defaultC23 );
  real64 & c33 = getReference< real64 >( viewKeyStruct::defaultC33 );
  real64 & c44 = getReference< real64 >( viewKeyStruct::defaultC44 );
  real64 & c55 = getReference< real64 >( viewKeyStruct::defaultC55 );
  real64 & c66 = getReference< real64 >( viewKeyStruct::defaultC66 );

/**
  real64 const Et = m_defaultYoungsModulusTransverse;
  real64 const Ea = m_defaultYoungsModulusAxial;
  real64 const Nut = m_defaultPoissonTransverse;
  real64 const Nuat = m_defaultPoissonAxialTransverse;
  real64 const Gat = m_defaultShearModulusAxialTransverse;
  real64 const Nuta = Nuat * ( Et/Ea );

  real64 const delta = ( 1 + Nut ) * ( 1 - Nut - 2 * Nuta * Nuat ) / ( Et * Et * Ea );
  real64 const c11Default = ( 1.0 - Nuta * Nuat ) / ( Ea * Et * delta );
  real64 const c13Default = Nuat * ( 1.0 + Nut ) / ( Ea * Et * delta );
  real64 const c33Default = ( 1 - Nut * Nut ) / ( Et * Et * delta );
  real64 const c44Default = Gat;
  real64 const c66Default = 0.5 * Et / ( 1 + Nut );
*/

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c11 )->
    setApplyDefaultValue( c11 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c12 )->
    setApplyDefaultValue( c12 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c13 )->
    setApplyDefaultValue( c13 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c22 )->
    setApplyDefaultValue( c22 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c23 )->
    setApplyDefaultValue( c23 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c33 )->
    setApplyDefaultValue( c33 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c44 )->
    setApplyDefaultValue( c44 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c55 )->
    setApplyDefaultValue( c55 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c66 )->
    setApplyDefaultValue( c66 );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, ElasticOrthotropic, string const &, Group * const )
}
} /* namespace geosx */
