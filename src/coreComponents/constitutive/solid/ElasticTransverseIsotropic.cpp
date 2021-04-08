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
 *  @file ElasticTransverseIsotropic.cpp
 */

#include "ElasticTransverseIsotropic.hpp"

namespace geosx
{
using namespace dataRepository;

namespace constitutive
{

ElasticTransverseIsotropic::ElasticTransverseIsotropic( string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_defaultYoungsModulusTransverse(),
  m_defaultYoungsModulusAxial(),
  m_defaultPoissonTransverse(),
  m_defaultPoissonAxialTransverse(),
  m_defaultShearModulusAxialTransverse(),
  m_c11(),
  m_c13(),
  m_c33(),
  m_c44(),
  m_c66()
{
  registerWrapper( viewKeyStruct::defaultYoungsModulusTransverseString(), &m_defaultYoungsModulusTransverse ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Elastic Bulk Modulus Parameter" );

  registerWrapper( viewKeyStruct::defaultYoungsModulusAxialString(), &m_defaultYoungsModulusAxial ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Elastic Shear Modulus Parameter" );

  registerWrapper( viewKeyStruct::defaultPoissonRatioTransverseString(), &m_defaultPoissonTransverse ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Elastic Shear Modulus Parameter" );

  registerWrapper( viewKeyStruct::defaultPoissonRatioAxialTransverseString(), &m_defaultPoissonAxialTransverse ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Elastic Shear Modulus Parameter" );

  registerWrapper( viewKeyStruct::defaultShearModulusAxialTransverseString(), &m_defaultShearModulusAxialTransverse ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Elastic Shear Modulus Parameter" );


  registerWrapper( viewKeyStruct::c11String(), &m_c11 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Bulk Modulus Field" );

  registerWrapper( viewKeyStruct::c13String(), &m_c13 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Bulk Modulus Field" );

  registerWrapper( viewKeyStruct::c33String(), &m_c33 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Bulk Modulus Field" );

  registerWrapper( viewKeyStruct::c44String(), &m_c44 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Bulk Modulus Field" );

  registerWrapper( viewKeyStruct::c66String(), &m_c66 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Bulk Modulus Field" );
}


ElasticTransverseIsotropic::~ElasticTransverseIsotropic()
{}


void ElasticTransverseIsotropic::postProcessInput()
{
  SolidBase::postProcessInput();

  real64 const Et = m_defaultYoungsModulusTransverse;
  real64 const Ea = m_defaultYoungsModulusAxial;
  real64 const Nut = m_defaultPoissonTransverse;
  real64 const Nuat = m_defaultPoissonAxialTransverse;
  real64 const Gat = m_defaultShearModulusAxialTransverse;
  real64 const Nuta = Nuat * ( Et/Ea );

  real64 const delta = ( 1 + Nut ) * ( 1 - Nut - 2 * Nuta * Nuat ) / ( Et * Et * Ea );
  real64 const c11Default = ( 1.0 - Nuta * Nuat ) / ( Ea * Et * delta );
  real64 const c13Default = Nuat * ( 1.0 + Nut ) / ( Ea * Et * delta );
  real64 const c33Default = ( 1 - Nut * Nut ) / ( Ea * Ea * delta );
  real64 const c44Default = Gat;
  real64 const c66Default = 0.5 * Ea / ( 1 + Nut );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c11String() ).
    setApplyDefaultValue( c11Default );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c13String() ).
    setApplyDefaultValue( c13Default );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c33String() ).
    setApplyDefaultValue( c33Default );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c44String() ).
    setApplyDefaultValue( c44Default );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c66String() ).
    setApplyDefaultValue( c66Default );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, ElasticTransverseIsotropic, string const &, Group * const )
}
} /* namespace geosx */
