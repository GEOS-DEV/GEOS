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
 *  @file ElasticTransverseIsotropicPressureDependent.cpp
 */

#include "ElasticTransverseIsotropicPressureDependent.hpp"

namespace geos
{
using namespace dataRepository;

namespace constitutive
{

ElasticTransverseIsotropicPressureDependent::ElasticTransverseIsotropicPressureDependent( string const & name, Group * const parent ):
  ElasticTransverseIsotropic( name, parent ),
  m_defaultYoungModulusTransversePressureDerivative(),
  m_defaultYoungModulusAxialPressureDerivative(),
  m_defaultShearModulusAxialTransversePressureDerivative(),
  m_dc11dp(),
  m_dc13dp(),
  m_dc33dp(),
  m_dc44dp(),
  m_dc66dp()
{
  registerWrapper( viewKeyStruct::defaultYoungModulusTransversePressureDerivativeString(), &m_defaultYoungModulusTransversePressureDerivative ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Transverse Young's modulus pressure derivative" );

  registerWrapper( viewKeyStruct::defaultYoungModulusAxialPressureDerivativeString(), &m_defaultYoungModulusAxialPressureDerivative ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Axial Young's modulus pressure derivative" );

  registerWrapper( viewKeyStruct::defaultShearModulusAxialTransversePressureDerivativeString(), &m_defaultShearModulusAxialTransversePressureDerivative ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Axial transverse shear modulus pressure derivative" );

  registerWrapper< real64 >( viewKeyStruct::defaultdC11dpString()).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Elastic Stiffness Field C11 pressure derivative" );

  registerWrapper< real64 >( viewKeyStruct::defaultdC13dpString()).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Elastic Stiffness Field C13 pressure derivative" );

  registerWrapper< real64 >( viewKeyStruct::defaultdC33dpString()).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Elastic Stiffness Field C33 pressure derivative" );

 registerWrapper< real64 >( viewKeyStruct::defaultdC44dpString()).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Elastic Stiffness Field C44 pressure derivative" );

registerWrapper< real64 >( viewKeyStruct::defaultdC66dpString()).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Elastic Stiffness Field C66 pressure derivative" );

  registerWrapper( viewKeyStruct::dc11dpString(), &m_dc11dp ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C11 pressure derivative" );

  registerWrapper( viewKeyStruct::dc13dpString(), &m_dc13dp ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C13 pressure derivative" );

  registerWrapper( viewKeyStruct::dc33dpString(), &m_dc33dp ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C33 pressure derivative" );

  registerWrapper( viewKeyStruct::dc44dpString(), &m_dc44dp ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C44 pressure derivative" );

  registerWrapper( viewKeyStruct::dc66dpString(), &m_dc66dp ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C66 pressure derivative" );
}

ElasticTransverseIsotropicPressureDependent::~ElasticTransverseIsotropicPressureDependent()
{}

void ElasticTransverseIsotropicPressureDependent::postProcessInput()
{
  ElasticTransverseIsotropic::postProcessInput();
  
  real64 & dc11dp  = getReference< real64 >( viewKeyStruct::defaultdC11dpString() );
  real64 & dc13dp  = getReference< real64 >( viewKeyStruct::defaultdC13dpString() );
  real64 & dc33dp  = getReference< real64 >( viewKeyStruct::defaultdC33dpString() );
  real64 & dc44dp  = getReference< real64 >( viewKeyStruct::defaultdC44dpString() );
  real64 & dc66dp  = getReference< real64 >( viewKeyStruct::defaultdC66dpString() );

  real64 & dEtdp = m_defaultYoungModulusTransversePressureDerivative;
  real64 & dEadp = m_defaultYoungModulusAxialPressureDerivative;
  real64 & dGatdp = m_defaultShearModulusAxialTransversePressureDerivative;

  real64 Nuat = m_defaultPoissonRatioAxialTransverse;
  real64 Nut = m_defaultPoissonRatioTransverse;

  // CC: could the pressure derivatives every be negative? also should I optionally allow poisson ratios to change?
  if( dEtdp > 0.0 && dEadp > 0.0 && dGatdp )
  {
    dc11dp = ( Nuat * Nuat - 1 ) * dEtdp / ( ( Nut + 1 ) * ( 2*Nuat * Nuat + Nut - 1 ) );
    dc13dp = -( Nuat * dEtdp ) / ( 2 * Nuat * Nuat + Nut - 1 );
    dc33dp = ( ( Nut * Nut - 1 ) * dEadp ) / ( ( Nut + 1 ) * ( 2 * Nuat * Nuat + Nut - 1 ) );
    dc44dp = 2 * dGatdp;
    dc66dp = dEtdp / ( Nut + 1 );
  } 

  this->getWrapper< real64 >( viewKeyStruct::dc11dpString() ).
    setApplyDefaultValue( dc11dp );

  this->getWrapper< real64 >( viewKeyStruct::dc13dpString() ).
    setApplyDefaultValue( dc13dp );

  this->getWrapper< real64 >( viewKeyStruct::dc33dpString() ).
    setApplyDefaultValue( dc33dp );

  this->getWrapper< real64 >( viewKeyStruct::dc44dpString() ).
    setApplyDefaultValue( dc44dp );

  this->getWrapper< real64 >( viewKeyStruct::dc66dpString() ).
    setApplyDefaultValue( dc66dp );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ElasticTransverseIsotropicPressureDependent, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
