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
 *  @file ElasticCubic.cpp
 */

#include "ElasticCubic.hpp"

namespace geos
{
using namespace dataRepository;

namespace constitutive
{

ElasticCubic::ElasticCubic( string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_c11(),
  m_c12(),
  m_c44(),
  m_bulkModulus(),
  m_effectiveShearModulus(),
  m_materialDirection()
{
  registerWrapper< real64 >( viewKeyStruct::defaultC11String() ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Stiffness Parameter C11" );

  registerWrapper< real64 >( viewKeyStruct::defaultC12String() ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Stiffness Parameter C12" );

  registerWrapper< real64 >( viewKeyStruct::defaultC44String() ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Stiffness Parameter C44" );

  registerWrapper( viewKeyStruct::c11String(), &m_c11 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C11" );

  registerWrapper( viewKeyStruct::c12String(), &m_c12 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C12" );

  registerWrapper( viewKeyStruct::c44String(), &m_c44 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C44" );

  registerWrapper( viewKeyStruct::bulkModulusString(), &m_bulkModulus ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Effective bulk modulus for stress control and wavespeed calculations" );
  
  registerWrapper( viewKeyStruct::effectiveShearModulusString(), &m_effectiveShearModulus ).
    setInputFlag( InputFlags::FALSE).
    setDescription( "Effective shear modulus for stress control and wavespeed calculations");

  registerWrapper( viewKeyStruct::materialDirectionString(), &m_materialDirection ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Material direction" );
}

ElasticCubic::~ElasticCubic()
{}

void ElasticCubic::allocateConstitutiveData( dataRepository::Group & parent, 
                                                           localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_materialDirection.resize( 0, 3, 3 );
}

void ElasticCubic::postProcessInput()
{
  SolidBase::postProcessInput();

  real64 & c11  = getReference< real64 >( viewKeyStruct::defaultC11String() );
  real64 & c12  = getReference< real64 >( viewKeyStruct::defaultC12String() );
  real64 & c44  = getReference< real64 >( viewKeyStruct::defaultC44String() );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c11String() ).
    setApplyDefaultValue( c11 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c12String() ).
    setApplyDefaultValue( c12 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c44String() ).
    setApplyDefaultValue( c44 );

  real64 K = (c11+2*c12)/2;
  real64 Geff = fmax( c44, 0.5*(c11-c12) );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::bulkModulusString() ).
    setApplyDefaultValue( K );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::effectiveShearModulusString() ).
    setApplyDefaultValue( Geff );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ElasticCubic, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
