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
 *  @file ElasticAnisotropic.cpp
 */

#include "ElasticAnisotropic.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{

ElasticAnisotropic::ElasticAnisotropic( string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_defaultStiffness( 6, 6 )
{
  registerWrapper( viewKeyStruct::defaultStiffnessString, &m_defaultStiffness )->
    setInputFlag( InputFlags::REQUIRED )->
    setSizedFromParent( 0 )->
    setDescription( "Default Elastic Stiffness Tensor in Voigt notation (6x6 matrix)" );

  m_defaultStiffness.resize( 6, 6 );

  // These are temporary until we figure out how to read in multidim arrays from input.

  registerWrapper( viewKeyStruct::stiffnessString, &m_stiffness )->
    setApplyDefaultValue( 0 )->
    setDescription( "Fully Anisotropic Elastic Stiffness Field in Voigt notation (6x6 matrix)" );

  m_stiffness.resizeDimension< 1, 2 >( 6, 6 );
}

ElasticAnisotropic::~ElasticAnisotropic()
{}

void ElasticAnisotropic::allocateConstitutiveData( dataRepository::Group * const parent,
                                                         localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );

  for( localIndex k=0; k< m_stiffness.size( 0 ); ++k )
  {
    for( localIndex i=0; i<6; ++i )
    {
      for( localIndex j=0; j<6; ++j )
      {
        m_stiffness( k, i, j ) = m_defaultStiffness[i][j];
      }
    }
  }

}

void ElasticAnisotropic::postProcessInput()
{
  SolidBase::postProcessInput();

  // need to get this to work. Need 2d array default value for 3d array.
  // getWrapper<array1d<real64>>(viewKeyStruct::stiffnessString)->setDefaultValue(m_defaultStiffness);
}

void ElasticAnisotropic::setDefaultStiffness( real64 const c[6][6] )
{
  for( int i=0; i<6; ++i )
  {
    for( int j=0; j<6; ++j )
    {
      m_defaultStiffness( i, j ) = c[i][j];
    }
  }
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ElasticAnisotropic, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geosx */
