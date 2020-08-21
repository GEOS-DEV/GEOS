/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file SolidBase.cpp
 */

#include "SolidBase.hpp"

namespace geosx
{
using namespace dataRepository;

namespace constitutive
{
  
SolidBase::SolidBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_newStress( 0, 0, 6 ),
  m_oldStress( 0, 0, 6 ),
  m_density()
{
  registerWrapper( viewKeyStruct::stressString, &m_newStress )->
    setPlotLevel( PlotLevel::LEVEL_0 )->
    setApplyDefaultValue( 0 )-> // default to zero initial stress
    setDescription( "Current Material Stress" );
    
  registerWrapper( viewKeyStruct::stressString, &m_oldStress )->
    setApplyDefaultValue( 0 )-> // default to zero initial stress
    setDescription( "Previous Material Stress" );

  registerWrapper( viewKeyStruct::densityString, &m_density )->
    setApplyDefaultValue( -1 )-> // will be overwritten
    setDescription( "Material Density" );
  
  registerWrapper( viewKeyStruct::defaultDensityString, &m_defaultDensity )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Default Material Density" );
}


SolidBase::~SolidBase()
{}


void SolidBase::PostProcessInput()
{
  this->getWrapper< array2d< real64 > >( viewKeyStruct::densityString )->
    setApplyDefaultValue( m_defaultDensity );
}


void SolidBase::allocateConstitutiveData( dataRepository::Group * const parent,
                                          localIndex const numConstitutivePointsPerParentIndex )
{
  m_density.resize( 0, numConstitutivePointsPerParentIndex );
  m_newStress.resize( 0, numConstitutivePointsPerParentIndex, 6 );
  m_oldStress.resize( 0, numConstitutivePointsPerParentIndex, 6 );
  
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}


} /* namespace constitutive */
} /* namespace geosx */
