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
  m_density(),
  m_defaultDensity( 0 ),
  m_postProcessed( false )
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


void SolidBase::DeliverClone( string const & GEOSX_UNUSED_PARAM( name ),
                              Group * const GEOSX_UNUSED_PARAM( parent ),
                              std::unique_ptr< ConstitutiveBase > & clone ) const
{
  SolidBase * const newConstitutiveRelation = dynamic_cast< SolidBase * >(clone.get());

  newConstitutiveRelation->m_newStress = m_newStress;
  newConstitutiveRelation->m_oldStress = m_oldStress;
  newConstitutiveRelation->m_density = m_density;
  newConstitutiveRelation->m_defaultDensity = m_defaultDensity;
}


void SolidBase::AllocateConstitutiveData( dataRepository::Group * const parent,
                                          localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
  
  localIndex const numElems = parent->size();
  this->resize( numElems );
  
  m_density.resize( numElems, numConstitutivePointsPerParentIndex );
  m_density.setValues< serialPolicy >( m_defaultDensity );

  m_newStress.resize( numElems, numConstitutivePointsPerParentIndex, 6 );
  m_oldStress.resize( numElems, numConstitutivePointsPerParentIndex, 6 );
}


} /* namespace constitutive */
} /* namespace geosx */
