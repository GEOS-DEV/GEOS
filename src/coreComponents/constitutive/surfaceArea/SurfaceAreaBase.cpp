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
 * @file SurfaceAreaBase.cpp
 */

#include "constitutive/surfaceArea/SurfaceAreaBase.hpp"
#include "constitutive/surfaceArea/SurfaceAreaFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

SurfaceAreaBase::SurfaceAreaBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_numKineticReactions( 0 )
{
  registerWrapper( viewKeyStruct::defaultInitialSurfaceAreaString(), &m_defaultInitialSurfaceArea ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default initial reactive surface area of each mineral. "
                    "Its size defines the number of kinetic reactions handled by this model." );

  registerField< fields::surfacearea::surfaceArea >( &m_surfaceArea );
  registerField< fields::surfacearea::initialSurfaceArea >( &m_initialSurfaceArea );
}

std::unique_ptr< ConstitutiveBase > SurfaceAreaBase::deliverClone( string const & name, Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = ConstitutiveBase::deliverClone( name, parent );

  SurfaceAreaBase & newConstitutiveRelation = dynamicCast< SurfaceAreaBase & >( *clone );

  newConstitutiveRelation.m_numKineticReactions = m_numKineticReactions;

  return clone;
}

void SurfaceAreaBase::postInputInitialization()
{
  ConstitutiveBase::postInputInitialization();

  m_numKineticReactions = LvArray::integerConversion< integer >( m_defaultInitialSurfaceArea.size() );
}

void SurfaceAreaBase::allocateConstitutiveData( Group & parent,
                                                localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  resizeFields( parent.size(), numConstitutivePointsPerParentIndex );
}

void SurfaceAreaBase::resizeFields( localIndex const size, localIndex const numPts )
{
  integer const numKineticReactions = this->numKineticReactions();

  m_initialSurfaceArea.resize( size, numPts, numKineticReactions );
  m_surfaceArea.resize( size, numPts, numKineticReactions );

  // Apply the default values here so that they can still be overwritten element by element
  // with a FieldSpecification on <modelName>_initialSurfaceArea before initializeState is called.
  for( localIndex ei = 0; ei < size; ++ei )
  {
    for( localIndex q = 0; q < numPts; ++q )
    {
      for( integer r = 0; r < numKineticReactions; ++r )
      {
        m_initialSurfaceArea[ei][q][r] = m_defaultInitialSurfaceArea[r];
        m_surfaceArea[ei][q][r] = m_defaultInitialSurfaceArea[r];
      }
    }
  }
}

void SurfaceAreaBase::initializeState() const
{
  integer const numKineticReactions = this->numKineticReactions();

  for( localIndex ei = 0; ei < m_surfaceArea.size( 0 ); ++ei )
  {
    for( localIndex q = 0; q < m_surfaceArea.size( 1 ); ++q )
    {
      for( integer r = 0; r < numKineticReactions; ++r )
      {
        m_surfaceArea[ei][q][r] = m_initialSurfaceArea[ei][q][r];
      }
    }
  }
}

}
} /* namespace geos */
