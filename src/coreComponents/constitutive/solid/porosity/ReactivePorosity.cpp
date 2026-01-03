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
 * @file ReactivePorosity.cpp
 */

#include "ReactivePorosity.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

ReactivePorosity::ReactivePorosity( string const & name, Group * const parent ):
  PorosityBase( name, parent )
{
  registerWrapper( viewKeyStruct::defaultInitialVolumeFractionsString(), &m_defaultInitialVolumeFractions ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default initial volume fractions" );

  registerWrapper( viewKeyStruct::initialVolumeFractionsString(), &m_initialVolumeFractions ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Initial volume fractions" );

  registerWrapper( viewKeyStruct::volumeFractionsString(), &m_volumeFractions ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Current volume fractions" );

  registerWrapper( viewKeyStruct::volumeFractions_nString(), &m_volumeFractions_n ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Volume fractions at last time step" );

  registerWrapper( viewKeyStruct::molarWeightsString(), &m_molarWeights ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Mineral molar weights" );

  registerWrapper( viewKeyStruct::mineralDensitiesString(), &m_mineralDensities ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Mineral densities" );
}

std::unique_ptr< ConstitutiveBase > ReactivePorosity::deliverClone( string const & name, Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = ConstitutiveBase::deliverClone( name, parent );

  ReactivePorosity & newConstitutiveRelation = dynamicCast< ReactivePorosity & >( *clone );

  newConstitutiveRelation.m_numKineticReactions = m_numKineticReactions;

  return clone;
}

void ReactivePorosity::postInputInitialization()
{
  PorosityBase::postInputInitialization();

  GEOS_THROW_IF_NE_MSG( m_molarWeights.size(), m_defaultInitialVolumeFractions.size(),
                        GEOS_FMT( "{}: mismatch in number of components in porosity model attribute '{}' and '{}",
                                  getFullName(), viewKeyStruct::molarWeightsString(), viewKeyStruct::initialVolumeFractionsString() ),
                        InputError );

  GEOS_THROW_IF_NE_MSG( m_mineralDensities.size(), m_defaultInitialVolumeFractions.size(),
                        GEOS_FMT( "{}: mismatch in number of components in porosity model attribute '{}' and '{}",
                                  getFullName(), viewKeyStruct::mineralDensitiesString(), viewKeyStruct::initialVolumeFractionsString() ),
                        InputError );

  m_numKineticReactions = m_defaultInitialVolumeFractions.size();
}

void ReactivePorosity::allocateConstitutiveData( dataRepository::Group & parent,
                                                 localIndex const numConstitutivePointsPerParentIndex )
{
  PorosityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  resizeFields( parent.size(), numConstitutivePointsPerParentIndex );

}

void ReactivePorosity::resizeFields( localIndex const size, localIndex const numPts )
{
  integer const numKineticReactions = this->numKineticReactions();

  m_initialVolumeFractions.resize( size, numPts, numKineticReactions );
  m_volumeFractions.resize( size, numPts, numKineticReactions );
  m_volumeFractions_n.resize( size, numPts, numKineticReactions );
}

void ReactivePorosity::saveConvergedState() const
{
  PorosityBase::saveConvergedState();

  m_volumeFractions_n.setValues< parallelDevicePolicy<> >( m_volumeFractions.toViewConst() );
}

void ReactivePorosity::initializeState() const
{
  integer const numKineticReactions = this->numKineticReactions();

  for( localIndex ei = 0; ei < m_newPorosity.size( 0 ); ++ei )
  {
    for( localIndex q = 0; q < m_newPorosity.size( 1 ); ++q )
    {
      m_newPorosity[ei][q] = m_referencePorosity[ei];
    }
  }

  PorosityBase::initializeState();

  for( localIndex ei = 0; ei < m_volumeFractions.size( 0 ); ++ei )
  {
    for( localIndex q = 0; q < m_volumeFractions.size( 1 ); ++q )
    {
      for( integer r = 0; r < numKineticReactions; ++r )
      {
        m_volumeFractions[ei][q][r] = m_defaultInitialVolumeFractions[r];
        m_initialVolumeFractions[ei][q][r] = m_defaultInitialVolumeFractions[r];
        m_volumeFractions_n[ei][q][r] = m_defaultInitialVolumeFractions[r];
      }
    }
  }
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactivePorosity, string const &, Group * const )
}
} /* namespace geos */
