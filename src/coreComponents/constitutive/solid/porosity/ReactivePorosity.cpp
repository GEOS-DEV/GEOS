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
  registerWrapper( viewKeyStruct::initialVolumeFractionsString(), &m_initialVolumeFractions ).
    setInputFlag( InputFlags::REQUIRED ).
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

  // Hard code for now 
  m_numKineticReactions = 1; 
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

  m_volumeFractions.resize( size, numPts, numKineticReactions );
  m_volumeFractions_n.resize( size, numPts, numKineticReactions );
}

void ReactivePorosity::postInputInitialization()
{
  PorosityBase::postInputInitialization();
  
  integer const numKineticReactions = this->numKineticReactions();

  GEOS_THROW_IF_NE_MSG( m_initialVolumeFractions.size(), numKineticReactions,
                        GEOS_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), viewKeyStruct::initialVolumeFractionsString() ),
                        InputError );

  GEOS_THROW_IF_NE_MSG( m_molarWeights.size(), numKineticReactions,
                        GEOS_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), viewKeyStruct::molarWeightsString() ),
                        InputError );

  GEOS_THROW_IF_NE_MSG( m_mineralDensities.size(), numKineticReactions,
                        GEOS_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), viewKeyStruct::mineralDensitiesString() ),
                        InputError );
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
        m_volumeFractions[ei][q][r] = m_initialVolumeFractions[r]; 
        m_volumeFractions_n[ei][q][r] = m_initialVolumeFractions[r]; 
      }
    }
  }
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactivePorosity, string const &, Group * const )
}
} /* namespace geos */
