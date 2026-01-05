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
 *  @file CoupledCohesiveZone.cpp
 */

#include "CoupledCohesiveZone.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

CoupledCohesiveZone::CoupledCohesiveZone( string const & name, Group * const parent ):
  CohesiveZoneBase( name, parent ),
  m_characteristicNormalDisplacement( 1.0 ),
  m_characteristicTangentialDisplacement( 1.0 ),
  m_maxNormalStress( 1.0 ),
  m_maxShearStress( 1.0 ),
  m_maxNormalDisplacement( DBL_MAX ),
  m_maxTangentialDisplacement( DBL_MAX ),
  m_damage()
{
  // register constants
  registerWrapper( viewKeyStruct::characteristicNormalDisplacementString(), &m_characteristicNormalDisplacement ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Characteristic normal displacement" );

  registerWrapper( viewKeyStruct::characteristicTangentialDisplacementString(), &m_characteristicTangentialDisplacement ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Characteristic tangential displacement" );

  registerWrapper( viewKeyStruct::maxNormalStressString(), &m_maxNormalStress ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum normal stress" );

  registerWrapper( viewKeyStruct::maxShearStressString(), &m_maxShearStress ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum shear stress" );  

  registerWrapper( viewKeyStruct::maxNormalDisplacementString(), &m_maxNormalDisplacement ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_maxNormalDisplacement ).
    setDescription( "Maximum normal displacement" );

  registerWrapper( viewKeyStruct::maxTangentialDisplacementString(), &m_maxTangentialDisplacement ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_maxTangentialDisplacement ).
    setDescription( "Maximum tangential displacement" );

  // register fields
  registerWrapper( viewKeyStruct::damageString(), &m_damage ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Damage" );
}


CoupledCohesiveZone::~CoupledCohesiveZone()
{}


void CoupledCohesiveZone::allocateConstitutiveData( dataRepository::Group & parent,
                                                      localIndex const numConstitutivePointsPerParentIndex )
{
  CohesiveZoneBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_damage.resize( 0 );
}


void CoupledCohesiveZone::postInputInitialization()
{
  CohesiveZoneBase::postInputInitialization();

  GEOS_THROW_IF( m_characteristicNormalDisplacement < 0.0, "Characteristic normal displacement must be a positive number.", InputError );
  GEOS_THROW_IF( m_characteristicTangentialDisplacement < 0.0, "Characteristic tangential displacement must be a positive number.", InputError );
  GEOS_THROW_IF( m_maxNormalStress < 0.0, "Maximum normal stress must be a positive number.", InputError );
  GEOS_THROW_IF( m_maxShearStress < 0.0, "Maximum shear stress must be a positive number.", InputError );
  GEOS_THROW_IF( m_maxNormalDisplacement < 0.0, "Maximum normal displacement must be a positive number.", InputError );
  GEOS_THROW_IF( m_maxTangentialDisplacement < 0.0, "Maximum tangential displacement must be a positive number.", InputError );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, CoupledCohesiveZone, std::string const &, Group * const )
}
} /* namespace geos */
