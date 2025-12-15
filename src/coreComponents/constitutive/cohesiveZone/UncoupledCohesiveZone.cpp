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
 *  @file UncoupledCohesiveZone.cpp
 */

#include "UncoupledCohesiveZone.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

UncoupledCohesiveZone::UncoupledCohesiveZone( string const & name, Group * const parent ):
  CohesiveZoneBase( name, parent ),
  m_normalForceConstant( 1.0 ),
  m_shearForceConstant( 1.0 ),
  m_maxNormalDisplacement( DBL_MAX ),
  m_maxTangentialDisplacement( DBL_MAX ),
  m_damage()
{
  // register constants
  registerWrapper( viewKeyStruct::normalForceConstantString(), &m_normalForceConstant ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Normal force constant" );

  registerWrapper( viewKeyStruct::shearForceConstantString(), &m_shearForceConstant ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear force constant" );

  registerWrapper( viewKeyStruct::maxNormalDisplacementString(), &m_maxNormalDisplacement ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum normal displacement" );

  registerWrapper( viewKeyStruct::maxTangentialDisplacementString(), &m_maxTangentialDisplacement ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum tangential displacement" );

  // register fields
  registerWrapper( viewKeyStruct::damageString(), &m_damage ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Damage" );
}


UncoupledCohesiveZone::~UncoupledCohesiveZone()
{}


void UncoupledCohesiveZone::allocateConstitutiveData( dataRepository::Group & parent,
                                                      localIndex const numConstitutivePointsPerParentIndex )
{
  CohesiveZoneBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_damage.resize( 0 );
}


void UncoupledCohesiveZone::postInputInitialization()
{
  CohesiveZoneBase::postInputInitialization();

  GEOS_THROW_IF( m_normalForceConstant < 0.0, "Normal force constant must be a positive number.", InputError );
  GEOS_THROW_IF( m_shearForceConstant < 0.0, "Shear force constant must be a positive number.", InputError );
  GEOS_THROW_IF( m_maxNormalDisplacement < 0.0, "Maximum normal displacement must be a positive number.", InputError );
  GEOS_THROW_IF( m_maxTangentialDisplacement < 0.0, "Maximum tangential displacement must be a positive number.", InputError );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, UncoupledCohesiveZone, std::string const &, Group * const )
}
} /* namespace geos */
