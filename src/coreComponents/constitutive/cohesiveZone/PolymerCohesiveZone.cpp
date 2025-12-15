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
 *  @file PolymerCohesiveZone.cpp
 */

#include "PolymerCohesiveZone.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

PolymerCohesiveZone::PolymerCohesiveZone( string const & name, Group * const parent ):
  CohesiveZoneBase( name, parent ),
  m_thickness( 1.0 ),
  m_bulkModulus( 1.0 ),
  m_shearModulus( 1.0 ),
  m_yieldStrength0( 1.0 ),
  m_r0( 0.0 ),
  m_r1( 1.0 ),
  m_r2( 0.0 ),
  m_Gr( 0.0 ),
  m_maxStretch( DBL_MAX ),
  m_damage()
{
  // register constants
  registerWrapper( viewKeyStruct::thicknessString(), &m_thickness ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Cohesive zone thickness" );

  registerWrapper( viewKeyStruct::bulkModulusString(), &m_bulkModulus ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Bulk modulus" );

  registerWrapper( viewKeyStruct::shearModulusString(), &m_shearModulus ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear modulus" );

  registerWrapper( viewKeyStruct::yieldStrength0String(), &m_yieldStrength0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Initial yield strength" );

  registerWrapper( viewKeyStruct::r0String(), &m_r0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "r0 softening parameter" );

  registerWrapper( viewKeyStruct::r1String(), &m_r1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "r1 softening parameter" );

  registerWrapper( viewKeyStruct::r2String(), &m_r2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "r2 softening parameter" );

  registerWrapper( viewKeyStruct::GrString(), &m_Gr ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Gr hardening parameter" );

  registerWrapper( viewKeyStruct::maxStretchString(), &m_maxStretch ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum stretch" );

  // register fields
  registerWrapper( viewKeyStruct::damageString(), &m_damage ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Damage" );
}


PolymerCohesiveZone::~PolymerCohesiveZone()
{}


void PolymerCohesiveZone::allocateConstitutiveData( dataRepository::Group & parent,
                                                      localIndex const numConstitutivePointsPerParentIndex )
{
  CohesiveZoneBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_damage.resize( 0 );
}


void PolymerCohesiveZone::postInputInitialization()
{
  CohesiveZoneBase::postInputInitialization();

  GEOS_THROW_IF( m_thickness < 0.0, "Thickness must be a positive number.", InputError );
  GEOS_THROW_IF( m_bulkModulus < 0.0, "Bulk modulus must be a positive number.", InputError );
  GEOS_THROW_IF( m_shearModulus < 0.0, "Shear modulus must be a positive number.", InputError );
  GEOS_THROW_IF( m_yieldStrength0 < 0.0, "Initial yield strength must be a positive number.", InputError );
  GEOS_THROW_IF( m_maxStretch < 0.0, "Maximum stretch must be a positive number.", InputError );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, PolymerCohesiveZone, std::string const &, Group * const )
}
} /* namespace geos */
