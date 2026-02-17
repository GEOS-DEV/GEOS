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
  m_damage(),
  m_previousLambda(),
  m_previousPlasticStrain()
{
  // register constants
  registerWrapper( viewKeyStruct::thicknessString(), &m_thickness ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Cohesive zone thickness" );

  registerWrapper( viewKeyStruct::bulkModulusString(), &m_bulkModulus ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Bulk modulus" );

  registerWrapper( viewKeyStruct::bulkModulusAString(), &m_bulkModulusA ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Bulk modulus A" );

  registerWrapper( viewKeyStruct::bulkModulusBString(), &m_bulkModulusB ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Bulk modulus B" );

  registerWrapper( viewKeyStruct::bulkModulusT0String(), &m_bulkModulusT0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Bulk modulus T0" );

  registerWrapper( viewKeyStruct::shearModulusString(), &m_shearModulus ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear modulus" );

  registerWrapper( viewKeyStruct::shearModulusAString(), &m_shearModulusA ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear modulus A" );

  registerWrapper( viewKeyStruct::shearModulusBString(), &m_shearModulusB ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear modulus B" );

  registerWrapper( viewKeyStruct::shearModulusT0String(), &m_shearModulusT0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear modulus T0" );

  registerWrapper( viewKeyStruct::yieldStrength0String(), &m_yieldStrength0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Initial yield strength" );

  registerWrapper( viewKeyStruct::yieldStrengthAString(), &m_yieldStrengthA ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Initial yield strength A" );

  registerWrapper( viewKeyStruct::yieldStrengthBString(), &m_yieldStrengthB ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Initial yield strength B" );

  registerWrapper( viewKeyStruct::yieldStrengthT0String(), &m_yieldStrengthT0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Initial yield strength T0" );

  registerWrapper( viewKeyStruct::r0String(), &m_r0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "r0 softening parameter" );

  registerWrapper( viewKeyStruct::r0AString(), &m_r0A ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "r0 softening parameter A" );

  registerWrapper( viewKeyStruct::r0BString(), &m_r0B ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "r0 softening parameter B" );

  registerWrapper( viewKeyStruct::r0T0String(), &m_r0T0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "r0 softening parameter T0" );

  registerWrapper( viewKeyStruct::r1String(), &m_r1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "r1 softening parameter" );

  registerWrapper( viewKeyStruct::r2String(), &m_r2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "r2 softening parameter" );

  registerWrapper( viewKeyStruct::GrString(), &m_Gr ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Gr hardening parameter" );

  registerWrapper( viewKeyStruct::GrAString(), &m_GrA ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Gr hardening parameter A" );

  registerWrapper( viewKeyStruct::GrBString(), &m_GrB ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Gr hardening parameter B" );

  registerWrapper( viewKeyStruct::GrT0String(), &m_GrT0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Gr hardening parameter T0" );

  registerWrapper( viewKeyStruct::maxStretchString(), &m_maxStretch ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum stretch" );

  registerWrapper( viewKeyStruct::maxStretchAString(), &m_maxStretchA ).
  setInputFlag( InputFlags::REQUIRED ).
  setDescription( "Maximum stretch A" );

  registerWrapper( viewKeyStruct::maxStretchBString(), &m_maxStretchB ).
  setInputFlag( InputFlags::REQUIRED ).
  setDescription( "Maximum stretch B" );

  registerWrapper( viewKeyStruct::maxStretchT0String(), &m_maxStretchT0 ).
  setInputFlag( InputFlags::REQUIRED ).
  setDescription( "Maximum stretch T0" );

  // register fields
  registerWrapper( viewKeyStruct::damageString(), &m_damage ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Damage" );

  registerWrapper( viewKeyStruct::temperatureString(), &m_temperature ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Temperature" );

  registerWrapper( viewKeyStruct::previousLambdaString(), &m_previousLambda ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Previous stretch" );

  registerWrapper( viewKeyStruct::previousPlasticStrainString(), &m_previousPlasticStrain ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Previous plastic strain" );

}


PolymerCohesiveZone::~PolymerCohesiveZone()
{}


void PolymerCohesiveZone::allocateConstitutiveData( dataRepository::Group & parent,
                                                      localIndex const numConstitutivePointsPerParentIndex )
{
  CohesiveZoneBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_damage.resize( 0 );
  m_temperature.resize( 0 );
  m_previousLambda.resize( 0 );
  m_previousPlasticStrain( 0 );
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
