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
 * @file BlackOilFluid.cpp
 */

#include "BlackOilFluid.hpp"

#include "codingUtilities/Utilities.hpp"
#include "common/Path.hpp"

#include "pvt/pvt.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

BlackOilFluid::BlackOilFluid( string const & name, Group * const parent )
  : MultiFluidPVTPackageWrapper( name, parent )
{
  getWrapperBase( viewKeyStruct::componentMolarWeightString() ).setInputFlag( InputFlags::REQUIRED );
  getWrapperBase( viewKeyStruct::phaseNamesString() ).setInputFlag( InputFlags::REQUIRED );

  registerWrapper( viewKeyStruct::surfaceDensitiesString(), &m_surfaceDensities ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "List of surface densities for each phase" );

  registerWrapper( viewKeyStruct::tableFilesString(), &m_tableFiles ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "List of filenames with input PVT tables" );

}

BlackOilFluid::~BlackOilFluid()
{}

std::unique_ptr< ConstitutiveBase >
BlackOilFluid::deliverClone( string const & name,
                             Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase >
  clone = MultiFluidPVTPackageWrapper::deliverClone( name, parent );
  BlackOilFluid & fluid = dynamicCast< BlackOilFluid & >( *clone );

  fluid.createFluid();
  return clone;
}

void BlackOilFluid::postProcessInput()
{
  // TODO maybe use different names?
  m_componentNames = m_phaseNames;

  MultiFluidPVTPackageWrapper::postProcessInput();

  localIndex const NP = numFluidPhases();

  #define BOFLUID_CHECK_INPUT_LENGTH( data, expected, attr ) \
    if( LvArray::integerConversion< localIndex >((data).size()) != LvArray::integerConversion< localIndex >( expected )) \
    { \
      GEOSX_ERROR( "BlackOilFluid: invalid number of entries in " \
                   << (attr) << " attribute (" \
                   << (data).size() << "given, " \
                   << (expected) << " expected)" ); \
    }

  BOFLUID_CHECK_INPUT_LENGTH( m_surfaceDensities, NP, viewKeyStruct::surfaceDensitiesString() )
  BOFLUID_CHECK_INPUT_LENGTH( m_tableFiles, NP, viewKeyStruct::surfaceDensitiesString() )

  #undef BOFLUID_CHECK_INPUT_LENGTH
}

void BlackOilFluid::createFluid()
{
  std::vector< pvt::PHASE_TYPE > phases( m_phaseTypes.begin(), m_phaseTypes.end() );
  std::vector< string > tableFiles( m_tableFiles.begin(), m_tableFiles.end() );
  std::vector< double > densities( m_surfaceDensities.begin(), m_surfaceDensities.end() );
  std::vector< double > molarWeights( m_componentMolarWeight.begin(), m_componentMolarWeight.end() );

  m_fluid = pvt::MultiphaseSystemBuilder::buildLiveOil( phases, tableFiles, densities, molarWeights );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, BlackOilFluid, string const &, Group * const )
} // namespace constitutive

} // namespace geosx
