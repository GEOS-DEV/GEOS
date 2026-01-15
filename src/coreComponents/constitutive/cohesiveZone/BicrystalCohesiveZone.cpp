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
 *  @file BicrystalCohesiveZone.cpp
 */

#include "BicrystalCohesiveZone.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

BicrystalCohesiveZone::BicrystalCohesiveZone( string const & name, Group * const parent ):
  CoupledCohesiveZone( name, parent )
{
  // // register constants
  // registerWrapper( viewKeyStruct::characteristicNormalDisplacementString(), &m_characteristicNormalDisplacement ).
  //   setInputFlag( InputFlags::REQUIRED ).
  //   setDescription( "Characteristic normal displacement" );

  // register fields
  registerWrapper( viewKeyStruct::misorientationString(), &m_misorientation ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Misorientation rotation matrix" );
}


BicrystalCohesiveZone::~BicrystalCohesiveZone()
{}


void BicrystalCohesiveZone::allocateConstitutiveData( dataRepository::Group & parent,
                                                      localIndex const numConstitutivePointsPerParentIndex )
{
  CoupledCohesiveZone::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_misorientation.resize( 0, 3, 3 );
}


void BicrystalCohesiveZone::postInputInitialization()
{
  CoupledCohesiveZone::postInputInitialization();

  // GEOS_THROW_IF( m_characteristicNormalDisplacement < 0.0, "Characteristic normal displacement must be a positive number.", InputError );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, BicrystalCohesiveZone, std::string const &, Group * const )
}
} /* namespace geos */
