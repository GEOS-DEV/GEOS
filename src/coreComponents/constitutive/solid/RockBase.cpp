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
 * @file RockBase.cpp
 */

#include "RockBase.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


RockBase::RockBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_newPorosity(),
  m_oldPorosity(),
  m_dPorosity_dPressure(),
  m_compressibility(),
  m_grainBulkModulus(),
  m_grainDensity()
{
  registerWrapper( viewKeyStruct::newPorosityString(), &m_newPorosity ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setApplyDefaultValue( -1.0 ); // will be overwritten

  registerWrapper( viewKeyStruct::oldPorosityString(), &m_oldPorosity ).
    setApplyDefaultValue( -1.0 );// will be overwritten

  registerWrapper( viewKeyStruct::dPorosity_dPressureString(), &m_dPorosity_dPressure ).
    setApplyDefaultValue( 0.0 );// will be overwritten

  registerWrapper( viewKeyStruct::compressibilityString(), &m_compressibility ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Solid compressibility" );

  registerWrapper( viewKeyStruct::grainBulkModulusString(), &m_grainBulkModulus ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Grains bulk modulus" );

  registerWrapper( viewKeyStruct::grainDensityString(), &m_grainDensity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Grains density" );
}

RockBase::~RockBase() = default;

std::unique_ptr< ConstitutiveBase >
RockBase::deliverClone( string const & name,
                        Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = ConstitutiveBase::deliverClone( name, parent );

  return clone;
}

void RockBase::allocateConstitutiveData( dataRepository::Group & parent,
                                             localIndex const numConstitutivePointsPerParentIndex )
{
  m_newPorosity.resize( 0, numConstitutivePointsPerParentIndex );
  m_oldPorosity.resize( 0, numConstitutivePointsPerParentIndex );
  m_dPorosity_dPressure.resize( 0, numConstitutivePointsPerParentIndex );

  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void RockBase::postProcessInput()
{}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, RockBase, string const &, Group * const )
}
} /* namespace geosx */
