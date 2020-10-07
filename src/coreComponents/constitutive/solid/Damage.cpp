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
 * @file Damage.cpp
 */

#include "Damage.hpp"
#include "ElasticIsotropic.hpp"
#include "ElasticTransverseIsotropic.hpp"

namespace geosx
{

using namespace dataRepository;
namespace constitutive
{

template< typename BASE >
Damage< BASE >::Damage( string const & name, Group * const parent ):
  BASE( name, parent ),
  m_damage(),
  m_strainEnergyDensity()
{

  this->registerWrapper( viewKeyStruct::damageString, &m_damage )->
    setApplyDefaultValue( 0.0 )->
    setPlotLevel( PlotLevel::LEVEL_0 )->
    setDescription( "Material Damage Variable" );

  this->registerWrapper( viewKeyStruct::strainEnergyDensityString, &m_strainEnergyDensity )->
    setApplyDefaultValue( 0.0 )->
    setPlotLevel( PlotLevel::LEVEL_0 )->
    setDescription( "Stress Deviator" );

}

template< typename BASE >
Damage< BASE >::~Damage()
{}

template< typename BASE >
void Damage< BASE >::PostProcessInput()
{
  BASE::PostProcessInput();
}

template< typename BASE >
void Damage< BASE >::allocateConstitutiveData( dataRepository::Group * const parent,
                                               localIndex const numConstitutivePointsPerParentIndex )
{
  m_damage.resize( 0, numConstitutivePointsPerParentIndex );
  m_strainEnergyDensity.resize( 0, numConstitutivePointsPerParentIndex );
  BASE::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

typedef Damage< ElasticIsotropic > DamageElasticIsotropic;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DamageElasticIsotropic, string const &, Group * const )

}
} /* namespace geosx */
