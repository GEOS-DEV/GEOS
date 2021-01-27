
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

#include "LinearElasticAnisotropic.hpp"
#include "LinearElasticIsotropic.hpp"
#include "LinearElasticTransverseIsotropic.hpp"

namespace geosx
{

using namespace dataRepository;
namespace constitutive
{

template< typename BASE >
Damage< BASE >::Damage( std::string const & name, Group * const parent ):
  BASE( name, parent ),
  m_damage(),
  m_strainEnergyDensity(),
  m_lengthScale(),
  m_criticalFractureEnergy(),
  m_criticalStrainEnergy()
{

  this->registerWrapper( viewKeyStruct::damageString, &m_damage )->
    setApplyDefaultValue( 0.0 )->
    setPlotLevel( PlotLevel::LEVEL_0 )->
    setDescription( "Material Damage Variable" );

  this->registerWrapper( viewKeyStruct::strainEnergyDensityString, &m_strainEnergyDensity )->
    setApplyDefaultValue( 0.0 )->
    setPlotLevel( PlotLevel::LEVEL_0 )->
    setDescription( "Stress Deviator" );

  this->registerWrapper( viewKeyStruct::lengthScaleString, &m_lengthScale )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "lenght scale l in the phase-field equation" );

  this->registerWrapper( viewKeyStruct::criticalFractureEnergyString, &m_criticalFractureEnergy )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "critical fracture energy" );

  this->registerWrapper( viewKeyStruct::criticalStrainEnergyString, &m_criticalStrainEnergy )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "material critical stress in a 1d tension test" );
}


template< typename BASE >
void Damage< BASE >::postProcessInput()
{
  BASE::postProcessInput();
}

template< typename BASE >
void Damage< BASE >::allocateConstitutiveData( dataRepository::Group * const parent,
                                               localIndex const numConstitutivePointsPerParentIndex )
{
  m_damage.resize( 0, numConstitutivePointsPerParentIndex );
  m_strainEnergyDensity.resize( 0, numConstitutivePointsPerParentIndex );
  BASE::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

typedef Damage< LinearElasticIsotropic > DamageLinearElasticIsotropic;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DamageLinearElasticIsotropic, std::string const &, Group * const )

}
} /* namespace geosx */
