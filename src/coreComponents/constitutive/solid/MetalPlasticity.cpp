/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file MetalPlasticity.cpp
 */

#include "MetalPlasticity.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

MetalPlasticity::MetalPlasticity( string const & name, Group * const parent ):
  ElasticIsotropic( name, parent ),
  m_defaultYieldStrength(),
  m_yieldStrength(),
  m_deformationGradient(),
  m_velocityGradient(),
  m_plasticStrain(),
  m_temperature()
{
  registerWrapper( viewKeyStruct::defaultYieldStrengthString(), &m_defaultYieldStrength ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default yield strength parameter" );

  registerWrapper( viewKeyStruct::yieldStrengthString(), &m_yieldStrength ).
    setApplyDefaultValue( -1 ).
    setDescription( "Yield strength field" );

  registerWrapper( viewKeyStruct::deformationGradientString(), &m_deformationGradient ).
    setApplyDefaultValue( 1.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Deformation gradient field" );

  registerWrapper( viewKeyStruct::velocityGradientString(), &m_velocityGradient ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Velocity gradient field" );

  registerWrapper( viewKeyStruct::plasticStrainString(), &m_plasticStrain ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Plastic strain field" );

  registerWrapper( viewKeyStruct::temperatureString(), &m_temperature ).
    setApplyDefaultValue( 300.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Array of quadrature point temperature values" );
}

MetalPlasticity::~MetalPlasticity()
{}

void MetalPlasticity::allocateConstitutiveData( dataRepository::Group & parent,
                                          localIndex const numConstitutivePointsPerParentIndex )
{
  ElasticIsotropic::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_deformationGradient.resize( 0, 3, 3 );
  m_velocityGradient.resize( 0, 3, 3 );
  m_plasticStrain.resize( 0, numConstitutivePointsPerParentIndex, 6 );
  m_yieldStrength.resize( 0 );
}


void MetalPlasticity::postInputInitialization()
{
  ElasticIsotropic::postInputInitialization();

  this->getWrapper< array1d< real64 > >( viewKeyStruct::yieldStrengthString() ).
    setApplyDefaultValue( m_defaultYieldStrength );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, MetalPlasticity, string const &, Group * const )
}
} /* namespace geos */
