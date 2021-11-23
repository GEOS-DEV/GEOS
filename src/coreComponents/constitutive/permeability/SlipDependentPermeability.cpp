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
 * @file SlipDependentPermeability.cpp
 */

#include "SlipDependentPermeability.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


SlipDependentPermeability::SlipDependentPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::shearDispThresholdString(), &m_shearDispThreshold ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Threshold of shear displacement." );

  registerWrapper( viewKeyStruct::maxPermMultiplierString(), &m_maxPermMultiplier ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum permeability multiplier." );

  registerWrapper( viewKeyStruct::dPerm_dDispJumpString(), &m_dPerm_dDispJump ).
    setDescription( "Derivative of the permeability w.r.t. the displacement jump." );

  registerWrapper( viewKeyStruct::initialPermeabilityString(), &m_initialPermeability ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( " initial permeability of the rock." ).
    setApplyDefaultValue( -1.0 );   // will be overwritten
}

std::unique_ptr< ConstitutiveBase >
SlipDependentPermeability::deliverClone( string const & name,
                                         Group * const parent ) const
{
  return ConstitutiveBase::deliverClone( name, parent );
}

void SlipDependentPermeability::allocateConstitutiveData( dataRepository::Group & parent,
                                                          localIndex const numConstitutivePointsPerParentIndex )
{
// NOTE: enforcing 1 quadrature point
  m_initialPermeability.resize( 0, 1, 3 );
  m_dPerm_dDispJump.resize( 0, 1, 3, 3 );

  PermeabilityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}


void SlipDependentPermeability::initializeState() const
{
  arrayView3d< real64 const > permeability     = m_permeability;
  arrayView3d< real64 >       initialPermeability = m_initialPermeability;

  forAll< parallelDevicePolicy<> >( permeability.size( 0 ), [=] GEOSX_HOST_DEVICE ( localIndex const k )
  {
    LvArray::tensorOps::copy< 3 >( initialPermeability[k][0], permeability[k][0] );
  } );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, SlipDependentPermeability, string const &, Group * const )

} /* namespace constitutive */
} /* namespace geosx */
