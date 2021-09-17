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
 * @file ProppantPermeability.cpp
 */

#include "ProppantPermeability.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


ProppantPermeability::ProppantPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent ),
  m_permeabilityMultiplier(),
  m_proppantDiameter(),
  m_maxProppantConcentration(),
  m_proppantPackPermeability()
{
  registerWrapper( viewKeyStruct::maxProppantConcentrationString(), &m_maxProppantConcentration ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Maximum proppant concentration." );

  registerWrapper( viewKeyStruct::proppantDiameterString(), &m_proppantDiameter ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Proppant diameter." );

  registerWrapper( viewKeyStruct::dPerm_dApertureString(), &m_dPerm_dAperture );
  registerWrapper( viewKeyStruct::permeabilityMultiplierString(), &m_permeabilityMultiplier );
  registerWrapper( viewKeyStruct::proppantPackPermeabilityString(), &m_proppantPackPermeability );
}

std::unique_ptr< ConstitutiveBase >
ProppantPermeability::deliverClone( string const & name,
                                    Group * const parent ) const
{
  return ConstitutiveBase::deliverClone( name, parent );
}

void ProppantPermeability::postProcessInput()
{
  real64 const oneMinusMaxConcentration = ( 1.0 - m_maxProppantConcentration );
  m_proppantPackPermeability  = m_proppantDiameter * m_proppantDiameter / 180.0;
  m_proppantPackPermeability *= ( oneMinusMaxConcentration * oneMinusMaxConcentration * oneMinusMaxConcentration )
                                / ( m_maxProppantConcentration * m_maxProppantConcentration );
}

void ProppantPermeability::allocateConstitutiveData( dataRepository::Group & parent,
                                                     localIndex const numConstitutivePointsPerParentIndex )
{
  // NOTE: enforcing 1 quadrature point
  m_dPerm_dAperture.resize( 0, 1, 3 );
  m_permeabilityMultiplier.resize( 0, 1, 3 );
  PermeabilityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, ProppantPermeability, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geosx */
