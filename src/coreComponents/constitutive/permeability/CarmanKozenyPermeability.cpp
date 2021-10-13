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
 * @file CarmanKozenyPermeability.cpp
 */

#include "CarmanKozenyPermeability.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


CarmanKozenyPermeability::CarmanKozenyPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::particleDiameterString(), &m_particleDiameter ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Diameter of the spherical particles." );

  registerWrapper( viewKeyStruct::sphericityString(), &m_sphericity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Sphericity of the particles." );

  registerWrapper( viewKeyStruct::dPerm_dPorosityString(), &m_dPerm_dPorosity );
}

std::unique_ptr< ConstitutiveBase >
CarmanKozenyPermeability::deliverClone( string const & name,
                                        Group * const parent ) const
{
  return PermeabilityBase::deliverClone( name, parent );
}

void CarmanKozenyPermeability::allocateConstitutiveData( dataRepository::Group & parent,
                                                         localIndex const numConstitutivePointsPerParentIndex )
{
  // NOTE: enforcing 1 quadrature point
  m_dPerm_dPorosity.resize( 0, 1, 3 );
  PermeabilityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CarmanKozenyPermeability, string const &, Group * const )

}
} /* namespace geosx */
