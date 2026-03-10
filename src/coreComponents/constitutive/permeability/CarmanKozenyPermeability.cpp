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
 * @file CarmanKozenyPermeability.cpp
 */

#include "CarmanKozenyPermeability.hpp"

namespace geos
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

  registerWrapper( viewKeyStruct::anisotropyString(), &m_anisotropy ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( { 1.0, 1.0, 1.0 } ).
    setDescription( "Anisotropy factors for three permeability components." );

  registerWrapper( viewKeyStruct::dPerm_dPorosityString(), &m_dPerm_dPorosity );
}

std::unique_ptr< ConstitutiveBase >
CarmanKozenyPermeability::deliverClone( string const & name,
                                        Group * const parent ) const
{
  return PermeabilityBase::deliverClone( name, parent );
}

void CarmanKozenyPermeability::allocateConstitutiveData( Group & parent,
                                                         localIndex const numPts )
{
  // NOTE: enforcing 1 quadrature point
  m_dPerm_dPorosity.resize( 0, 1, 3 );

  PermeabilityBase::allocateConstitutiveData( parent, numPts );
}

void CarmanKozenyPermeability::initializeState() const
{
  localIndex const numE = m_permeability.size( 0 );
  integer constexpr numQuad = 1; // NOTE: enforcing 1 quadrature point

  auto permView = m_permeability.toView();
  real64 const permComponents[3] = { m_particleDiameter,
                                     m_particleDiameter,
                                     m_particleDiameter };

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    for( localIndex q = 0; q < numQuad; ++q )
    {
      for( integer dim=0; dim < 3; ++dim )
      {
        // The default value is -1 so if it still -1 it needs to be set to something physical
        if( permView[ei][q][dim] < 0 )
        {
          permView[ei][q][dim] =  permComponents[dim];
        }
      }
    }
  } );
}

void CarmanKozenyPermeability::postInputInitialization()
{}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CarmanKozenyPermeability, string const &, Group * const )

}
} /* namespace geos */
