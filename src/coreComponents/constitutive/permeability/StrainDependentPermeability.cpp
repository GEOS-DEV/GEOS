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
 * @file StrainDependentPermeability.cpp
 */

#include "StrainDependentPermeability.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{


StrainDependentPermeability::StrainDependentPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::referencePermeabilityComponentsString(), &m_referencePermeabilityComponents ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Reference xx, yy and zz components of a diagonal permeability tensor." );

  registerWrapper( viewKeyStruct::strainDependenceConstantsString(), &m_strainDependenceConstants ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Volumetric strain dependence coefficients for each permeability component." );

  registerWrapper( viewKeyStruct::referencePermeabilityString(), &m_referencePermeability ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Reference permeability field" );

  registerWrapper( viewKeyStruct::dPerm_dPorosityString(), &m_dPerm_dPorosity );
}

std::unique_ptr< ConstitutiveBase >
StrainDependentPermeability::deliverClone( string const & name,
                                           Group * const parent ) const
{
  return PermeabilityBase::deliverClone( name, parent );
}

void StrainDependentPermeability::allocateConstitutiveData( Group & parent,
                                                            localIndex const numPts )
{
  // NOTE: enforcing 1 quadrature point
  m_dPerm_dPorosity.resize( 0, 1, 3 );

  m_referencePermeability.resize( 0, 1, 3 );

  PermeabilityBase::allocateConstitutiveData( parent, numPts );

  integer constexpr numQuad = 1; // NOTE: enforcing 1 quadrature point

  for( localIndex ei = 0; ei < parent.size(); ++ei )
  {
    for( localIndex q = 0; q < numQuad; ++q )
    {
      m_referencePermeability[ei][q][0] = m_referencePermeabilityComponents[0];
      m_referencePermeability[ei][q][1] = m_referencePermeabilityComponents[1];
      m_referencePermeability[ei][q][2] = m_referencePermeabilityComponents[2];
    }
  }
}

void StrainDependentPermeability::initializeState() const
{
  localIndex const numE = m_permeability.size( 0 );
  integer constexpr numQuad = 1; // NOTE: enforcing 1 quadrature point

  auto permView = m_permeability.toView();
  real64 const permComponents[3] = { m_referencePermeabilityComponents[0],
                                     m_referencePermeabilityComponents[1],
                                     m_referencePermeabilityComponents[2] };

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    for( localIndex q = 0; q < numQuad; ++q )
    {
      for( integer dim=0; dim < 3; ++dim )
      {
        // The default value is -1 so if it still -1 it needs to be set to something physical
        if( permView[ei][q][dim] < 0 )
        {
          permView[ei][q][dim] = permComponents[dim];
        }
      }
    }
  } );
}

void StrainDependentPermeability::postInputInitialization()
{}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, StrainDependentPermeability, string const &, Group * const )

}
} /* namespace geos */
