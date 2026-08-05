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
 * @file PowerLawPermeability.cpp
 */

#include "PowerLawPermeability.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{


PowerLawPermeability::PowerLawPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::referencePermeabilityComponentsString(), &m_referencePermeabilityComponents ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Reference xx, yy and zz components of a diagonal permeability tensor, "
                    "reached at the reference porosity." );

  registerWrapper( viewKeyStruct::referencePorosityString(), &m_referencePorosity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference porosity at which the permeability is the reference permeability." );

  registerWrapper( viewKeyStruct::criticalPorosityString(), &m_criticalPorosity ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Critical porosity below which the pore space no longer percolates." );

  registerWrapper( viewKeyStruct::exponentString(), &m_exponent ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Exponent of the power law." );

  registerWrapper( viewKeyStruct::minPermeabilityString(), &m_minPermeability ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Residual permeability reached when the pore space stops percolating." );

  registerWrapper( viewKeyStruct::referencePermeabilityString(), &m_referencePermeability ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Reference permeability field" );

  registerWrapper( viewKeyStruct::dPerm_dPorosityString(), &m_dPerm_dPorosity );
}

void PowerLawPermeability::postInputInitialization()
{
  GEOS_THROW_IF_LE_MSG( m_referencePorosity, m_criticalPorosity,
                        GEOS_FMT( "{}: '{}' must be strictly larger than '{}'",
                                  getFullName(),
                                  viewKeyStruct::referencePorosityString(),
                                  viewKeyStruct::criticalPorosityString() ),
                        InputError );

  GEOS_THROW_IF_LT_MSG( m_criticalPorosity, 0.0,
                        GEOS_FMT( "{}: invalid value of attribute '{}'",
                                  getFullName(), viewKeyStruct::criticalPorosityString() ),
                        InputError );

  GEOS_THROW_IF_LT_MSG( m_exponent, 0.0,
                        GEOS_FMT( "{}: invalid value of attribute '{}'",
                                  getFullName(), viewKeyStruct::exponentString() ),
                        InputError );

  GEOS_THROW_IF_LT_MSG( m_minPermeability, 0.0,
                        GEOS_FMT( "{}: invalid value of attribute '{}'",
                                  getFullName(), viewKeyStruct::minPermeabilityString() ),
                        InputError );

  for( integer dim = 0; dim < 3; ++dim )
  {
    GEOS_THROW_IF_LT_MSG( m_referencePermeabilityComponents[dim], m_minPermeability,
                          GEOS_FMT( "{}: component {} of '{}' must be larger than '{}'",
                                    getFullName(), dim,
                                    viewKeyStruct::referencePermeabilityComponentsString(),
                                    viewKeyStruct::minPermeabilityString() ),
                          InputError );
  }
}

void PowerLawPermeability::allocateConstitutiveData( Group & parent,
                                                     localIndex const numPts )
{
  // NOTE: enforcing 1 quadrature point
  m_dPerm_dPorosity.resize( 0, 1, 3 );
  m_referencePermeability.resize( 0, 1, 3 ); // 0 to resize and assign default value later

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

void PowerLawPermeability::initializeState() const
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

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PowerLawPermeability, string const &, Group * const )

}
} /* namespace geos */
