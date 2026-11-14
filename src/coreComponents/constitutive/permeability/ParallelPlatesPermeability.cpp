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
 * @file ParallelPlatesPermeability.cpp
 */

#include "ParallelPlatesPermeability.hpp"

#include "constitutive/permeability/PermeabilityFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{


ParallelPlatesPermeability::ParallelPlatesPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent ),
  m_updateTransversalComponent( true )
{
  registerWrapper( viewKeyStruct::transversalPermeabilityString(), &m_transversalPermeability ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Default value of the permeability normal to the surface. If not specified the permeability is updated using the cubic law. " );

  registerWrapper( viewKeyStruct::jointRoughnessCoefficientString(), &m_jointRoughnessCoefficient ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Joint roughness coefficient (JRC), between 1 and 20. "
                    "If specified, the aperture used in the cubic law is reduced following Barton et al. (1985): "
                    "e_h = E^2 / JRC^2.5 (apertures in microns), only where e_h <= E. "
                    "If not specified or negative (default: -1), the correction is disabled and the aperture is used directly." );

  registerField< fields::permeability::dPerm_dDispJump >( &m_dPerm_dDispJump );
}

void ParallelPlatesPermeability::allocateConstitutiveData( Group & parent,
                                                           localIndex const numPts )
{
  // NOTE: enforcing 1 quadrature point
  m_dPerm_dDispJump.resize( 0, 1, 3, 3 );

  PermeabilityBase::allocateConstitutiveData( parent, numPts );
}

void ParallelPlatesPermeability::postInputInitialization()
{
  PermeabilityBase::postInputInitialization();

  if( m_transversalPermeability > -1 )
  {
    m_updateTransversalComponent = false;
  }

  // a negative value (default: -1) means the JRC correction is not used
  if( m_jointRoughnessCoefficient >= 0 )
  {
    GEOS_THROW_IF( m_jointRoughnessCoefficient < 1 || m_jointRoughnessCoefficient > 20,
                   GEOS_FMT( "Joint Roughness Coefficient must be between 1 and 20. Provided value: {}",
                             m_jointRoughnessCoefficient ),
                   InputError, getDataContext() );
  }
}

void ParallelPlatesPermeability::initializeState() const
{
  localIndex const numE = m_permeability.size( 0 );
  integer constexpr numQuad = 1; // NOTE: enforcing 1 quadrature point

  auto permView = m_permeability.toView();

  real64 const transversalPerm = m_transversalPermeability;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    for( localIndex q = 0; q < numQuad; ++q )
    {
      permView[ei][q][2] = transversalPerm;
    }
  } );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ParallelPlatesPermeability, string const &, Group * const )

}
} /* namespace geos */
