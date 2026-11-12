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
    setDescription( "xx, yy and zz components of the diagonal permeability tensor at the initial state (zero volumetric strain). "
                    "Only used where the permeability is not set by a field specification." );

  registerWrapper( viewKeyStruct::permeabilityExponentsString(), &m_permeabilityExponents ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Exponents n of the porosity ratio for the xx, yy and zz permeability components: "
                    "k = k_0 ( phi / phi_0 )^n, with phi = 1 - (1 - phi_0) exp( -volumetricStrain ) (Chin et al., 2000)." );

  registerWrapper( viewKeyStruct::referencePermeabilityString(), &m_referencePermeability ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Permeability at the initial state" );

  registerWrapper( viewKeyStruct::dPerm_dVolStrainString(), &m_dPerm_dVolStrain ).
    setDescription( "Derivative of the permeability w.r.t. the volumetric strain" );

  registerWrapper( viewKeyStruct::volumetricStrainString(), &m_volStrain ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Total volumetric strain since the initial state (positive in extension)" );

  registerWrapper( viewKeyStruct::volumetricStrain_nString(), &m_volStrain_n ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Total volumetric strain at the last converged state" );

  registerWrapper( viewKeyStruct::permeability_nString(), &m_permeability_n ).
    setDescription( "Permeability at the last converged state" );
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
  // NOTE: enforcing 1 quadrature point for the permeability
  m_dPerm_dVolStrain.resize( 0, 1, 3 );
  m_referencePermeability.resize( 0, 1, 3 );
  m_permeability_n.resize( 0, 1, 3 );

  // the volumetric strain is stored at the quadrature points
  m_volStrain.resize( 0, numPts );
  m_volStrain_n.resize( 0, numPts );

  PermeabilityBase::allocateConstitutiveData( parent, numPts );
}

void StrainDependentPermeability::initializeState() const
{
  localIndex const numE = m_permeability.size( 0 );
  integer constexpr numQuad = 1; // NOTE: enforcing 1 quadrature point

  arrayView3d< real64 > const permView = m_permeability.toView();
  arrayView3d< real64 > const refPermView = m_referencePermeability.toView();
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
        // the initial permeability is the reference permeability of the model (zero volumetric strain)
        refPermView[ei][q][dim] = permView[ei][q][dim];
      }
    }
  } );

  m_permeability_n.setValues< parallelDevicePolicy<> >( m_permeability.toViewConst() );
}

void StrainDependentPermeability::saveConvergedState() const
{
  m_volStrain_n.setValues< parallelDevicePolicy<> >( m_volStrain.toViewConst() );
  m_permeability_n.setValues< parallelDevicePolicy<> >( m_permeability.toViewConst() );
}

void StrainDependentPermeability::ignoreConvergedState() const
{
  m_volStrain.setValues< parallelDevicePolicy<> >( m_volStrain_n.toViewConst() );
  m_permeability.setValues< parallelDevicePolicy<> >( m_permeability_n.toViewConst() );
}

void StrainDependentPermeability::postInputInitialization()
{
  PermeabilityBase::postInputInitialization();

  for( integer dim = 0; dim < 3; ++dim )
  {
    GEOS_THROW_IF( m_referencePermeabilityComponents[dim] <= 0.0,
                   GEOS_FMT( "{} must be strictly positive. Provided value for component {}: {}",
                             viewKeyStruct::referencePermeabilityComponentsString(), dim, m_referencePermeabilityComponents[dim] ),
                   InputError, getDataContext() );
  }
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, StrainDependentPermeability, string const &, Group * const )

}
} /* namespace geos */
