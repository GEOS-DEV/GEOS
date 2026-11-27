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
 * @file BartonBandisPermeability.cpp
 */

#include "BartonBandisPermeability.hpp"

#include "constitutive/permeability/PermeabilityFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{


BartonBandisPermeability::BartonBandisPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent ),
  m_updateTransversalComponent( true )
{
  registerWrapper( viewKeyStruct::transversalPermeabilityString(), &m_transversalPermeability ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Default value of the permeability normal to the surface. If not specified the permeability is updated using the cubic law. " );
    
  /// TODO: must become a required parameter.
  registerWrapper( viewKeyStruct::apertureZeroString(), &m_aperture0 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1e-6 ).
    setDescription( "Reference hydraulic aperture. It is the aperture at zero normal stress." );
    
  registerWrapper( viewKeyStruct::biotCoefficientString(), &m_biotCoefficient ).
    setApplyDefaultValue( 1.0 ). 
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Biot coefficient." );
  registerWrapper( viewKeyStruct::poissonRatioString(), &m_poissonRatio ).
    setApplyDefaultValue( 0.3 ). 
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Poisson ratio." );
  registerWrapper( viewKeyStruct::normalStiffnessString(), &m_normalStiffness ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Normal stiffness: Kni." );
  registerWrapper( viewKeyStruct::referencePressureString(), &m_referencePressure ).
    setApplyDefaultValue( 1.0e5 ). 
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference pressure: p_0." );
  registerWrapper( viewKeyStruct::referenceTotalStressString(), &m_referenceTotalStress ).
    setApplyDefaultValue( { 85.0e6, 85.0e6, 105.0e6 } ). 
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Total stress at reference state: sigmaT_0." );
}


void BartonBandisPermeability::postInputInitialization()
{
  PermeabilityBase::postInputInitialization();

  if( m_transversalPermeability > -1 )
  {
    m_updateTransversalComponent = false;
  }
}

void BartonBandisPermeability::initializeState() const
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

REGISTER_CATALOG_ENTRY( ConstitutiveBase, BartonBandisPermeability, string const &, Group * const )

}
} /* namespace geos */
