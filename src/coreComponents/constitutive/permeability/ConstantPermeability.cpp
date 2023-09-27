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
 * @file ConstantPermeability.cpp
 */

#include "ConstantPermeability.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{


ConstantPermeability::ConstantPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent ),
  m_initialPermeability()
{
  registerWrapper( viewKeyStruct::permeabilityComponentsString(), &m_permeabilityComponents ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "xx, yy and zz components of a diagonal permeability tensor." );

  registerWrapper( viewKeyStruct::pressureDependenceConstantString(), &m_pressureDependenceConstant ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Pressure dependence constant for the permeability." );

  registerWrapper( viewKeyStruct::defaultReferencePressureString(), &m_defaultReferencePressure ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default reference pressure" );

  registerWrapper( viewKeyStruct::referencePressureString(), &m_referencePressure ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Reference pressure" );

  registerWrapper( viewKeyStruct::initialPermeabilityString(), &m_initialPermeability ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Initial permeability" );
}

std::unique_ptr< ConstitutiveBase >
ConstantPermeability::deliverClone( string const & name,
                                    Group * const parent ) const
{
  return PermeabilityBase::deliverClone( name, parent );
}

void ConstantPermeability::allocateConstitutiveData( dataRepository::Group & parent,
                                                     localIndex const numConstitutivePointsPerParentIndex )
{
  m_initialPermeability.resize( 0, 1, 3 );

  PermeabilityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  // integer const numQuad = 1; // NOTE: enforcing 1 quadrature point

  // for( localIndex ei = 0; ei < parent.size(); ++ei )
  // {
  //   for( localIndex q = 0; q < numQuad; ++q )
  //   {
  //     m_permeability[ei][q][0] =  m_permeabilityComponents[0];
  //     m_permeability[ei][q][1] =  m_permeabilityComponents[1];
  //     m_permeability[ei][q][2] =  m_permeabilityComponents[2];

  //     m_initialPermeability[ei][q][0] =  m_permeabilityComponents[0];
  //     m_initialPermeability[ei][q][1] =  m_permeabilityComponents[1];
  //     m_initialPermeability[ei][q][2] =  m_permeabilityComponents[2];
  //   }
  // }
}

void ConstantPermeability::initializeState() const
{
  localIndex const numE = m_permeability.size( 0 );
  integer constexpr numQuad = 1; // NOTE: enforcing 1 quadrature point

  auto permView = m_permeability.toView();
  auto initialPerm = m_initialPermeability.toView();
  real64 const permComponents[3] = { m_permeabilityComponents[0],
                                     m_permeabilityComponents[1],
                                     m_permeabilityComponents[2] };

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    for( localIndex q = 0; q < numQuad; ++q )
    {
      for( integer dim=0; dim < 3; ++dim )
      {
        // The default value is -1 so if it still -1 it needs to be set to something physical
        if( permView[ei][q][0] < 0 )
        {
          permView[ei][q][dim] =  permComponents[dim];
        }
      }
    }
  } );

  m_initialPermeability.setValues< parallelDevicePolicy<> >( m_permeability.toViewConst() );
}

void ConstantPermeability::postProcessInput()
{
  // set results as array default values
  this->getWrapper< array1d< real64 > >( viewKeyStruct::referencePressureString() ).
    setApplyDefaultValue( m_defaultReferencePressure );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ConstantPermeability, string const &, Group * const )

}
} /* namespace geos */
