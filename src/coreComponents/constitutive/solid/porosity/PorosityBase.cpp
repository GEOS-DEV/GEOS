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
 * @file PorosityBase.cpp
 */

#include "PorosityBase.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


PorosityBase::PorosityBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_newPorosity(),
  m_oldPorosity(),
  m_dPorosity_dPressure(),
  m_initialPorosity(),
  m_referencePorosity(),
  m_defaultReferencePorosity()
{
  registerWrapper( viewKeyStruct::newPorosityString(), &m_newPorosity ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setApplyDefaultValue( 1.0 ); // will be overwritten but it's important for newly created faceElements.

  registerWrapper( viewKeyStruct::oldPorosityString(), &m_oldPorosity ).
    setApplyDefaultValue( 1.0 ); // will be overwritten but it's important for newly created faceElements.

  registerWrapper( viewKeyStruct::dPorosity_dPressureString(), &m_dPorosity_dPressure ).
    setApplyDefaultValue( 0.0 ); // will be overwritten

  registerWrapper( viewKeyStruct::initialPorosityString(), &m_initialPorosity ).
    setApplyDefaultValue( 0.0 ); // will be overwritten

  registerWrapper( viewKeyStruct::defaultRefererencePorosityString(), &m_defaultReferencePorosity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default value of the reference porosity" );

  registerWrapper( viewKeyStruct::referencePorosityString(), &m_referencePorosity ).
    setApplyDefaultValue( 1.0 );
}

void PorosityBase::allocateConstitutiveData( dataRepository::Group & parent,
                                             localIndex const numConstitutivePointsPerParentIndex )
{
  m_newPorosity.resize( 0, numConstitutivePointsPerParentIndex );
  m_oldPorosity.resize( 0, numConstitutivePointsPerParentIndex );
  m_dPorosity_dPressure.resize( 0, numConstitutivePointsPerParentIndex );
  m_initialPorosity.resize( 0, numConstitutivePointsPerParentIndex );

  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void PorosityBase::postProcessInput()
{
  getWrapper< array1d< real64 > >( viewKeyStruct::referencePorosityString() )
    .setApplyDefaultValue( m_defaultReferencePorosity );
}


void PorosityBase::saveConvergedState() const
{
  localIndex const numE = numElem();
  localIndex const numQ = numQuad();

  arrayView2d< real64 const > newPorosity = m_newPorosity;
  arrayView2d< real64 >       oldPorosity = m_oldPorosity;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOSX_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numQ; ++q )
    {
      oldPorosity[k][q] = newPorosity[k][q];
    }
  } );
}

void PorosityBase::initializeState() const
{
  localIndex const numE = numElem();
  localIndex const numQ = numQuad();

  arrayView2d< real64 const > newPorosity     = m_newPorosity;
  arrayView2d< real64 >       oldPorosity     = m_oldPorosity;
  arrayView2d< real64 >       initialPorosity = m_initialPorosity;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOSX_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numQ; ++q )
    {
      oldPorosity[k][q]     = newPorosity[k][q];
      initialPorosity[k][q] = newPorosity[k][q];
    }
  } );
}

}
} /* namespace geosx */
