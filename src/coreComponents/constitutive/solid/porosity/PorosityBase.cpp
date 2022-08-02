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

#include "constitutive/solid/porosity/PorosityBase.hpp"
#include "constitutive/solid/porosity/PorosityExtrinsicData.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


PorosityBase::PorosityBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_newPorosity(),
  m_porosity_n(),
  m_dPorosity_dPressure(),
  m_initialPorosity(),
  m_referencePorosity(),
  m_defaultReferencePorosity()
{
  registerWrapper( viewKeyStruct::defaultReferencePorosityString(), &m_defaultReferencePorosity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default value of the reference porosity" );

  registerExtrinsicData( extrinsicMeshData::porosity::porosity{}, &m_newPorosity );

  registerExtrinsicData( extrinsicMeshData::porosity::porosity_n{}, &m_porosity_n );

  registerExtrinsicData( extrinsicMeshData::porosity::dPorosity_dPressure{}, &m_dPorosity_dPressure );

  registerExtrinsicData( extrinsicMeshData::porosity::initialPorosity{}, &m_initialPorosity );

  registerExtrinsicData( extrinsicMeshData::porosity::referencePorosity{}, &m_referencePorosity );
}

void PorosityBase::allocateConstitutiveData( dataRepository::Group & parent,
                                             localIndex const numConstitutivePointsPerParentIndex )
{
  m_newPorosity.resize( 0, numConstitutivePointsPerParentIndex );
  m_porosity_n.resize( 0, numConstitutivePointsPerParentIndex );
  m_dPorosity_dPressure.resize( 0, numConstitutivePointsPerParentIndex );
  m_initialPorosity.resize( 0, numConstitutivePointsPerParentIndex );

  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void PorosityBase::postProcessInput()
{
  getExtrinsicData< extrinsicMeshData::porosity::referencePorosity >().
    setApplyDefaultValue( m_defaultReferencePorosity );
}


void PorosityBase::saveConvergedState() const
{
  localIndex const numE = numElem();
  localIndex const numQ = numQuad();

  arrayView2d< real64 const > newPorosity = m_newPorosity;
  arrayView2d< real64 >       porosity_n  = m_porosity_n;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOSX_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numQ; ++q )
    {
      porosity_n[k][q] = newPorosity[k][q];
    }
  } );
}

void PorosityBase::initializeState() const
{
  localIndex const numE = numElem();
  localIndex const numQ = numQuad();

  arrayView2d< real64 const > newPorosity     = m_newPorosity;
  arrayView2d< real64 >       porosity_n      = m_porosity_n;
  arrayView2d< real64 >       initialPorosity = m_initialPorosity;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOSX_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numQ; ++q )
    {
      porosity_n[k][q]      = newPorosity[k][q];
      initialPorosity[k][q] = newPorosity[k][q];
    }
  } );
}

}
} /* namespace geosx */
