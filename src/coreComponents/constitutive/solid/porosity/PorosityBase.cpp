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
 * @file PorosityBase.cpp
 */

#include "constitutive/solid/porosity/PorosityBase.hpp"
#include "constitutive/solid/porosity/PorosityFields.hpp"
#include "mesh/ElementSubRegionBase.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{


PorosityBase::PorosityBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::defaultReferencePorosityString(), &m_defaultReferencePorosity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default value of the reference porosity" );
}

void PorosityBase::allocateConstitutiveData( Group & parent,
                                             localIndex const numPts )
{
  auto subregion = dynamic_cast< ElementSubRegionBase * >( &parent ); // TODO remove

  subregion->registerField< fields::porosity::porosity >( &m_newPorosity ).
    reference().resizeDimension< 1 >( numPts );

  subregion->registerField< fields::porosity::porosity_n >( &m_porosity_n ).
    reference().resizeDimension< 1 >( numPts );

  subregion->registerField< fields::porosity::dPorosity_dPressure >( &m_dPorosity_dPressure ).
    reference().resizeDimension< 1 >( numPts );

  subregion->registerField< fields::porosity::dPorosity_dTemperature >( &m_dPorosity_dTemperature ).
    reference().resizeDimension< 1 >( numPts );

  subregion->registerField< fields::porosity::initialPorosity >( &m_initialPorosity ).
    reference().resizeDimension< 1 >( numPts );

  subregion->registerField< fields::porosity::referencePorosity >( &m_referencePorosity ).
    setApplyDefaultValue( m_defaultReferencePorosity );

  ConstitutiveBase::allocateConstitutiveData( parent, numPts );
}

void PorosityBase::postInputInitialization()
{
//  getField< fields::porosity::referencePorosity >().
//    setApplyDefaultValue( m_defaultReferencePorosity );
}

void PorosityBase::scaleReferencePorosity( arrayView1d< real64 const > scalingFactors ) const
{
  localIndex const numE = numElem();

  arrayView1d< real64 > referencePorosity = m_referencePorosity;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    referencePorosity[k] *= scalingFactors[k];
  } );
}

void PorosityBase::saveConvergedState() const
{
  m_porosity_n.setValues< parallelDevicePolicy<> >( m_newPorosity.toViewConst() );
}

void PorosityBase::ignoreConvergedState() const
{
  m_newPorosity.setValues< parallelDevicePolicy<> >( m_porosity_n.toViewConst() );
}

void PorosityBase::initializeState() const
{
  m_porosity_n.setValues< parallelDevicePolicy<> >( m_newPorosity.toViewConst() );
  m_initialPorosity.setValues< parallelDevicePolicy<> >( m_newPorosity.toViewConst() );
}

}
} /* namespace geos */
