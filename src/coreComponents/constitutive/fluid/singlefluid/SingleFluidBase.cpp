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
 * @file SingleFluidBase.cpp
 */

#include "SingleFluidBase.hpp"

#include "SingleFluidFields.hpp"
#include "mesh/ElementSubRegionBase.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

SingleFluidBase::SingleFluidBase( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent ),
  m_numDOF( 1 )
{}

void SingleFluidBase::postInputInitialization()
{
  ConstitutiveBase::postInputInitialization();

  // for fracture elements, set the default value - TODO check why this is needed
  getField< fields::singlefluid::density_n >().
    setDefaultValue( defaultDensity() );
}

void SingleFluidBase::saveConvergedState() const
{
  localIndex const numElem = m_density.value.size( 0 );
  localIndex const numGauss = m_density.value.size( 1 );

  SingleFluidProp::ViewTypeConst const density = m_density.toViewConst();
  SingleFluidProp::ViewTypeConst const internalEnergy = m_internalEnergy.toViewConst();

  arrayView2d< real64, singlefluid::USD_FLUID > const density_n = m_density_n.toView();
  arrayView2d< real64, singlefluid::USD_FLUID > const internalEnergy_n = m_internalEnergy_n.toView();

  forAll< parallelDevicePolicy<> >( numElem, [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numGauss; ++q )
    {
      density_n[k][q] = density.value[k][q];
      internalEnergy_n[k][q] = internalEnergy.value[k][q];
    }
  } );
}

void SingleFluidBase::allocateConstitutiveData( dataRepository::Group & parent,
                                                localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  auto subRegion = dynamic_cast< ElementSubRegionBase * >( &parent );

  // density
  subRegion->registerField< fields::singlefluid::density >( &m_density.value ).
    reference().resizeDimension< 1 >( numConstitutivePointsPerParentIndex );
  subRegion->registerField< fields::singlefluid::dDensity >( &m_density.derivs ).
    reference().resizeDimension< 1, 2 >( numConstitutivePointsPerParentIndex, m_numDOF );
  subRegion->registerField< fields::singlefluid::density_n >( &m_density_n ).
    reference().resizeDimension< 1 >( numConstitutivePointsPerParentIndex );

  // viscosity
  subRegion->registerField< fields::singlefluid::viscosity >( &m_viscosity.value ).
    reference().resizeDimension< 1 >( numConstitutivePointsPerParentIndex );
  subRegion->registerField< fields::singlefluid::dViscosity >( &m_viscosity.derivs ).
    reference().resizeDimension< 1, 2 >( numConstitutivePointsPerParentIndex, m_numDOF );

  // internal energy
  subRegion->registerField< fields::singlefluid::internalEnergy >( &m_internalEnergy.value ).
    reference().resizeDimension< 1 >( numConstitutivePointsPerParentIndex );
  subRegion->registerField< fields::singlefluid::dInternalEnergy >( &m_internalEnergy.derivs ).
    reference().resizeDimension< 1, 2 >( numConstitutivePointsPerParentIndex, m_numDOF );
  subRegion->registerField< fields::singlefluid::internalEnergy_n >( &m_internalEnergy_n ).
    reference().resizeDimension< 1 >( numConstitutivePointsPerParentIndex );

  // enthalpy
  subRegion->registerField< fields::singlefluid::enthalpy >( &m_enthalpy.value ).
    reference().resizeDimension< 1 >( numConstitutivePointsPerParentIndex );
  subRegion->registerField< fields::singlefluid::dEnthalpy >( &m_enthalpy.derivs ).
    reference().resizeDimension< 1, 2 >( numConstitutivePointsPerParentIndex, m_numDOF );
}

} //namespace constitutive

} //namespace geos
