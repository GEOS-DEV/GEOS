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
 * @file SolidBase.cpp
 */

#include "SolidBase.hpp"
#include "SolidFields.hpp"
#include "mesh/ElementSubRegionBase.hpp"

namespace geos
{
using namespace dataRepository;

namespace constitutive
{

SolidBase::SolidBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::defaultDensityString(), &m_defaultDensity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default Material Density" );

  registerWrapper( viewKeyStruct::defaultThermalExpansionCoefficientString(), &m_defaultThermalExpansionCoefficient ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Linear Thermal Expansion Coefficient of the Solid Rock Frame" );
}


void SolidBase::postInputInitialization()
{
  getWrapper< array2d< real64 > >( fields::solid::solidDensity::key() ).
    setApplyDefaultValue( m_defaultDensity );

  getWrapper< array1d< real64 > >( fields::solid::thermalExpansionCoefficient::key() ).
    setApplyDefaultValue( m_defaultThermalExpansionCoefficient );
}


void SolidBase::allocateConstitutiveData( Group & parent,
                                          localIndex const numPts )
{
  auto subregion = dynamic_cast< ElementSubRegionBase * >( &parent ); // TODO remove

  string const voightLabels[6] = { "XX", "YY", "ZZ", "YZ", "XZ", "XY" };

  subregion->registerField< fields::solid::stress >( getName(), &m_newStress ).
    setDimLabels( 2, voightLabels ).
    reference().resizeDimension< 1, 2 >( numPts, 6 );

  subregion->registerField< fields::solid::oldStress >( getName(), &m_oldStress ).
    setDimLabels( 2, voightLabels ).
    reference().resizeDimension< 1, 2 >( numPts, 6 );

  subregion->registerField< fields::solid::solidDensity >( getName(), &m_density ).
    reference().resizeDimension< 1 >( numPts );

  subregion->registerField< fields::solid::thermalExpansionCoefficient >( getName(), &m_thermalExpansionCoefficient );

  subregion->registerField< fields::solid::thermalExpansionCoefficient >( getName(), &m_thermalExpansionCoefficient );

  ConstitutiveBase::allocateConstitutiveData( parent, numPts );
}


void SolidBase::saveConvergedState() const
{
  localIndex const numE = numElem();
  localIndex const numQ = numQuad();

  arrayView3d< real64 const, solid::STRESS_USD > newStress = m_newStress;
  arrayView3d< real64, solid::STRESS_USD > oldStress = m_oldStress;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numQ; ++q )
    {
      LvArray::tensorOps::copy< 6 >( oldStress[k][q], newStress[k][q] );
    }
  } );
}


} /* namespace constitutive */
} /* namespace geos */
