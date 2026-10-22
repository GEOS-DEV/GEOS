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

  // register fields

  string const voightLabels[6] = { "XX", "YY", "ZZ", "YZ", "XZ", "XY" };

  registerField< fields::solid::stress >( &m_newStress ).
    setDimLabels( 2, voightLabels );

  registerField< fields::solid::oldStress >( &m_oldStress ).
    setDimLabels( 2, voightLabels );

  registerField< fields::solid::density >( &m_density );

  registerField< fields::solid::thermalExpansionCoefficient >( &m_thermalExpansionCoefficient );
}


void SolidBase::postInputInitialization()
{
  getField< fields::solid::density >().
    setApplyDefaultValue( m_defaultDensity );

  getField< fields::solid::thermalExpansionCoefficient >().
    setApplyDefaultValue( m_defaultThermalExpansionCoefficient );
}


void SolidBase::allocateConstitutiveData( Group & parent, localIndex const numPts )
{
  // 0 to resize and assign default value later
  m_density.resize( 0, numPts );
  m_newStress.resize( 0, numPts, 6 );
  m_oldStress.resize( 0, numPts, 6 );

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
