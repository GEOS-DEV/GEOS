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

  registerWrapper( viewKeyStruct::defaultAnelasticStrainRateString(), &m_defaultAnelasticStrainRate ).
    setApplyDefaultValue( R1Tensor{} ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default anelastic strain rate components" );

  registerWrapper( viewKeyStruct::enableAnelasticStrainString(), &m_enableAnelasticStrain ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Whether to enable stress modification due to anelastic strain. Can be 0 or 1" );

  // register fields

  string const voightLabels[6] = { "XX", "YY", "ZZ", "YZ", "XZ", "XY" };
  string const normalStrainLabels[3] = { "XX", "YY", "ZZ" };

  registerField< fields::solid::stress >( &m_newStress ).
    setDimLabels( 2, voightLabels );

  registerField< fields::solid::oldStress >( &m_oldStress ).
    setDimLabels( 2, voightLabels );

  registerField< fields::solid::density >( &m_density );

  registerField< fields::solid::thermalExpansionCoefficient >( &m_thermalExpansionCoefficient );

  registerField< fields::solid::anelasticStrainRate >( &m_anelasticStrainRate ).
    setDimLabels( 1, normalStrainLabels );
  registerField< fields::solid::newAnelasticStrain >( &m_newAnelasticStrain ).
    setDimLabels( 1, normalStrainLabels );
  registerField< fields::solid::oldAnelasticStrain >( &m_oldAnelasticStrain ).
    setDimLabels( 1, normalStrainLabels );
}


void SolidBase::postInputInitialization()
{
  getField< fields::solid::density >().
    setApplyDefaultValue( m_defaultDensity );

  getField< fields::solid::thermalExpansionCoefficient >().
    setApplyDefaultValue( m_defaultThermalExpansionCoefficient );

  getField< fields::solid::anelasticStrainRate >().
    setApplyDefaultValue( 0.0 );

  GEOS_ERROR_IF( m_enableAnelasticStrain == 0 &&
                 LvArray::tensorOps::l2NormSquared< 3 >( m_defaultAnelasticStrainRate ) > 0.0,
                 getDataContext() << ": enableAnelasticStrain flag must be 1 if a nonzero"
                                     " AnelasticStrainRate is used" );
}


void SolidBase::allocateConstitutiveData( Group & parent, localIndex const numPts )
{
  // 0 to resize and assign default value later
  m_density.resize( 0, numPts );
  m_newStress.resize( 0, numPts, 6 );
  m_oldStress.resize( 0, numPts, 6 );
  m_anelasticStrainRate.resize( 0, 3 );
  m_newAnelasticStrain.resize( 0, 3 );
  m_oldAnelasticStrain.resize( 0, 3 );

  ConstitutiveBase::allocateConstitutiveData( parent, numPts );

  for( localIndex ei = 0; ei < parent.size(); ++ei )
  {
    for( integer dim = 0; dim < 3; ++dim )
    {
      m_anelasticStrainRate[ei][dim] = m_defaultAnelasticStrainRate[dim];
    }
  }
}


void SolidBase::saveConvergedState() const
{
  localIndex const numE = numElem();
  localIndex const numQ = numQuad();

  arrayView3d< real64 const, solid::STRESS_USD > newStress = m_newStress;
  arrayView3d< real64, solid::STRESS_USD > oldStress = m_oldStress;

  arrayView2d< real64 const > newAnelasticStrain = m_newAnelasticStrain;
  arrayView2d< real64 > oldAnelasticStrain = m_oldAnelasticStrain;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    LvArray::tensorOps::copy< 3 >( oldAnelasticStrain[k], newAnelasticStrain[k] );

    for( localIndex q = 0; q < numQ; ++q )
    {
      LvArray::tensorOps::copy< 6 >( oldStress[k][q], newStress[k][q] );
    }
  } );
}


} /* namespace constitutive */
} /* namespace geos */
