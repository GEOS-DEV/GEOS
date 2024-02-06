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
 * @file SolidBase.cpp
 */

#include "SolidBase.hpp"

namespace geos
{
using namespace dataRepository;

namespace constitutive
{

SolidBase::SolidBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_newStress( 0, 0, 6 ),
  m_oldStress( 0, 0, 6 ),
  m_density(),
  m_thermalExpansionCoefficient()
{
  string const voightLabels[6] = { "XX", "YY", "ZZ", "YZ", "XZ", "XY" };

  registerWrapper( viewKeyStruct::stressString(), &m_newStress ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setApplyDefaultValue( 0 ). // default to zero initial stress
    setDescription( "Current Material Stress" ).
    setDimLabels( 2, voightLabels );

  registerWrapper( viewKeyStruct::oldStressString(), &m_oldStress ).
    setApplyDefaultValue( 0 ). // default to zero initial stress
    setDescription( "Previous Material Stress" );

  registerWrapper( viewKeyStruct::densityString(), &m_density ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setApplyDefaultValue( -1 ). // will be overwritten
    setDescription( "Material Density" );

  registerWrapper( viewKeyStruct::defaultDensityString(), &m_defaultDensity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default Material Density" );

#if !defined(GEOS_TEMP_MINIMUM_ALLOCATION_FLAG)
  registerWrapper( viewKeyStruct::defaultThermalExpansionCoefficientString(), &m_defaultThermalExpansionCoefficient ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Linear Thermal Expansion Coefficient of the Solid Rock Frame" );

  registerWrapper( viewKeyStruct::thermalExpansionCoefficientString(), &m_thermalExpansionCoefficient ).
    setApplyDefaultValue( -1.0 ). // will be overwritten
    setDescription( "Linear Thermal Expansion Coefficient Field" );
#endif
}


SolidBase::~SolidBase()
{}


void SolidBase::postProcessInput()
{
  this->getWrapper< array2d< real64 > >( viewKeyStruct::densityString() ).
    setApplyDefaultValue( m_defaultDensity );

#if !defined(GEOS_TEMP_MINIMUM_ALLOCATION_FLAG)
  this->getWrapper< array1d< real64 > >( viewKeyStruct::thermalExpansionCoefficientString() ).
    setApplyDefaultValue( m_defaultThermalExpansionCoefficient );
#endif
}


void SolidBase::allocateConstitutiveData( dataRepository::Group & parent,
                                          localIndex const numConstitutivePointsPerParentIndex )
{
#if !defined(GEOS_TEMP_MINIMUM_ALLOCATION_FLAG)
  m_density.resize( 0, numConstitutivePointsPerParentIndex );
  m_newStress.resize( 0, numConstitutivePointsPerParentIndex, 6 );
  m_oldStress.resize( 0, numConstitutivePointsPerParentIndex, 6 );
#endif
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
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
