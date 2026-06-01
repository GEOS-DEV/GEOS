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
 * @file SlurryFluidBase.cpp
 */

#include "SlurryFluidBase.hpp"

#include "SlurryFluidFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

SlurryFluidBase::SlurryFluidBase( string const & name, Group * const parent ):
  SingleFluidBase( name, parent ),
  m_isNewtonianFluid( true )
{

  registerWrapper( viewKeyStruct::componentNamesString(), &m_componentNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "List of fluid component names" );

  registerWrapper( viewKeyStruct::defaultComponentDensityString(), &m_defaultComponentDensity ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Default value for the component density." );

  registerWrapper( viewKeyStruct::defaultCompressibilityString(), &m_defaultComponentCompressibility ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Default value for the component compressibility." );

  registerWrapper( viewKeyStruct::defaultComponentViscosityString(), &m_defaultComponentViscosity ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default value for the component viscosity." );

  registerWrapper( viewKeyStruct::flowBehaviorIndexString(), &m_nIndices ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flow behavior index" );

  registerWrapper( viewKeyStruct::flowConsistencyIndexString(), &m_Ks ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flow consistency index" );

  // these would be in dDensity
  registerField< fields::slurryfluid::dDensity_dProppantConcentration >( &m_dDensity_dProppantConc );
  registerField< fields::slurryfluid::dDensity_dComponentConcentration >( &m_dDensity_dCompConc );


  registerField< fields::slurryfluid::fluidDensity >( &m_fluidDensity.value );
  registerField< fields::slurryfluid::dFluidDensity_dPressure >( &m_dFluidDens_dPres );
  registerField< fields::slurryfluid::dFluidDensity_dComponentConcentration >( &m_dFluidDens_dCompConc );

  registerField< fields::slurryfluid::fluidViscosity >( &m_fluidViscosity );
  registerField< fields::slurryfluid::dFluidViscosity_dPressure >( &m_dFluidVisc_dPres );
  registerField< fields::slurryfluid::dFluidViscosity_dComponentConcentration >( &m_dFluidVisc_dCompConc );

  registerField< fields::slurryfluid::componentDensity >( &m_componentDensity );
  registerField< fields::slurryfluid::dComponentDensity_dPressure >( &m_dCompDens_dPres );
  registerField< fields::slurryfluid::dComponentDensity_dComponentConcentration >( &m_dCompDens_dCompConc );

  registerField< fields::slurryfluid::dViscosity_dProppantConcentration >( &m_dViscosity_dProppantConc );
  registerField< fields::slurryfluid::dViscosity_dComponentConcentration >( &m_dViscosity_dCompConc );
}

void SlurryFluidBase::postInputInitialization()
{
  SingleFluidBase::postInputInitialization();

  localIndex const NC = numFluidComponents();

  GEOS_ERROR_IF( m_defaultComponentDensity.size() != NC,
                 "The number of default density values is not the same as the component number",
                 getDataContext() );

  GEOS_ERROR_IF( m_defaultComponentCompressibility.size() != NC,
                 "The number of default compressibility values is not the same as the component number",
                 getDataContext() );

  GEOS_ERROR_IF( m_defaultComponentViscosity.size() != NC,
                 "The number of default viscosity values is not the same as the component number",
                 getDataContext() );

}

localIndex SlurryFluidBase::numFluidComponents() const
{
  return LvArray::integerConversion< localIndex >( m_componentNames.size());
}

void SlurryFluidBase::allocateConstitutiveData( Group & parent,
                                                localIndex const numPts )
{
  localIndex const NC = numFluidComponents();
  m_numDOF = 2 + NC;  // pressure,proppantconc, NC compconc

  // These are also sized in m_dDenisty in base class , only dP and dT are populated
  // Future dev should incorporate concentration derivatives in dDensity
  m_dDensity_dProppantConc.resize( 0, numPts );
  m_dDensity_dCompConc.resize( 0, numPts, NC );

  m_componentDensity.resize( 0, numPts, NC );
  m_dCompDens_dPres.resize( 0, numPts, NC );
  m_dCompDens_dCompConc.resize( 0, numPts, NC, NC );

  m_fluidDensity.value.resize( 0, numPts );
  m_dFluidDens_dPres.resize( 0, numPts );
  m_dFluidDens_dCompConc.resize( 0, numPts, NC );

  m_fluidViscosity.resize( 0, numPts );
  m_dFluidVisc_dPres.resize( 0, numPts );
  m_dFluidVisc_dCompConc.resize( 0, numPts, NC );

  m_dViscosity_dProppantConc.resize( 0, numPts );
  m_dViscosity_dCompConc.resize( 0, numPts, NC );

  SingleFluidBase::allocateConstitutiveData( parent, numPts );
}


} //namespace constitutive

} //namespace geos
