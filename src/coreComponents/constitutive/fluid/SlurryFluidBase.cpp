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
 * @file SlurryFluidBase.cpp
 */

#include "SlurryFluidBase.hpp"

#include "SlurryFluidExtrinsicData.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

SlurryFluidBase::SlurryFluidBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_isNewtonianFluid( true )
{

  registerWrapper( viewKeyStruct::componentNamesString(), &m_componentNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of fluid component names" );

  registerWrapper( viewKeyStruct::defaultDensityString(), &m_defaultDensity ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default value for density." );

  registerWrapper( viewKeyStruct::defaultCompressibilityString(), &m_defaultCompressibility ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default value for compressibility." );

  registerWrapper( viewKeyStruct::defaultViscosityString(), &m_defaultViscosity ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default value for viscosity." );

  registerWrapper( viewKeyStruct::flowBehaviorIndexString(), &m_nIndices ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flow behavior index" );

  registerWrapper( viewKeyStruct::flowConsistencyIndexString(), &m_Ks ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flow consistency index" );

  registerExtrinsicData( extrinsicMeshData::slurryfluid::density{}, &m_density );
  registerExtrinsicData( extrinsicMeshData::slurryfluid::dDensity_dPressure{}, &m_dDens_dPres );
  registerExtrinsicData( extrinsicMeshData::slurryfluid::dDensity_dProppantConcentration{}, &m_dDens_dProppantConc );
  registerExtrinsicData( extrinsicMeshData::slurryfluid::dDensity_dComponentConcentration{}, &m_dDens_dCompConc );

  registerExtrinsicData( extrinsicMeshData::slurryfluid::fluidDensity{}, &m_fluidDensity );
  registerExtrinsicData( extrinsicMeshData::slurryfluid::dFluidDensity_dPressure{}, &m_dFluidDens_dPres );
  registerExtrinsicData( extrinsicMeshData::slurryfluid::dFluidDensity_dComponentConcentration{}, &m_dFluidDens_dCompConc );

  registerExtrinsicData( extrinsicMeshData::slurryfluid::fluidViscosity{}, &m_fluidViscosity );
  registerExtrinsicData( extrinsicMeshData::slurryfluid::dFluidViscosity_dPressure{}, &m_dFluidVisc_dPres );
  registerExtrinsicData( extrinsicMeshData::slurryfluid::dFluidViscosity_dComponentConcentration{}, &m_dFluidVisc_dCompConc );

  registerExtrinsicData( extrinsicMeshData::slurryfluid::componentDensity{}, &m_componentDensity );
  registerExtrinsicData( extrinsicMeshData::slurryfluid::dComponentDensity_dPressure{}, &m_dCompDens_dPres );
  registerExtrinsicData( extrinsicMeshData::slurryfluid::dComponentDensity_dComponentConcentration{}, &m_dCompDens_dCompConc );

  registerExtrinsicData( extrinsicMeshData::slurryfluid::viscosity{}, &m_viscosity );
  registerExtrinsicData( extrinsicMeshData::slurryfluid::dViscosity_dPressure{}, &m_dVisc_dPres );
  registerExtrinsicData( extrinsicMeshData::slurryfluid::dViscosity_dProppantConcentration{}, &m_dVisc_dProppantConc );
  registerExtrinsicData( extrinsicMeshData::slurryfluid::dViscosity_dComponentConcentration{}, &m_dVisc_dCompConc );
}

SlurryFluidBase::~SlurryFluidBase() = default;

void SlurryFluidBase::postProcessInput()
{
  ConstitutiveBase::postProcessInput();

  localIndex const NC = numFluidComponents();

  GEOSX_ERROR_IF( m_defaultDensity.size() != NC,
                  "The number of default density values is not the same as the component number" );

  GEOSX_ERROR_IF( m_defaultCompressibility.size() != NC,
                  "The number of default compressibility values is not the same as the component number" );

  GEOSX_ERROR_IF( m_defaultViscosity.size() != NC,
                  "The number of default viscosity values is not the same as the component number" );

}

localIndex SlurryFluidBase::numFluidComponents() const
{
  return LvArray::integerConversion< localIndex >( m_componentNames.size());
}

void SlurryFluidBase::allocateConstitutiveData( Group & parent,
                                                localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent.size() );

  localIndex const NC = numFluidComponents();

  m_density.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dDens_dPres.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dDens_dProppantConc.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dDens_dCompConc.resize( parent.size(), numConstitutivePointsPerParentIndex, NC );

  m_componentDensity.resize( parent.size(), numConstitutivePointsPerParentIndex, NC );
  m_dCompDens_dPres.resize( parent.size(), numConstitutivePointsPerParentIndex, NC );
  m_dCompDens_dCompConc.resize( parent.size(), numConstitutivePointsPerParentIndex, NC, NC );


  m_fluidDensity.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dFluidDens_dPres.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dFluidDens_dCompConc.resize( parent.size(), numConstitutivePointsPerParentIndex, NC );

  m_fluidViscosity.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dFluidVisc_dPres.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dFluidVisc_dCompConc.resize( parent.size(), numConstitutivePointsPerParentIndex, NC );

  m_viscosity.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dVisc_dPres.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dVisc_dProppantConc.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dVisc_dCompConc.resize( parent.size(), numConstitutivePointsPerParentIndex, NC );

}


} //namespace constitutive

} //namespace geosx
