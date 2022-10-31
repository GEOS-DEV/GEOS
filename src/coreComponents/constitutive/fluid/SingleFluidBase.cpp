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
 * @file SingleFluidBase.cpp
 */

#include "SingleFluidBase.hpp"

#include "SingleFluidExtrinsicData.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

SingleFluidBase::SingleFluidBase( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent )
{
  registerExtrinsicData( extrinsicMeshData::singlefluid::density{}, &m_density );
  registerExtrinsicData( extrinsicMeshData::singlefluid::dDensity_dPressure{}, &m_dDensity_dPressure );
  registerExtrinsicData( extrinsicMeshData::singlefluid::dDensity_dTemperature{}, &m_dDensity_dTemperature );
  registerExtrinsicData( extrinsicMeshData::singlefluid::density_n{}, &m_density_n );

  registerExtrinsicData( extrinsicMeshData::singlefluid::viscosity{}, &m_viscosity );
  registerExtrinsicData( extrinsicMeshData::singlefluid::dViscosity_dPressure{}, &m_dViscosity_dPressure );
  registerExtrinsicData( extrinsicMeshData::singlefluid::dViscosity_dTemperature{}, &m_dViscosity_dTemperature );

  registerExtrinsicData( extrinsicMeshData::singlefluid::internalEnergy{}, &m_internalEnergy );
  registerExtrinsicData( extrinsicMeshData::singlefluid::internalEnergy_n{}, &m_internalEnergy_n );
  registerExtrinsicData( extrinsicMeshData::singlefluid::dInternalEnergy_dPressure{}, &m_dInternalEnergy_dPressure );
  registerExtrinsicData( extrinsicMeshData::singlefluid::dInternalEnergy_dTemperature{}, &m_dInternalEnergy_dTemperature );

  registerExtrinsicData( extrinsicMeshData::singlefluid::enthalpy{}, &m_enthalpy );
  registerExtrinsicData( extrinsicMeshData::singlefluid::dEnthalpy_dPressure{}, &m_dEnthalpy_dPressure );
  registerExtrinsicData( extrinsicMeshData::singlefluid::dEnthalpy_dTemperature{}, &m_dEnthalpy_dTemperature );

}

void SingleFluidBase::postProcessInput()
{
  ConstitutiveBase::postProcessInput();

  // for fracture elements, set the default value
  getField< extrinsicMeshData::singlefluid::density_n >().
    setDefaultValue( defaultDensity() );
}

void SingleFluidBase::initializeState() const
{
  saveConvergedState();
}

void SingleFluidBase::saveConvergedState() const
{
  m_density_n.setValues< parallelDevicePolicy<> >( m_density.toViewConst() );
  m_internalEnergy_n.setValues< parallelDevicePolicy<> >( m_internalEnergy.toViewConst() );
}

//START_SPHINX_INCLUDE_00
void SingleFluidBase::allocateConstitutiveData( Group & parent,
                                                localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  resize( parent.size() );

  m_density.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dDensity_dPressure.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dDensity_dTemperature.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_density_n.resize( parent.size(), numConstitutivePointsPerParentIndex );

  m_viscosity.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dViscosity_dPressure.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dViscosity_dTemperature.resize( parent.size(), numConstitutivePointsPerParentIndex );

  m_internalEnergy.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_internalEnergy_n.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dInternalEnergy_dPressure.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dInternalEnergy_dTemperature.resize( parent.size(), numConstitutivePointsPerParentIndex );

  m_enthalpy.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dEnthalpy_dPressure.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dEnthalpy_dTemperature.resize( parent.size(), numConstitutivePointsPerParentIndex );
}
//END_SPHINX_INCLUDE_00

} //namespace constitutive

} //namespace geosx
