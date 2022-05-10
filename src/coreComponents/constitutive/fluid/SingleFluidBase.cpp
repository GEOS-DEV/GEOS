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
  registerExtrinsicData( extrinsicMeshData::singlefluid::initialDensity{}, &m_initialDensity );
  registerExtrinsicData( extrinsicMeshData::singlefluid::density_n{}, &m_density_n );

  registerExtrinsicData( extrinsicMeshData::singlefluid::viscosity{}, &m_viscosity );
  registerExtrinsicData( extrinsicMeshData::singlefluid::dViscosity_dPressure{}, &m_dViscosity_dPressure );

}

void SingleFluidBase::postProcessInput()
{
  ConstitutiveBase::postProcessInput();

  // for fracture elements, set the default value
  getExtrinsicData< extrinsicMeshData::singlefluid::density_n >().
    setDefaultValue( defaultDensity() );
}

void SingleFluidBase::initializeState() const
{
  arrayView2d< real64 const > density        = m_density;
  arrayView2d< real64 >       initialDensity = m_initialDensity;
  arrayView2d< real64 >       density_n     = m_density_n;
  initialDensity.setValues< parallelDevicePolicy<> >( density );
  density_n.setValues< parallelDevicePolicy<> >( density );
}

void SingleFluidBase::saveConvergedState() const
{
  arrayView2d< real64 const > density   = m_density;
  arrayView2d< real64 >       density_n = m_density_n;
  density_n.setValues< parallelDevicePolicy<> >( density );
}

//START_SPHINX_INCLUDE_00
void SingleFluidBase::allocateConstitutiveData( Group & parent,
                                                localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  resize( parent.size() );

  m_density.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dDensity_dPressure.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_initialDensity.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_density_n.resize( parent.size(), numConstitutivePointsPerParentIndex );

  m_viscosity.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dViscosity_dPressure.resize( parent.size(), numConstitutivePointsPerParentIndex );
}
//END_SPHINX_INCLUDE_00

} //namespace constitutive

} //namespace geosx
