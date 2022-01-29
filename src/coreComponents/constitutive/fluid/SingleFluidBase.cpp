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

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

SingleFluidBase::SingleFluidBase( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::densityString(), &m_density ).setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dDens_dPresString(), &m_dDensity_dPressure );

  registerWrapper( viewKeyStruct::initialDensityString(), &m_initialDensity );

  registerWrapper( viewKeyStruct::viscosityString(), &m_viscosity ).setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dVisc_dPresString(), &m_dViscosity_dPressure );
}

void SingleFluidBase::postProcessInput()
{
  ConstitutiveBase::postProcessInput();
}

void SingleFluidBase::initializeState() const
{
  localIndex const numE = m_density.size( 0 );
  localIndex const numG = m_density.size( 1 );

  arrayView2d< real64 const > newDensity     = m_density;
  arrayView2d< real64 >       initialDensity = m_initialDensity;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOSX_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numG; ++q )
    {
      initialDensity[k][q] = newDensity[k][q];
    }
  } );
}

void SingleFluidBase::allocateConstitutiveData( Group & parent,
                                                localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  resize( parent.size() );

  m_density.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dDensity_dPressure.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_initialDensity.resize( parent.size(), numConstitutivePointsPerParentIndex );

  m_viscosity.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dViscosity_dPressure.resize( parent.size(), numConstitutivePointsPerParentIndex );
}

} //namespace constitutive

} //namespace geosx
