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
 * @file SolidInternalEnergy.cpp
 */

#include "SolidInternalEnergy.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

SolidInternalEnergy::SolidInternalEnergy( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_internalEnergy(),
  m_dInternalEnergy_dTemperature(),
  m_specificHeatCapacity(),
  m_referenceTemperature(),
  m_referenceInternalEnergy()
{
  registerWrapper( viewKeyStruct::internalEnergyString(), &m_internalEnergy ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Internal energy of the solid" );

  registerWrapper( viewKeyStruct::oldInternalEnergyString(), &m_oldInternalEnergy ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Internal energy of the solid at the previous time-step" );

  registerWrapper( viewKeyStruct::dInternalEnergy_dTemperatureString(), &m_dInternalEnergy_dTemperature ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Derivative of the solid internal energy w.r.t. temperature" );

  registerWrapper( viewKeyStruct::specificHeatCapacityString(), &m_specificHeatCapacity ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Solid specific heat capacity" );

  registerWrapper( viewKeyStruct::referenceTemperatureString(), &m_referenceTemperature ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Reference temperature" );

  registerWrapper( viewKeyStruct::referenceInternalEnergyString(), &m_referenceInternalEnergy ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Internal energy at the reference temperature" );
}

void SolidInternalEnergy::allocateConstitutiveData( Group & parent,
                                                    localIndex const numConstitutivePointsPerParentIndex )
{
  m_internalEnergy.resize( 0, 1 );
  m_dInternalEnergy_dTemperature.resize( 0, 1 );
  m_oldInternalEnergy.resize( 0, 1 );

  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void SolidInternalEnergy::saveConvergedState() const
{
  arrayView2d< real64 const > internalEnergy = m_internalEnergy;
  arrayView2d< real64 >       oldInternalEnergy = m_oldInternalEnergy;

  forAll< parallelDevicePolicy<> >( internalEnergy.size( 0 ), [=] GEOSX_HOST_DEVICE ( localIndex const k )
  {
    oldInternalEnergy[k][0] = internalEnergy[k][0];
  } );
}

}

}
