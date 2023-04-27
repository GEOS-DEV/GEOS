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

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

SolidInternalEnergy::SolidInternalEnergy( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_internalEnergy(),
  m_dInternalEnergy_dTemperature(),
  m_volumetricHeatCapacity(),
  m_referenceTemperature(),
  m_referenceInternalEnergy()
{
  registerWrapper( viewKeyStruct::internalEnergyString(), &m_internalEnergy ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Internal energy of the solid per unit volume [J/m^3]" );

  registerWrapper( viewKeyStruct::oldInternalEnergyString(), &m_internalEnergy_n ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Internal energy of the solid per unit volume at the previous time-step [J/m^3]" );

  registerWrapper( viewKeyStruct::dInternalEnergy_dTemperatureString(), &m_dInternalEnergy_dTemperature ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Derivative of the solid internal energy w.r.t. temperature [J/(m^3.K)]" );

  registerWrapper( viewKeyStruct::volumetricHeatCapacityString(), &m_volumetricHeatCapacity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Solid volumetric heat capacity [J/(kg.K)]" );

  registerWrapper( viewKeyStruct::referenceTemperatureString(), &m_referenceTemperature ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference temperature [K]" );

  registerWrapper( viewKeyStruct::referenceInternalEnergyString(), &m_referenceInternalEnergy ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Internal energy at the reference temperature [J/kg]" );
}

void SolidInternalEnergy::allocateConstitutiveData( Group & parent,
                                                    localIndex const numConstitutivePointsPerParentIndex )
{
  m_internalEnergy.resize( 0, 1 );
  m_dInternalEnergy_dTemperature.resize( 0, 1 );
  m_internalEnergy_n.resize( 0, 1 );

  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void SolidInternalEnergy::saveConvergedState() const
{
  arrayView2d< real64 const > internalEnergy   = m_internalEnergy;
  arrayView2d< real64 >       internalEnergy_n = m_internalEnergy_n;

  forAll< parallelDevicePolicy<> >( internalEnergy.size( 0 ), [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    internalEnergy_n[k][0] = internalEnergy[k][0];
  } );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, SolidInternalEnergy, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
