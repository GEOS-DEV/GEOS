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
 * @file ThermalBase.cpp
 */

#include "ThermalBase.hpp"

namespace geos
{
using namespace dataRepository;

namespace constitutive
{

ThermalBase::ThermalBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_newTemperature( 0, 0 ),
  m_oldTemperature( 0, 0 ),
  m_defaultHeatCapacity( DBL_MAX ),
  m_heatcapacity( 0 )
{
  registerWrapper( viewKeyStruct::temperatureString(), &m_newTemperature ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setApplyDefaultValue( 0 ). // default to zero initial 0 K
    setDescription( "Current Material Temperature" );

  registerWrapper( viewKeyStruct::oldTemperatureString(), &m_oldTemperature ).
    setApplyDefaultValue( 0 ). // default to zero initial 0 K
    setDescription( "Previous Material Temperature" );

  registerWrapper( viewKeyStruct::heatCapacityString(), &m_heatcapacity ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setApplyDefaultValue( -1 ). // will be overwritten
    setDescription( "Material Heat Capacity" );

  registerWrapper( viewKeyStruct::defaultHeatCapacityString(), &m_defaultHeatCapacity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default Material Heat Capacity" );
}


ThermalBase::~ThermalBase()
{}


void ThermalBase::postInputInitialization()
{
  this->getWrapper< array2d< real64 > >( viewKeyStruct::heatCapacityString() ).
    setApplyDefaultValue( m_defaultHeatCapacity );
}


void ThermalBase::allocateConstitutiveData( dataRepository::Group & parent,
                                              localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_heatcapacity.resize( 0, numConstitutivePointsPerParentIndex );
  m_newTemperature.resize( 0, numConstitutivePointsPerParentIndex );
  m_oldTemperature.resize( 0, numConstitutivePointsPerParentIndex );
}


void ThermalBase::saveConvergedState() const
{
  localIndex const numE = numElem();
  localIndex const numQ = numQuad();

  arrayView3d< real64 const > newTemperature = m_newTemperature;
  arrayView3d< real64 > oldTemperature = m_oldTemperature;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numQ; ++q )
    {
      LvArray::tensorOps::copy< 6 >( oldTemperature[k][q], newTemperature[k][q] );
    }
  } );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermalBase, string const &, Group * const )
} /* namespace constitutive */
} /* namespace geos */
