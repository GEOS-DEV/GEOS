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
 * @file BiotReactivePorosity.cpp
 */

#include "BiotReactivePorosity.hpp"
#include "PorosityFields.hpp"
#include "constitutive/solid/SolidBase.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

BiotReactivePorosity::BiotReactivePorosity( string const & name, Group * const parent ):
  ReactivePorosityBase( name, parent )
{
  registerWrapper( viewKeyStruct::defaultGrainBulkModulusString(), &m_defaultGrainBulkModulus ).
    setInputFlag( InputFlags::REQUIRED ).
    setApplyDefaultValue( -1.0 ).
    setDescription( "Default grain bulk modulus" );

  registerWrapper( viewKeyStruct::defaultMineralBulkModulusString(), &m_defaultMineralBulkModulus ).
    setInputFlag( InputFlags::REQUIRED ).
    setApplyDefaultValue( -1.0 ).
    setDescription( "Default mineral bulk modulus" );

  registerWrapper( viewKeyStruct::mineralBulkModulusString(), &m_mineralBulkModulus ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Mineral bulk modulus" );

  registerWrapper( viewKeyStruct::mineralPressureString(), &m_mineralPressure ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Current mineral pressure" );

  registerWrapper( viewKeyStruct::mineralPressure_nString(), &m_mineralPressure_n ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Mineral pressure at last time step" );

  registerField< fields::porosity::biotCoefficient >( &m_biotCoefficient );

  registerField< fields::porosity::grainBulkModulus >( &m_grainBulkModulus );
}

void BiotReactivePorosity::postInputInitialization()
{
  ReactivePorosityBase::postInputInitialization();

  // set results as array default values
  getWrapper< array1d< real64 > >( fields::porosity::grainBulkModulus::key() ).
    setApplyDefaultValue( m_defaultGrainBulkModulus );

  getWrapper< array1d< real64 > >( viewKeyStruct::mineralBulkModulusString() ).
    setApplyDefaultValue( m_defaultMineralBulkModulus );
}

void BiotReactivePorosity::initializeState() const
{
  ReactivePorosityBase::initializeState();

  m_mineralPressure_n.setValues< parallelDevicePolicy<> >( m_mineralPressure.toViewConst() );
}

void BiotReactivePorosity::saveConvergedState() const
{
  ReactivePorosityBase::saveConvergedState();

  m_mineralPressure_n.setValues< parallelDevicePolicy<> >( m_mineralPressure.toViewConst() );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, BiotReactivePorosity, string const &, Group * const )
} /* namespace constitutive */
} /* namespace geos */
