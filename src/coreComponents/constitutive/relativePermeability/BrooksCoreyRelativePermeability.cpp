/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file BrooksCoreyRelativePermeability.cpp
 */

#include "BrooksCoreyRelativePermeability.hpp"

#include <cmath>

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

//START_SPHINX_INCLUDE_00

BrooksCoreyRelativePermeability::BrooksCoreyRelativePermeability( string const & name,
                                                                  Group * const parent )
  : RelativePermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::defaultPhaseMinVolumeFractionString(), &m_defaultPhaseMinVolumeFraction ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default minimum volume fraction value for each phase" );

  registerWrapper( viewKeyStruct::defaultPhaseRelPermExponentString(), &m_defaultPhaseRelPermExponent ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default minimum rel perm power law exponent for each phase" );

  registerWrapper( viewKeyStruct::defaultPhaseRelPermMaxValueString(), &m_defaultPhaseRelPermMaxValue ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default maximum rel perm value for each phase" );

  registerWrapper( viewKeyStruct::phaseMinVolumeFractionString(), &m_phaseMinVolumeFraction ).
    setSizedFromParent( 1 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Minimum volume fraction value for each phase" );

  registerWrapper( viewKeyStruct::phaseRelPermExponentString(), &m_phaseRelPermExponent ).
    setSizedFromParent( 1 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Minimum rel perm power law exponent for each phase" );

  registerWrapper( viewKeyStruct::phaseRelPermMaxValueString(), &m_phaseRelPermMaxValue ).
    setSizedFromParent( 1 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Maximum rel perm value for each phase" );

  registerWrapper( viewKeyStruct::volFracScaleString(), &m_volFracScale ).
    setSizedFromParent( 1 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Factor used to scale the phase capillary pressure, defined as: one minus the sum of the phase minimum volume fractions." );

}

BrooksCoreyRelativePermeability::~BrooksCoreyRelativePermeability()
{}


void BrooksCoreyRelativePermeability::postProcessInput()
{
  RelativePermeabilityBase::postProcessInput();

  localIndex const numPhases = numFluidPhases();

  #define COREY_CHECK_INPUT_LENGTH( data, expected, attr ) \
    if( LvArray::integerConversion< localIndex >((data).size()) != LvArray::integerConversion< localIndex >( expected )) \
    { \
      GEOSX_ERROR( "BrooksCoreyRelativePermeability: invalid number of entries in " \
                   << (attr) << " attribute (" \
                   << (data).size() << "given, " \
                   << (expected) << " expected)" ); \
    }

  COREY_CHECK_INPUT_LENGTH( m_defaultPhaseMinVolumeFraction, numPhases, viewKeyStruct::defaultPhaseMinVolumeFractionString() )
  COREY_CHECK_INPUT_LENGTH( m_defaultPhaseRelPermExponent, numPhases, viewKeyStruct::defaultPhaseRelPermExponentString() )
  COREY_CHECK_INPUT_LENGTH( m_defaultPhaseRelPermMaxValue, numPhases, viewKeyStruct::defaultPhaseRelPermMaxValueString() )

#undef COREY_CHECK_INPUT_LENGTH

  for( localIndex ip = 0; ip < numPhases; ++ip )
  {
    GEOSX_ERROR_IF( m_defaultPhaseMinVolumeFraction[ip] < 0.0 || m_defaultPhaseMinVolumeFraction[ip] > 1.0,
                    "BrooksCoreyRelativePermeability: invalid min volume fraction value: " << m_defaultPhaseMinVolumeFraction[ip] );

    GEOSX_ERROR_IF( m_defaultPhaseRelPermExponent[ip] < 0.0,
                    "BrooksCoreyRelativePermeability: invalid exponent value: " << m_defaultPhaseRelPermExponent[ip] );

    GEOSX_ERROR_IF( m_defaultPhaseRelPermMaxValue[ip] < 0.0 || m_defaultPhaseRelPermMaxValue[ip] > 1.0,
                    "BrooksCoreyRelativePermeability: invalid maximum value: " << m_defaultPhaseRelPermMaxValue[ip] );
  }

}

void BrooksCoreyRelativePermeability::resizeFields( localIndex const size, localIndex const numPts )
{
  RelativePermeabilityBase::resizeFields( size, numPts );

  m_phaseMinVolumeFraction.resize( size, numFluidPhases() );
  m_phaseRelPermExponent.resize( size, numFluidPhases() );
  m_phaseRelPermMaxValue.resize( size, numFluidPhases() );
  m_volFracScale.resize( size );

  // TODO: forAll
  for( localIndex ei = 0; ei < size; ++ei )
  {
    m_volFracScale[ei] = 1.0;
    for( localIndex ip = 0; ip < numFluidPhases(); ++ip )
    {
      m_phaseMinVolumeFraction[ei][ip] = m_defaultPhaseMinVolumeFraction[ip];
      m_phaseRelPermExponent[ei][ip] = m_defaultPhaseRelPermExponent[ip];
      m_phaseRelPermMaxValue[ei][ip] = m_defaultPhaseRelPermMaxValue[ip];
      m_volFracScale[ei] -= m_phaseMinVolumeFraction[ei][ip];
    }
    GEOSX_ERROR_IF( m_volFracScale[ei] < 0.0, "BrooksCoreyRelativePermeability: sum of min volume fractions exceeds 1.0" );
  }
}

void BrooksCoreyRelativePermeability::initializePostInitialConditionsPreSubGroups()
{
  // TODO: forAll
  for( localIndex ei = 0; ei < size(); ++ei )
  {
    m_volFracScale[ei] = 1.0;
    for( localIndex ip = 0; ip < numFluidPhases(); ++ip )
    {
      m_volFracScale[ei] -= m_phaseMinVolumeFraction[ei][ip];
    }
  }
}

BrooksCoreyRelativePermeability::KernelWrapper BrooksCoreyRelativePermeability::createKernelWrapper()
{
  return KernelWrapper( m_phaseMinVolumeFraction,
                        m_phaseRelPermExponent,
                        m_phaseRelPermMaxValue,
                        m_volFracScale,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_phaseRelPerm,
                        m_dPhaseRelPerm_dPhaseVolFrac );
}

//START_SPHINX_INCLUDE_01
REGISTER_CATALOG_ENTRY( ConstitutiveBase, BrooksCoreyRelativePermeability, string const &, Group * const )

} // namespace constitutive

} // namespace geosx
