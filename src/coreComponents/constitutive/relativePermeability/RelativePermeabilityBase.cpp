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
 * @file RelativePermeabilityBase.cpp
 */

#include "RelativePermeabilityBase.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

constexpr integer RelativePermeabilityBase::PhaseType::GAS;
constexpr integer RelativePermeabilityBase::PhaseType::OIL;
constexpr integer RelativePermeabilityBase::PhaseType::WATER;

namespace
{

std::unordered_map< string, integer > const phaseDict =
{
  { "gas", RelativePermeabilityBase::PhaseType::GAS   },
  { "oil", RelativePermeabilityBase::PhaseType::OIL   },
  { "water", RelativePermeabilityBase::PhaseType::WATER }
};

}


RelativePermeabilityBase::RelativePermeabilityBase( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseNamesString(), &m_phaseNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "List of fluid phases" );

  registerWrapper( viewKeyStruct::phaseTypesString(), &m_phaseTypes ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::phaseOrderString(), &m_phaseOrder ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::phaseRelPermString(), &m_phaseRelPerm ).setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dPhaseRelPerm_dPhaseVolFractionString(), &m_dPhaseRelPerm_dPhaseVolFrac );
}

RelativePermeabilityBase::~RelativePermeabilityBase()
{}


void RelativePermeabilityBase::postProcessInput()
{
  ConstitutiveBase::postProcessInput();

  localIndex const numPhases = numFluidPhases();

  GEOSX_THROW_IF( numPhases < 2, "RelativePermeabilityBase: number of fluid phases should be at least 2", InputError );

  GEOSX_THROW_IF( numPhases > PhaseType::MAX_NUM_PHASES,
                  "RelativePermeabilityBase: number of fluid phases exceeds the maximum of " << PhaseType::MAX_NUM_PHASES,
                  InputError );

  m_phaseTypes.resize( numPhases );
  m_phaseOrder.resize( PhaseType::MAX_NUM_PHASES );
  m_phaseOrder.setValues< serialPolicy >( -1 );

  for( localIndex ip = 0; ip < numPhases; ++ip )
  {
    auto it = phaseDict.find( m_phaseNames[ip] );
    GEOSX_THROW_IF( it == phaseDict.end(), "RelativePermeabilityBase: phase not supported: " << m_phaseNames[ip], InputError );
    integer const phaseIndex = it->second;
    GEOSX_THROW_IF( phaseIndex >= PhaseType::MAX_NUM_PHASES, "RelativePermeabilityBase: invalid phase index " << phaseIndex, InputError );

    m_phaseTypes[ip] = phaseIndex;
    m_phaseOrder[phaseIndex] = LvArray::integerConversion< integer >( ip );
  }

  // call to correctly set member array tertiary sizes on the 'main' material object
  resizeFields( 0, 0 );
}

void RelativePermeabilityBase::resizeFields( localIndex const size, localIndex const numPts )
{
  localIndex const numPhases = numFluidPhases();

  m_phaseRelPerm.resize( size, numPts, numPhases );
  m_dPhaseRelPerm_dPhaseVolFrac.resize( size, numPts, numPhases, numPhases );
}

void RelativePermeabilityBase::allocateConstitutiveData( dataRepository::Group & parent,
                                                         localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
  resizeFields( parent.size(), numConstitutivePointsPerParentIndex );
}

} // namespace constitutive

} // namespace geosx
