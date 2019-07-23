/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
  * @file RelativePermeabilityBase.cpp
  */

#include "RelativePermeabilityBase.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

namespace constitutive
{

constexpr integer RelativePermeabilityBase::PhaseType::GAS;
constexpr integer RelativePermeabilityBase::PhaseType::OIL;
constexpr integer RelativePermeabilityBase::PhaseType::WATER;

namespace
{

std::unordered_map<string, integer> const phaseDict =
{
  { "gas",   RelativePermeabilityBase::PhaseType::GAS   },
  { "oil",   RelativePermeabilityBase::PhaseType::OIL   },
  { "water", RelativePermeabilityBase::PhaseType::WATER }
};

}


RelativePermeabilityBase::RelativePermeabilityBase( std::string const & name, ManagedGroup * const parent )
  : ConstitutiveBase( name, parent )
{
  RegisterViewWrapper( viewKeyStruct::phaseNamesString, &m_phaseNames, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("List of fluid phases");

  RegisterViewWrapper( viewKeyStruct::phaseTypesString, &m_phaseTypes, false )->
    setSizedFromParent(0);

  RegisterViewWrapper( viewKeyStruct::phaseOrderString, &m_phaseOrder, false )->
    setSizedFromParent(0);

  RegisterViewWrapper( viewKeyStruct::phaseRelPermString, &m_phaseRelPerm, false )->setPlotLevel( PlotLevel::LEVEL_0 );
  RegisterViewWrapper( viewKeyStruct::dPhaseRelPerm_dPhaseVolFractionString, &m_dPhaseRelPerm_dPhaseVolFrac, false );
}

RelativePermeabilityBase::~RelativePermeabilityBase()
{

}


void RelativePermeabilityBase::PostProcessInput()
{
  ConstitutiveBase::PostProcessInput();

  localIndex const NP = numFluidPhases();

  GEOS_ERROR_IF( NP < 2, "RelativePermeabilityBase: number of fluid phases should be at least 2" );

  GEOS_ERROR_IF( NP > PhaseType::MAX_NUM_PHASES,
                 "RelativePermeabilityBase: number of fluid phases exceeds the maximum of " << PhaseType::MAX_NUM_PHASES );

  m_phaseTypes.resize( NP );
  m_phaseOrder.resize( PhaseType::MAX_NUM_PHASES );
  m_phaseOrder = -1;

  for (localIndex ip = 0; ip < NP; ++ip)
  {
    auto it = phaseDict.find( m_phaseNames[ip] );
    GEOS_ERROR_IF( it == phaseDict.end(), "RelativePermeabilityBase: phase not supported: " << m_phaseNames[ip] );
    integer const phaseIndex = it->second;
    GEOS_ERROR_IF( phaseIndex >= PhaseType::MAX_NUM_PHASES, "RelativePermeabilityBase: invalid phase index " << phaseIndex );

    m_phaseTypes[ip] = phaseIndex;
    m_phaseOrder[phaseIndex] = integer_conversion<integer>(ip);
  }

  // call to correctly set member array tertiary sizes on the 'main' material object
  ResizeFields( 0, 0 );
}

void RelativePermeabilityBase::ResizeFields( localIndex const size, localIndex const numPts )
{
  localIndex const NP = numFluidPhases();

  m_phaseRelPerm.resize( size, numPts, NP );
  m_dPhaseRelPerm_dPhaseVolFrac.resize( size, numPts, NP, NP );
}

void RelativePermeabilityBase::AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                                         localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
  ResizeFields( parent->size(), numConstitutivePointsPerParentIndex );
}

localIndex RelativePermeabilityBase::numFluidPhases() const
{
  return integer_conversion<localIndex>(m_phaseNames.size());
}

string const & RelativePermeabilityBase::phaseName( localIndex ip ) const
{
  GEOS_ERROR_IF( ip >= numFluidPhases(), "Index " << ip << " exceeds number of fluid phases" );
  return m_phaseNames[ip];
}

} // namespace constitutive

} // namespace geosx
