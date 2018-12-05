/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
  RegisterViewWrapper( viewKeyStruct::phaseNamesString, &m_phaseNames, false );
  RegisterViewWrapper( viewKeyStruct::phaseTypesString, &m_phaseTypes, false );
  RegisterViewWrapper( viewKeyStruct::phaseOrderString, &m_phaseOrder, false );

  RegisterViewWrapper( viewKeyStruct::phaseRelPermString, &m_phaseRelPerm, false );
  RegisterViewWrapper( viewKeyStruct::dPhaseRelPerm_dPhaseVolFractionString, &m_dPhaseRelPerm_dPhaseVolFrac, false );
}

RelativePermeabilityBase::~RelativePermeabilityBase()
{

}

void RelativePermeabilityBase::FillDocumentationNode()
{
  DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName( this->GetCatalogName() );
  docNode->setSchemaType( "Node" );
  docNode->setShortDescription( "Relative permeability model" );

  docNode->AllocateChildNode( viewKeyStruct::phaseNamesString,
                              viewKeyStruct::phaseNamesString,
                              -1,
                              "string_array",
                              "string_array",
                              "List of fluid phases",
                              "List of fluid phases",
                              "REQUIRED",
                              "",
                              1,
                              1,
                              0 );
}

void RelativePermeabilityBase::ReadXML_PostProcess()
{
  ConstitutiveBase::ReadXML_PostProcess();

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
}

void RelativePermeabilityBase::AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                                         localIndex const numPts )
{
  ConstitutiveBase::AllocateConstitutiveData(parent, numPts);

  localIndex const size = parent->size();
  localIndex const NP = numFluidPhases();

  m_phaseRelPerm.resize( size, numPts, NP );
  m_dPhaseRelPerm_dPhaseVolFrac.resize( size, numPts, NP, NP );
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
