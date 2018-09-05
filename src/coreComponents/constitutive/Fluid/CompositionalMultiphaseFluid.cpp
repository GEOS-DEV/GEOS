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
  * @file CompositionalMultiphaseFluid.cpp
  */

#include "CompositionalMultiphaseFluid.hpp"

#ifdef USE_PVT_PACKAGE
#include "MultiphaseSystem/CompositionalMultiphaseSystem.hpp"
#endif

using namespace PVTPackage;

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

namespace constitutive
{

CompositionalMultiphaseFluid::CompositionalMultiphaseFluid(std::string const & name, ManagedGroup * const parent)
  : ConstitutiveBase(name, parent),
    m_fluid(nullptr)
{
  RegisterViewWrapper( viewKeys.phases.Key(), &m_phases, false );
  RegisterViewWrapper( viewKeys.equationsOfState.Key(), &m_equationsOfState, false );
  RegisterViewWrapper( viewKeys.componentNames.Key(), &m_componentNames, false );

  RegisterViewWrapper( viewKeys.componentCriticalPressure.Key(), &m_componentCriticalPressure, false );
  RegisterViewWrapper( viewKeys.componentCriticalTemperature.Key(), &m_componentCriticalTemperature, false );
  RegisterViewWrapper( viewKeys.componentAcentricFactor.Key(), &m_componentAcentricFactor, false );
  RegisterViewWrapper( viewKeys.componentMolarWeigth.Key(), &m_componentMolarWeigth, false );
  RegisterViewWrapper( viewKeys.componentVolumeShift.Key(), &m_componentVolumeShift, false );
  RegisterViewWrapper( viewKeys.componentBinaryCoeff.Key(), &m_componentBinaryCoeff, false );

  RegisterViewWrapper( viewKeys.phaseVolumeFraction.Key(), &m_phaseVolumeFraction, false );
  RegisterViewWrapper( viewKeys.phaseDensity.Key(), &m_phaseDensity, false );
  RegisterViewWrapper( viewKeys.phaseComponentMoleFraction.Key(), &m_phaseCompMoleFraction, false );
  RegisterViewWrapper( viewKeys.phaseComponentDensity.Key(), &m_phaseCompDensity, false );

  RegisterViewWrapper( viewKeys.dPhaseVolumeFraction_dPhasePressure.Key(), &m_dPhaseVolumeFraction_dPhasePressure, false );
  RegisterViewWrapper( viewKeys.dPhaseVolumeFraction_dGlobalCompMoleFraction.Key(), &m_dPhaseVolumeFraction_dGlobalCompMoleFraction, false );
  RegisterViewWrapper( viewKeys.dPhaseDensity_dPhasePressure.Key(), &m_dPhaseDensity_dPhasePressure, false );
  RegisterViewWrapper( viewKeys.dPhaseDensity_dGlobalCompMoleFraction.Key(), &m_dPhaseDensity_dGlobalCompMoleFraction, false );
  RegisterViewWrapper( viewKeys.dPhaseCompMoleFraction_dPhasePressure.Key(), &m_dPhaseCompMoleFraction_dPhasePressure, false );
  RegisterViewWrapper( viewKeys.dPhaseCompMoleFraction_dGlobalCompMoleFraction.Key(), &m_dPhaseCompMoleFraction_dGlobalCompMoleFraction, false );
  RegisterViewWrapper( viewKeys.dPhaseCompDensity_dPhasePressure.Key(), &m_dPhaseCompDensity_dPhasePressure, false );
  RegisterViewWrapper( viewKeys.dPhaseCompDensity_dGlobalCompMoleFraction.Key(), &m_dPhaseCompDensity_dGlobalCompMoleFraction, false );
}

CompositionalMultiphaseFluid::~CompositionalMultiphaseFluid()
{
  delete m_fluid;
}

std::unique_ptr<ConstitutiveBase>
CompositionalMultiphaseFluid::DeliverClone(string const & name, ManagedGroup * const parent) const
{
  auto clone = std::make_unique<CompositionalMultiphaseFluid>( name, parent );

  clone->m_phases           = this->m_phases;
  clone->m_equationsOfState = this->m_equationsOfState;
  clone->m_componentNames   = this->m_componentNames;

  clone->m_componentCriticalPressure    = this->m_componentCriticalPressure;
  clone->m_componentCriticalTemperature = this->m_componentCriticalTemperature;
  clone->m_componentAcentricFactor      = this->m_componentAcentricFactor;
  clone->m_componentMolarWeigth         = this->m_componentMolarWeigth;
  clone->m_componentVolumeShift         = this->m_componentVolumeShift;
  clone->m_componentBinaryCoeff         = this->m_componentBinaryCoeff;

  clone->createFluid();

  return clone;
}

void CompositionalMultiphaseFluid::AllocateConstitutiveData(dataRepository::ManagedGroup * const parent,
                                                            localIndex const numPts)
{
  ConstitutiveBase::AllocateConstitutiveData( parent, numPts );

  localIndex const size = parent->size();
  localIndex const np = numFluidPhases();
  localIndex const nc = numFluidComponents();

  this->resize( size );

  m_phaseVolumeFraction.resize( size, numPts, np );
  m_phaseDensity.resize( size, numPts, np );
  m_phaseCompMoleFraction.resize( size, numPts, np, nc );
  m_phaseCompDensity.resize( size, numPts, np, nc );

  m_dPhaseVolumeFraction_dPhasePressure.resize( size, numPts, np );
  m_dPhaseVolumeFraction_dGlobalCompMoleFraction.resize( size, numPts, np, nc );
  m_dPhaseDensity_dPhasePressure.resize( size, numPts, np );
  m_dPhaseDensity_dGlobalCompMoleFraction.resize( size, numPts, np, nc );
  m_dPhaseCompMoleFraction_dPhasePressure.resize( size, numPts, np, nc );
  m_dPhaseCompMoleFraction_dGlobalCompMoleFraction.resize( size, numPts, np, nc, nc );
  m_dPhaseCompDensity_dPhasePressure.resize( size, numPts, np, nc );
  m_dPhaseCompDensity_dGlobalCompMoleFraction.resize( size, numPts, np, nc, nc );
}

void CompositionalMultiphaseFluid::FillDocumentationNode()
{
  DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName( CatalogName() );
  docNode->setSchemaType( "Node" );
  docNode->setShortDescription( "Compositional multiphase fluid model" );

  docNode->AllocateChildNode( viewKeys.phases.Key(),
                              viewKeys.phases.Key(),
                              -1,
                              "real64",
                              "real64",
                              "List of fluid phases",
                              "List of fluid phases",
                              "REQUIRED",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.equationsOfState.Key(),
                              viewKeys.equationsOfState.Key(),
                              -1,
                              "real64",
                              "real64",
                              "List of equation of state types for each phase",
                              "List of equation of state types for each phase",
                              "REQUIRED",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.componentNames.Key(),
                              viewKeys.componentNames.Key(),
                              -1,
                              "real64",
                              "real64",
                              "List of component names",
                              "List of component names",
                              "REQUIRED",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.componentCriticalPressure.Key(),
                              viewKeys.componentCriticalPressure.Key(),
                              -1,
                              "real64",
                              "real64",
                              "Component critical pressures",
                              "Component critical pressures",
                              "REQUIRED",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.componentCriticalTemperature.Key(),
                              viewKeys.componentCriticalTemperature.Key(),
                              -1,
                              "real64",
                              "real64",
                              "Component critical temperatures",
                              "Component critical temperatures",
                              "REQUIRED",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.componentAcentricFactor.Key(),
                              viewKeys.componentAcentricFactor.Key(),
                              -1,
                              "real64",
                              "real64",
                              "Component acentric factors",
                              "Component acentric factors",
                              "REQUIRED",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.componentMolarWeigth.Key(),
                              viewKeys.componentMolarWeigth.Key(),
                              -1,
                              "real64",
                              "real64",
                              "Component molar weights",
                              "Component molar weights",
                              "REQUIRED",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.componentVolumeShift.Key(),
                              viewKeys.componentVolumeShift.Key(),
                              -1,
                              "real64",
                              "real64",
                              "Component volume shifts",
                              "Component volume shifts",
                              "",
                              "",
                              1,
                              1,
                              0 );
}

void CompositionalMultiphaseFluid::ReadXML_PostProcess()
{
#define COMPFLUID_CHECK_INPUT_LENGTH(len, expected, attr) \
  if (integer_conversion<localIndex>(len) != integer_conversion<localIndex>(expected)) \
  { \
    GEOS_ERROR("CompositionalMultiphaseFluid: invalid number of entries in " \
               << (attr) << " attribute (" \
               << (len) << "given, " \
               << (expected) << " expected)"); \
  }

  COMPFLUID_CHECK_INPUT_LENGTH(m_equationsOfState.size(), numFluidPhases(),
                               viewKeys.equationsOfState.Key())

  COMPFLUID_CHECK_INPUT_LENGTH(m_componentCriticalPressure.size(), numFluidComponents(),
                               viewKeys.componentCriticalPressure.Key())

  COMPFLUID_CHECK_INPUT_LENGTH(m_componentCriticalTemperature.size(), numFluidComponents(),
                               viewKeys.componentCriticalTemperature.Key())

  COMPFLUID_CHECK_INPUT_LENGTH(m_componentAcentricFactor.size(), numFluidComponents(),
                               viewKeys.componentAcentricFactor.Key())

  COMPFLUID_CHECK_INPUT_LENGTH(m_componentMolarWeigth.size(), numFluidComponents(),
                               viewKeys.componentMolarWeigth.Key())

  if (m_componentVolumeShift.empty())
  {
    m_componentVolumeShift.resize(numFluidComponents());
    m_componentVolumeShift = 0.0;
  }

  COMPFLUID_CHECK_INPUT_LENGTH(m_componentVolumeShift.size(), numFluidComponents(),
                               viewKeys.componentVolumeShift.Key())

  if (m_componentBinaryCoeff.empty())
  {
    m_componentBinaryCoeff.resize(numFluidComponents() * numFluidComponents());
    m_componentBinaryCoeff = 0.0;
  }

  COMPFLUID_CHECK_INPUT_LENGTH(m_componentBinaryCoeff.size(), numFluidComponents() * numFluidComponents(),
                               viewKeys.componentBinaryCoeff.Key())

#undef CHECK_INPUT_LENGTH
}

void CompositionalMultiphaseFluid::createFluid()
{
#ifdef USE_PVT_PACKAGE
  localIndex const numComp  = numFluidComponents();
  localIndex const numPhase = numFluidPhases();

  std::vector<PHASE_TYPE> phases(numPhase);
  std::vector<EOS_TYPE> eos(numPhase);

  static std::unordered_map<string, PHASE_TYPE> const phaseNameDict =
    {
      { "gas",   PHASE_TYPE::GAS },
      { "oil",   PHASE_TYPE::OIL },
      { "water", PHASE_TYPE::LIQUID_WATER_RICH }
    };

  static std::unordered_map<string, EOS_TYPE> const eosNameDict =
    {
      { "PR",   EOS_TYPE::PENG_ROBINSON },
      { "SRK",  EOS_TYPE::REDLICH_KWONG_SOAVE }
    };

  for (localIndex ip = 0; ip < numPhase; ++ip)
  {
    if (phaseNameDict.find(m_phases[ip]) != phaseNameDict.end())
    {
      phases.push_back(phaseNameDict.at(m_phases[ip]));
    }
    else
    {
      GEOS_ERROR("CompositionalMultiphaseFluid: invalid phase label: " << m_phases[ip]);
    }

    if (eosNameDict.find(m_equationsOfState[ip]) != eosNameDict.end())
    {
      eos.push_back(eosNameDict.at(m_equationsOfState[ip]));
    }
    else
    {
      GEOS_ERROR("CompositionalMultiphaseFluid: invalid EOS label: " << m_equationsOfState[ip]);
    }
  }

  std::vector<std::string> components(m_componentNames.begin(), m_componentNames.end());
  std::vector<double> Pc(m_componentCriticalPressure.begin(), m_componentCriticalPressure.end());
  std::vector<double> Tc(m_componentCriticalTemperature.begin(), m_componentCriticalTemperature.end());
  std::vector<double> Mw(m_componentMolarWeigth.begin(), m_componentMolarWeigth.end());
  std::vector<double> Omega(m_componentAcentricFactor.begin(), m_componentAcentricFactor.end());

  const ComponentProperties CompProps(numComp, components, Mw, Tc, Pc, Omega);
  // TODO choose flash type
  m_fluid = new CompositionalMultiphaseSystem(phases, eos, COMPOSITIONAL_FLASH_TYPE::TRIVIAL, CompProps);
#else
  GEOS_ERROR("Cannot use compositional fluid model without PVTPackage. Rebuild with ENABLE_PVT_PACKAGE=ON");
#endif
}

void CompositionalMultiphaseFluid::InitializePostSubGroups(ManagedGroup * const group)
{
  createFluid();
}

localIndex CompositionalMultiphaseFluid::numFluidComponents()
{
  return integer_conversion<localIndex>(m_componentNames.size());
}

localIndex CompositionalMultiphaseFluid::numFluidPhases()
{
  return integer_conversion<localIndex>(m_phases.size());
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompositionalMultiphaseFluid, std::string const &, ManagedGroup * const )
} // namespace constitutive

} // namespace geosx