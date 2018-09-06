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

namespace
{
#ifdef USE_PVT_PACKAGE
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
#endif
}

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

  RegisterViewWrapper( viewKeys.phaseMoleFraction.Key(), &m_phaseMoleFraction, false );
  RegisterViewWrapper( viewKeys.phaseVolumeFraction.Key(), &m_phaseVolumeFraction, false );
  RegisterViewWrapper( viewKeys.phaseDensity.Key(), &m_phaseDensity, false );
  RegisterViewWrapper( viewKeys.phaseComponentMoleFraction.Key(), &m_phaseCompMoleFraction, false );
  RegisterViewWrapper( viewKeys.phaseComponentDensity.Key(), &m_phaseCompDensity, false );

  RegisterViewWrapper( viewKeys.dPhaseMoleFraction_dPhasePressure.Key(), &m_dPhaseMoleFraction_dPhasePressure, false );
  RegisterViewWrapper( viewKeys.dPhaseMoleFraction_dGlobalCompMoleFraction.Key(), &m_dPhaseMoleFraction_dGlobalCompMoleFraction, false );
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
  localIndex const numPhase = numFluidPhases();
  localIndex const numComp = numFluidComponents();

  this->resize( size );

  m_phaseMoleFraction.resize( size, numPts, numPhase );
  m_phaseVolumeFraction.resize( size, numPts, numPhase );
  m_phaseDensity.resize( size, numPts, numPhase );
  m_phaseCompMoleFraction.resize( size, numPts, numPhase, numComp );
  m_phaseCompDensity.resize( size, numPts, numPhase, numComp );

  m_dPhaseMoleFraction_dPhasePressure.resize( size, numPts, numPhase );
  m_dPhaseMoleFraction_dGlobalCompMoleFraction.resize( size, numPts, numPhase, numComp );
  m_dPhaseVolumeFraction_dPhasePressure.resize( size, numPts, numPhase );
  m_dPhaseVolumeFraction_dGlobalCompMoleFraction.resize( size, numPts, numPhase, numComp );
  m_dPhaseDensity_dPhasePressure.resize( size, numPts, numPhase );
  m_dPhaseDensity_dGlobalCompMoleFraction.resize( size, numPts, numPhase, numComp );
  m_dPhaseCompMoleFraction_dPhasePressure.resize( size, numPts, numPhase, numComp );
  m_dPhaseCompMoleFraction_dGlobalCompMoleFraction.resize( size, numPts, numPhase, numComp, numComp );
  m_dPhaseCompDensity_dPhasePressure.resize( size, numPts, numPhase, numComp );
  m_dPhaseCompDensity_dGlobalCompMoleFraction.resize( size, numPts, numPhase, numComp, numComp );
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

void CompositionalMultiphaseFluid::StateUpdatePointMultiphaseFluid(real64 const & pres,
                                                                   real64 const & temp,
                                                                   real64 const * composition,
                                                                   localIndex const k,
                                                                   localIndex const q)
{
  // 0. set array views to the element/point data to avoid awkward quadruple indexing
  arrayView1d<real64> phaseMoleFrac          = m_phaseMoleFraction[k][q];
  arrayView1d<real64> dPhaseMoleFrac_dPres   = m_dPhaseMoleFraction_dPhasePressure[k][q];
  arrayView2d<real64> dPhaseMoleFrac_dFeed   = m_dPhaseMoleFraction_dGlobalCompMoleFraction[k][q];

  arrayView1d<real64> phaseDens        = m_phaseDensity[k][q];
  arrayView1d<real64> dPhaseDens_dPres = m_dPhaseDensity_dPhasePressure[k][q];
  arrayView2d<real64> dPhaseDens_dFeed = m_dPhaseDensity_dGlobalCompMoleFraction[k][q];

  arrayView1d<real64> phaseVolFrac        = m_phaseVolumeFraction[k][q];
  arrayView1d<real64> dPhaseVolFrac_dPres = m_dPhaseVolumeFraction_dPhasePressure[k][q];
  arrayView2d<real64> dPhaseVolFrac_dFeed = m_dPhaseVolumeFraction_dGlobalCompMoleFraction[k][q];

  arrayView2d<real64> phaseCompMoleFrac        = m_phaseCompMoleFraction[k][q];
  arrayView2d<real64> dPhaseCompMoleFrac_dPres = m_dPhaseCompMoleFraction_dPhasePressure[k][q];
  arrayView3d<real64> dPhaseCompMoleFrac_dFeed = m_dPhaseCompMoleFraction_dGlobalCompMoleFraction[k][q];


#ifdef USE_PVT_PACKAGE
  localIndex const numPhase = numFluidPhases();
  localIndex const numComp = numFluidComponents();

  std::vector<double> feed(numComp);

  for (localIndex ic = 0; ic < numComp; ++ic)
    feed[ic] = composition[ic];

  // 1. Trigger PVTPackage compute and get back phase split

  m_fluid->Update(pres, temp, feed);
  MultiphaseSystemProperties const * split = m_fluid->get_MultiphaseSystemProperties();

  // 2. Extract phase split and properties from PVTPackage and compute volume fractions

  for (localIndex ip = 0; ip < numPhase; ++ip)
  {
    PHASE_TYPE phase = phaseNameDict.at(m_phases[ip]);
    PhaseProperties const * props = m_fluid->get_PhaseProperties(phase);
    std::vector<double> const & xcp = split->MoleComposition.at(phase);

    phaseMoleFrac[ip] = split->MoleFraction.at(phase);
    dPhaseMoleFrac_dPres[ip] = 0.0; // TODO

    phaseDens[ip] = props->MassDensity;
    dPhaseDens_dPres[ip] = 0.0; // TODO

    for (localIndex ic = 0; ic < numComp; ++ic)
    {
      dPhaseMoleFrac_dFeed[ip][ic] = 0.0; // TODO
      dPhaseDens_dFeed[ip][ic] = 0.0; // TODO

      phaseCompMoleFrac[ip][ic] = xcp[ic];
      dPhaseCompMoleFrac_dPres[ip][ic] = 0.0; // TODO

      for (localIndex jc = 0; jc < numComp; ++jc)
        dPhaseCompMoleFrac_dFeed[ip][ic][jc] = 0.0; // TODO
    }
  }

  // 3. Compute phase volume fractions (saturations) and derivatives, requires two passes for normalization

  real64 volumeFractionSum = 0.0;
  real64 dVolumeFractionSum_dPressure = 0.0;
  real64 dVolumeFractionSum_dCompMoleFrac[32]{}; // TODO limit number of components?

  for (localIndex ip = 0; ip < numPhase; ++ip)
  {
    PHASE_TYPE phase = phaseNameDict.at(m_phases[ip]);
    PhaseProperties const * props = m_fluid->get_PhaseProperties(phase);

    phaseVolFrac[ip] = phaseMoleFrac[ip] / props->MoleDensity;

    real64 const dPhaseMolarDensity_dPressure = 0.0; // TODO

    dPhaseVolFrac_dPres[ip] =
      (dPhaseMoleFrac_dPres[ip] - phaseVolFrac[ip] * dPhaseMolarDensity_dPressure) / props->MoleDensity;

    volumeFractionSum += phaseVolFrac[ip];
    dVolumeFractionSum_dPressure += dPhaseVolFrac_dPres[ip];

    for (localIndex ic = 0; ic < numComp; ++ic)
    {
      real64 const dPhaseMolarDens_dFeed = 0.0; // TODO

      dPhaseVolFrac_dFeed[ip][ic] =
        (dPhaseMoleFrac_dFeed[ip][ic] - phaseVolFrac[ip] * dPhaseMolarDens_dFeed) / props->MoleDensity;

      dVolumeFractionSum_dCompMoleFrac[ic] += dPhaseVolFrac_dFeed[ip][ic];
    }
  }
  // normalization pass
  for (localIndex ip = 0; ip < numPhase; ++ip)
  {
    phaseVolFrac[ip] /= volumeFractionSum;
    dPhaseVolFrac_dPres[ip] -= phaseVolFrac[ip] * dVolumeFractionSum_dPressure;
    dPhaseVolFrac_dPres[ip] /= volumeFractionSum;

    for (localIndex ic = 0; ic < numComp; ++ic)
    {
      dPhaseVolFrac_dFeed[ip][ic] -= phaseVolFrac[ip] * dVolumeFractionSum_dCompMoleFrac[ic];
      dPhaseVolFrac_dFeed[ip][ic] /= volumeFractionSum;
    }
  }

  // 4. Compute component-in-phase densities and derivatives

#endif
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompositionalMultiphaseFluid, std::string const &, ManagedGroup * const )
} // namespace constitutive

} // namespace geosx