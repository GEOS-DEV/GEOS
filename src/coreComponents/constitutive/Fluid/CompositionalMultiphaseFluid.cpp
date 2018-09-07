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
  RegisterViewWrapper( viewKeys.componentMolarWeight.Key(), &m_componentMolarWeigth, false );
  RegisterViewWrapper( viewKeys.componentVolumeShift.Key(), &m_componentVolumeShift, false );
  RegisterViewWrapper( viewKeys.componentBinaryCoeff.Key(), &m_componentBinaryCoeff, false );


  RegisterViewWrapper( viewKeys.phaseMoleFraction.Key(), &m_phaseMoleFraction, false );
  RegisterViewWrapper( viewKeys.dPhaseMoleFraction_dPressure.Key(), &m_dPhaseMoleFraction_dPressure, false );
  RegisterViewWrapper( viewKeys.dPhaseMoleFraction_dTemperature.Key(), &m_dPhaseMoleFraction_dTemperature, false );
  RegisterViewWrapper( viewKeys.dPhaseMoleFraction_dGlobalCompMoleFraction.Key(), &m_dPhaseMoleFraction_dGlobalCompMoleFraction, false );

  RegisterViewWrapper( viewKeys.phaseVolumeFraction.Key(), &m_phaseVolumeFraction, false );
  RegisterViewWrapper( viewKeys.dPhaseVolumeFraction_dPressure.Key(), &m_dPhaseVolumeFraction_dPressure, false );
  RegisterViewWrapper( viewKeys.dPhaseVolumeFraction_dTemperature.Key(), &m_dPhaseVolumeFraction_dTemperature, false );
  RegisterViewWrapper( viewKeys.dPhaseVolumeFraction_dGlobalCompMoleFraction.Key(), &m_dPhaseVolumeFraction_dGlobalCompMoleFraction, false );

  RegisterViewWrapper( viewKeys.phaseDensity.Key(), &m_phaseDensity, false );
  RegisterViewWrapper( viewKeys.dPhaseDensity_dPressure.Key(), &m_dPhaseDensity_dPressure, false );
  RegisterViewWrapper( viewKeys.dPhaseDensity_dTemperature.Key(), &m_dPhaseDensity_dTemperature, false );
  RegisterViewWrapper( viewKeys.dPhaseDensity_dGlobalCompMoleFraction.Key(), &m_dPhaseDensity_dGlobalCompMoleFraction, false );

  RegisterViewWrapper( viewKeys.phaseComponentMoleFraction.Key(), &m_phaseCompMoleFraction, false );
  RegisterViewWrapper( viewKeys.dPhaseCompMoleFraction_dPressure.Key(), &m_dPhaseCompMoleFraction_dPressure, false );
  RegisterViewWrapper( viewKeys.dPhaseCompMoleFraction_dTemperature.Key(), &m_dPhaseCompMoleFraction_dTemperature, false );
  RegisterViewWrapper( viewKeys.dPhaseCompMoleFraction_dGlobalCompMoleFraction.Key(), &m_dPhaseCompMoleFraction_dGlobalCompMoleFraction, false );

  RegisterViewWrapper( viewKeys.phaseComponentDensity.Key(), &m_phaseCompDensity, false );
  RegisterViewWrapper( viewKeys.dPhaseCompDensity_dPressure.Key(), &m_dPhaseCompDensity_dPressure, false );
  RegisterViewWrapper( viewKeys.dPhaseCompDensity_dTemperature.Key(), &m_dPhaseCompDensity_dTemperature, false );
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
  m_dPhaseMoleFraction_dPressure.resize( size, numPts, numPhase );
  m_dPhaseMoleFraction_dTemperature.resize( size, numPts, numPhase );
  m_dPhaseMoleFraction_dGlobalCompMoleFraction.resize( size, numPts, numPhase, numComp );

  m_phaseVolumeFraction.resize( size, numPts, numPhase );
  m_dPhaseVolumeFraction_dPressure.resize( size, numPts, numPhase );
  m_dPhaseVolumeFraction_dTemperature.resize( size, numPts, numPhase );
  m_dPhaseVolumeFraction_dGlobalCompMoleFraction.resize( size, numPts, numPhase, numComp );

  m_phaseDensity.resize( size, numPts, numPhase );
  m_dPhaseDensity_dPressure.resize( size, numPts, numPhase );
  m_dPhaseDensity_dTemperature.resize( size, numPts, numPhase );
  m_dPhaseDensity_dGlobalCompMoleFraction.resize( size, numPts, numPhase, numComp );

  m_phaseCompMoleFraction.resize( size, numPts, numPhase, numComp );
  m_dPhaseCompMoleFraction_dPressure.resize( size, numPts, numPhase, numComp );
  m_dPhaseCompMoleFraction_dTemperature.resize( size, numPts, numPhase, numComp );
  m_dPhaseCompMoleFraction_dGlobalCompMoleFraction.resize( size, numPts, numPhase, numComp, numComp );

  m_phaseCompDensity.resize( size, numPts, numPhase, numComp );
  m_dPhaseCompDensity_dPressure.resize( size, numPts, numPhase, numComp );
  m_dPhaseCompDensity_dTemperature.resize( size, numPts, numPhase, numComp );
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

  docNode->AllocateChildNode( viewKeys.componentMolarWeight.Key(),
                              viewKeys.componentMolarWeight.Key(),
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
  if (numFluidComponents() > MAX_NUM_COMPONENTS)
  {
    GEOS_ERROR("Number of components exceeds the maximum of " << MAX_NUM_COMPONENTS);
  }

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
                               viewKeys.componentMolarWeight.Key())

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
#ifdef USE_PVT_PACKAGE
  // 0. set array views to the element/point data to avoid awkward quadruple indexing
  VarContainer<1> phaseMoleFrac {
    m_phaseMoleFraction[k][q],
    m_dPhaseMoleFraction_dPressure[k][q],
    m_dPhaseMoleFraction_dTemperature[k][q],
    m_dPhaseMoleFraction_dGlobalCompMoleFraction[k][q]
  };

  VarContainer<1> phaseDens {
    m_phaseDensity[k][q],
    m_dPhaseDensity_dPressure[k][q],
    m_dPhaseDensity_dTemperature[k][q],
    m_dPhaseDensity_dGlobalCompMoleFraction[k][q]
  };

  VarContainer<1> phaseVolFrac {
    m_phaseVolumeFraction[k][q],
    m_dPhaseVolumeFraction_dPressure[k][q],
    m_dPhaseVolumeFraction_dTemperature[k][q],
    m_dPhaseVolumeFraction_dGlobalCompMoleFraction[k][q]
  };

  VarContainer<2> phaseCompMoleFrac {
    m_phaseCompMoleFraction[k][q],
    m_dPhaseCompMoleFraction_dPressure[k][q],
    m_dPhaseCompMoleFraction_dTemperature[k][q],
    m_dPhaseCompMoleFraction_dGlobalCompMoleFraction[k][q]
  };

  localIndex const numPhase = numFluidPhases();
  localIndex const numComp = numFluidComponents();

  std::vector<double> feed(numComp);

  for (localIndex ic = 0; ic < numComp; ++ic)
    feed[ic] = composition[ic];

  // 1. Trigger PVTPackage compute and get back phase split

  m_fluid->Update(pres, temp, feed);
  MultiphaseSystemProperties const * split = m_fluid->get_MultiphaseSystemProperties();

  // 2. Extract phase split, phase properties and derivatives from PVTPackage

  for (localIndex ip = 0; ip < numPhase; ++ip)
  {
    PHASE_TYPE phase = phaseNameDict.at(m_phases[ip]);
    PhaseProperties const * props = m_fluid->get_PhaseProperties(phase);
    std::vector<double> const & xcp = split->MoleComposition.at(phase);

    phaseMoleFrac.value[ip] = split->MoleFraction.at(phase);
    phaseMoleFrac.dPres[ip] = 0.0; // TODO
    phaseMoleFrac.dTemp[ip] = 0.0; // TODO

    phaseDens.value[ip] = props->MassDensity;
    phaseDens.dPres[ip] = 0.0; // TODO
    phaseDens.dTemp[ip] = 0.0; // TODO

    for (localIndex ic = 0; ic < numComp; ++ic)
    {
      phaseMoleFrac.dComp[ip][ic] = 0.0; // TODO
      phaseDens.dComp[ip][ic] = 0.0; // TODO

      phaseCompMoleFrac.value[ip][ic] = xcp[ic];
      phaseCompMoleFrac.dPres[ip][ic] = 0.0; // TODO
      phaseCompMoleFrac.dTemp[ip][ic] = 0.0; // TODO

      for (localIndex jc = 0; jc < numComp; ++jc)
        phaseCompMoleFrac.dComp[ip][ic][jc] = 0.0; // TODO
    }
  }

  // 3. Compute phase volume fractions and derivatives, requires two passes for normalization

  real64 volFracSum = 0.0;
  real64 dVolFracSum_dPres = 0.0;
  real64 dVolFracSum_dTemp = 0.0;
  real64 dVolFracSum_dComp[MAX_NUM_COMPONENTS]{}; // TODO limit number of components or invert loop structure?

  for (localIndex ip = 0; ip < numPhase; ++ip)
  {
    PHASE_TYPE phase = phaseNameDict.at(m_phases[ip]);
    PhaseProperties const * props = m_fluid->get_PhaseProperties(phase);

    phaseVolFrac.value[ip] = phaseMoleFrac.value[ip] / props->MoleDensity;

    real64 const dPhaseMolarDens_dPres = 0.0; // TODO
    real64 const dPhaseMolarDens_dTemp = 0.0; // TODO

    phaseVolFrac.dPres[ip] =
      (phaseMoleFrac.dPres[ip] - phaseVolFrac.value[ip] * dPhaseMolarDens_dPres) / props->MoleDensity;

    phaseVolFrac.dTemp[ip] =
      (phaseMoleFrac.dTemp[ip] - phaseVolFrac.value[ip] * dPhaseMolarDens_dTemp) / props->MoleDensity;

    volFracSum += phaseVolFrac.value[ip];
    dVolFracSum_dPres += phaseVolFrac.dPres[ip];
    dVolFracSum_dTemp += phaseVolFrac.dTemp[ip];

    for (localIndex ic = 0; ic < numComp; ++ic)
    {
      real64 const dPhaseMolarDens_dFeed = 0.0; // TODO

      phaseVolFrac.dComp[ip][ic] =
        (phaseMoleFrac.dComp[ip][ic] - phaseVolFrac.value[ip] * dPhaseMolarDens_dFeed) / props->MoleDensity;

      dVolFracSum_dComp[ic] += phaseVolFrac.dComp[ip][ic];
    }
  }
  // normalization pass
  for (localIndex ip = 0; ip < numPhase; ++ip)
  {
    phaseVolFrac.value[ip] /= volFracSum;
    phaseVolFrac.dPres[ip] = (phaseVolFrac.dPres[ip] - phaseVolFrac.value[ip] * dVolFracSum_dPres) / volFracSum;
    phaseVolFrac.dTemp[ip] = (phaseVolFrac.dTemp[ip] - phaseVolFrac.value[ip] * dVolFracSum_dTemp) / volFracSum;

    for (localIndex ic = 0; ic < numComp; ++ic)
    {
      phaseVolFrac.dComp[ip][ic] -= phaseVolFrac.value[ip] * dVolFracSum_dComp[ic];
      phaseVolFrac.dComp[ip][ic] /= volFracSum;
    }
  }

  // 4. Compute component-in-phase densities and derivatives

  for (localIndex ip = 0; ip < numPhase; ++ip)
  {
    PHASE_TYPE phase = phaseNameDict.at(m_phases[ip]);
    PhaseProperties const * props = m_fluid->get_PhaseProperties(phase);

    real64 const dPhaseMolarWeight_dPres = 0.0; // TODO
    real64 const dPhaseMolarWeight_dTemp = 0.0; // TODO

    for (localIndex ic = 0; ic < numComp; ++ic)
    {
      real64 const mw = m_componentMolarWeigth[ic];

      real64 const phaseCompMassDens = phaseCompMoleFrac.value[ip][ic] * mw / props->MolecularWeight;
      real64 const dPhaseCompMassDens_dPres =
        (phaseCompMoleFrac.dPres[ip][ic] * mw - phaseCompMassDens * dPhaseMolarWeight_dPres) / props->MolecularWeight;
      real64 const dPhaseCompMassDens_dTemp =
        (phaseCompMoleFrac.dTemp[ip][ic] * mw - phaseCompMassDens * dPhaseMolarWeight_dTemp) / props->MolecularWeight;

      phaseCompMoleFrac.value[ip][ic] = phaseCompMassDens * phaseDens.value[ip];
      phaseCompMoleFrac.dPres[ip][ic] =
        dPhaseCompMassDens_dPres * phaseDens.value[ip] + phaseCompMassDens * phaseDens.dPres[ip];
      phaseCompMoleFrac.dTemp[ip][ic] =
        dPhaseCompMassDens_dTemp * phaseDens.value[ip] + phaseCompMassDens * phaseDens.dTemp[ip];

      for (localIndex jc = 0; jc < numComp; ++jc)
      {
        real64 const dPhaseMolarWeight_dComp = 0.0; // TODO

        real64 const dPhaseCompMassDens_dComp =
          (phaseCompMoleFrac.dComp[ip][ic][jc] * mw - phaseCompMassDens * dPhaseMolarWeight_dComp) / props->MolecularWeight;

        phaseCompMoleFrac.dComp[ip][ic][jc] =
          dPhaseCompMassDens_dComp * phaseDens.value[ip] + phaseCompMassDens * phaseDens.dComp[ip][jc];
      }
    }
  }

#endif
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompositionalMultiphaseFluid, std::string const &, ManagedGroup * const )
} // namespace constitutive

} // namespace geosx