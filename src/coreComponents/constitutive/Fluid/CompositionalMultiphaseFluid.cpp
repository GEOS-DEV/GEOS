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

#include "codingUtilities/Utilities.hpp"

// PVTPackage includes
#include "MultiphaseSystem/CompositionalMultiphaseSystem.hpp"


using namespace PVTPackage;

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

namespace constitutive
{

namespace
{

std::unordered_map<string, PHASE_TYPE> const phaseNameDict =
  {
    { "gas",   PHASE_TYPE::GAS },
    { "oil",   PHASE_TYPE::OIL },
    { "water", PHASE_TYPE::LIQUID_WATER_RICH }
  };

std::unordered_map<string, EOS_TYPE> const eosNameDict =
  {
    { "PR",   EOS_TYPE::PENG_ROBINSON },
    { "SRK",  EOS_TYPE::REDLICH_KWONG_SOAVE }
  };

}

CompositionalMultiphaseFluid::CompositionalMultiphaseFluid( std::string const & name, ManagedGroup * const parent )
  : MultiFluidBase( name, parent ),
    m_fluid( nullptr )
{
  RegisterViewWrapper( viewKeys.equationsOfState.Key(), &m_equationsOfState, false );

  RegisterViewWrapper( viewKeys.componentCriticalPressure.Key(), &m_componentCriticalPressure, false );
  RegisterViewWrapper( viewKeys.componentCriticalTemperature.Key(), &m_componentCriticalTemperature, false );
  RegisterViewWrapper( viewKeys.componentAcentricFactor.Key(), &m_componentAcentricFactor, false );
  RegisterViewWrapper( viewKeys.componentMolarWeight.Key(), &m_componentMolarWeight, false );
  RegisterViewWrapper( viewKeys.componentVolumeShift.Key(), &m_componentVolumeShift, false );
  RegisterViewWrapper( viewKeys.componentBinaryCoeff.Key(), &m_componentBinaryCoeff, false );
}

CompositionalMultiphaseFluid::~CompositionalMultiphaseFluid()
{
  delete m_fluid;
}

std::unique_ptr<ConstitutiveBase>
CompositionalMultiphaseFluid::DeliverClone( string const & name, ManagedGroup * const parent ) const
{
  auto clone = std::make_unique<CompositionalMultiphaseFluid>( name, parent );

  clone->m_useMassFractions = this->m_useMassFractions;

  clone->m_phases           = this->m_phases;
  clone->m_equationsOfState = this->m_equationsOfState;
  clone->m_componentNames   = this->m_componentNames;

  clone->m_componentCriticalPressure    = this->m_componentCriticalPressure;
  clone->m_componentCriticalTemperature = this->m_componentCriticalTemperature;
  clone->m_componentAcentricFactor      = this->m_componentAcentricFactor;
  clone->m_componentMolarWeight         = this->m_componentMolarWeight;
  clone->m_componentVolumeShift         = this->m_componentVolumeShift;
  clone->m_componentBinaryCoeff         = this->m_componentBinaryCoeff;

  clone->createFluid();

  return clone;
}

void CompositionalMultiphaseFluid::FillDocumentationNode()
{
  MultiFluidBase::FillDocumentationNode();

  DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setShortDescription( "Compositional multiphase fluid model based on PVTPackage" );

  docNode->AllocateChildNode( viewKeys.equationsOfState.Key(),
                              viewKeys.equationsOfState.Key(),
                              -1,
                              "string_array",
                              "string_array",
                              "List of equation of state types for each phase",
                              "List of equation of state types for each phase",
                              "REQUIRED",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.componentCriticalPressure.Key(),
                              viewKeys.componentCriticalPressure.Key(),
                              -1,
                              "real64_array",
                              "real64_array",
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
                              "real64_array",
                              "real64_array",
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
                              "real64_array",
                              "real64_array",
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
                              "real64_array",
                              "real64_array",
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
                              "real64_array",
                              "real64_array",
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
  MultiFluidBase::ReadXML_PostProcess();

  localIndex const NC = numFluidComponents();
  localIndex const NP = numFluidPhases();

#define COMPFLUID_CHECK_INPUT_LENGTH( data, expected, attr ) \
  if (integer_conversion<localIndex>((data).size()) != integer_conversion<localIndex>(expected)) \
  { \
    GEOS_ERROR("CompositionalMultiphaseFluid: invalid number of entries in " \
               << (attr) << " attribute (" \
               << (data).size() << "given, " \
               << (expected) << " expected)"); \
  }

  COMPFLUID_CHECK_INPUT_LENGTH(m_equationsOfState, NP,
                               viewKeys.equationsOfState.Key())

  COMPFLUID_CHECK_INPUT_LENGTH(m_componentCriticalPressure, NC,
                               viewKeys.componentCriticalPressure.Key())

  COMPFLUID_CHECK_INPUT_LENGTH(m_componentCriticalTemperature, NC,
                               viewKeys.componentCriticalTemperature.Key())

  COMPFLUID_CHECK_INPUT_LENGTH(m_componentAcentricFactor, NC,
                               viewKeys.componentAcentricFactor.Key())

  COMPFLUID_CHECK_INPUT_LENGTH(m_componentMolarWeight, NC,
                               viewKeys.componentMolarWeight.Key())

  if (m_componentVolumeShift.empty())
  {
    m_componentVolumeShift.resize( NC );
    m_componentVolumeShift = 0.0;
  }

  COMPFLUID_CHECK_INPUT_LENGTH(m_componentVolumeShift, NC,
                               viewKeys.componentVolumeShift.Key())

  //if (m_componentBinaryCoeff.empty()) TODO
  {
    m_componentBinaryCoeff.resize( NC, NC );
    m_componentBinaryCoeff = 0.0;
  }

  COMPFLUID_CHECK_INPUT_LENGTH(m_componentBinaryCoeff, NC * NC,
                               viewKeys.componentBinaryCoeff.Key())

#undef CHECK_INPUT_LENGTH
}

void CompositionalMultiphaseFluid::createFluid()
{
  localIndex const numComp  = numFluidComponents();
  localIndex const numPhase = numFluidPhases();

  std::vector<PHASE_TYPE> phases(numPhase);
  std::vector<EOS_TYPE> eos(numPhase);

  for (localIndex ip = 0; ip < numPhase; ++ip)
  {
    if (phaseNameDict.find(m_phases[ip]) != phaseNameDict.end())
    {
      phases[ip] = phaseNameDict.at(m_phases[ip]);
    }
    else
    {
      GEOS_ERROR("CompositionalMultiphaseFluid: invalid phase label: " << m_phases[ip]);
    }

    if (eosNameDict.find(m_equationsOfState[ip]) != eosNameDict.end())
    {
      eos[ip] = eosNameDict.at(m_equationsOfState[ip]);
    }
    else
    {
      GEOS_ERROR("CompositionalMultiphaseFluid: invalid EOS label: " << m_equationsOfState[ip]);
    }
  }

  std::vector<std::string> components(m_componentNames.begin(), m_componentNames.end());
  std::vector<double> Pc(m_componentCriticalPressure.begin(), m_componentCriticalPressure.end());
  std::vector<double> Tc(m_componentCriticalTemperature.begin(), m_componentCriticalTemperature.end());
  std::vector<double> Mw(m_componentMolarWeight.begin(), m_componentMolarWeight.end());
  std::vector<double> Omega(m_componentAcentricFactor.begin(), m_componentAcentricFactor.end());

  const ComponentProperties CompProps(numComp, components, Mw, Tc, Pc, Omega);
  // TODO choose flash type
  m_fluid = new CompositionalMultiphaseSystem(phases, eos, COMPOSITIONAL_FLASH_TYPE::NEGATIVE_OIL_GAS, CompProps);

}

void CompositionalMultiphaseFluid::InitializePostSubGroups( ManagedGroup * const group )
{
  MultiFluidBase::InitializePostSubGroups( group );
  createFluid();
}

namespace
{

template<typename T, int DIM>
struct array_view_helper
{
  using type = array_view<T, DIM>;
};

#ifndef GEOSX_USE_ARRAY_BOUNDS_CHECK
// special tratment since there is no implicit conversion double * -> array_view<T, 1>
template<typename T>
struct array_view_helper<T, 1>
{
  using type = T *;
};
#endif

// an array view of DIM=0 decays to a reference to scalar
template<typename T>
struct array_view_helper<T, 0>
{
  using type = T &;
};

template<int DIM>
using real_array_view = typename array_view_helper<real64, DIM>::type;


// helper struct to represent a var and its derivatives
template<int DIM>
struct VarContainer
{
  real_array_view<DIM>   value; // variable value
  real_array_view<DIM>   dPres; // derivative w.r.t. pressure
  real_array_view<DIM>   dTemp; // derivative w.r.t. temperature
  real_array_view<DIM+1> dComp; // derivative w.r.t. composition
};

}

void CompositionalMultiphaseFluid::StateUpdatePointMultiphaseFluid( real64 const & pres,
                                                                    real64 const & temp,
                                                                    real64 const * composition,
                                                                    localIndex const k,
                                                                    localIndex const q )
{
  // 0. set array views to the element/point data to avoid awkward quadruple indexing
  VarContainer<1> phaseFrac {
    m_phaseFraction[k][q],
    m_dPhaseFraction_dPressure[k][q],
    m_dPhaseFraction_dTemperature[k][q],
    m_dPhaseFraction_dGlobalCompFraction[k][q]
  };

  VarContainer<1> phaseDens {
    m_phaseDensity[k][q],
    m_dPhaseDensity_dPressure[k][q],
    m_dPhaseDensity_dTemperature[k][q],
    m_dPhaseDensity_dGlobalCompFraction[k][q]
  };

  VarContainer<1> phaseVolFrac {
    m_phaseVolumeFraction[k][q],
    m_dPhaseVolumeFraction_dPressure[k][q],
    m_dPhaseVolumeFraction_dTemperature[k][q],
    m_dPhaseVolumeFraction_dGlobalCompFraction[k][q]
  };

  VarContainer<2> phaseCompFrac {
    m_phaseCompFraction[k][q],
    m_dPhaseCompFraction_dPressure[k][q],
    m_dPhaseCompFraction_dTemperature[k][q],
    m_dPhaseCompFraction_dGlobalCompFraction[k][q]
  };

  VarContainer<0> totalDens {
    m_totalDensity[k][q],
    m_dTotalDensity_dPressure[k][q],
    m_dTotalDensity_dTemperature[k][q],
    m_dTotalDensity_dGlobalCompFraction[k][q]
  };

  localIndex const NP = numFluidPhases();
  localIndex const NC = numFluidComponents();

  // 1. Convert input mass fractions to mole fractions and keep derivatives

  std::vector<double> compMoleFrac( NC );
  array2d<real64> dCompMoleFrac_dCompMassFrac;

  if (m_useMassFractions)
  {
    dCompMoleFrac_dCompMassFrac.resize( NC, NC );
    dCompMoleFrac_dCompMassFrac = 0.0;

    real64 totalMolality = 0.0;
    for (localIndex ic = 0; ic < NC; ++ic)
    {
      compMoleFrac[ic] = composition[ic] / m_componentMolarWeight[ic]; // this is molality (units of mol/kg)
      dCompMoleFrac_dCompMassFrac[ic][ic] = 1.0 / m_componentMolarWeight[ic];
      totalMolality += compMoleFrac[ic];
    }

    for (localIndex ic = 0; ic < NC; ++ic)
    {
      compMoleFrac[ic] /= totalMolality;

      for (localIndex jc = 0; jc < NC; ++jc)
      {
        dCompMoleFrac_dCompMassFrac[ic][jc] = -compMoleFrac[ic] / m_componentMolarWeight[jc];
        dCompMoleFrac_dCompMassFrac[ic][jc] /= totalMolality;
      }
    }
  }
  else
  {
    for (localIndex ic = 0; ic < NC; ++ic)
    {
      compMoleFrac[ic] = composition[ic];
    }
  }

  // 2. Trigger PVTPackage compute and get back phase split

  m_fluid->Update(pres, temp, compMoleFrac);
  MultiphaseSystemProperties const & split = m_fluid->get_MultiphaseSystemProperties();

  // 3. Extract phase split, phase properties and derivatives from PVTPackage

  for (localIndex ip = 0; ip < NP; ++ip)
  {
    PHASE_TYPE phase = phaseNameDict.at(m_phases[ip]);
    PhaseProperties const & props = m_fluid->get_PhaseProperties(phase);

    auto const & frac = split.PhaseMoleFraction.at(phase);
    auto const & comp = props.MoleComposition;
    auto const & dens = m_useMassFractions ? props.MassDensity : props.MoleDensity;

    phaseFrac.value[ip] = frac.value;
    phaseFrac.dPres[ip] = frac.dP;
    phaseFrac.dTemp[ip] = frac.dT;

    phaseDens.value[ip] = dens.value;
    phaseDens.dPres[ip] = dens.dP;
    phaseDens.dTemp[ip] = dens.dT;

    for (localIndex ic = 0; ic < NC; ++ic)
    {
      phaseFrac.dComp[ip][ic] = frac.dz[ic];
      phaseDens.dComp[ip][ic] = dens.dz[ic];

      phaseCompFrac.value[ip][ic] = comp.value[ic];
      phaseCompFrac.dPres[ip][ic] = comp.dP[ic];
      phaseCompFrac.dTemp[ip][ic] = comp.dT[ic];

      for (localIndex jc = 0; jc < NC; ++jc)
        phaseCompFrac.dComp[ip][ic][jc] = comp.dz[ic][jc];
    }
  }

  // 4. Compute phase volume fractions and derivatives, requires two passes for normalization

  real64 totalMolDensInv = 0.0;
  real64 dTotalMolDensInv_dP = 0.0;
  real64 dTotalMolDensInv_dT = 0.0;
  real64 dTotalMolDensInv_dC[MAX_NUM_COMPONENTS]{};

  for (localIndex ip = 0; ip < NP; ++ip)
  {
    PHASE_TYPE phase = phaseNameDict.at(m_phases[ip]);
    PhaseProperties const & props = m_fluid->get_PhaseProperties(phase);

    auto const & phaseMolarDens = props.MoleDensity;

    phaseVolFrac.value[ip] = phaseFrac.value[ip] / phaseMolarDens.value;

    phaseVolFrac.dPres[ip] =
      (phaseFrac.dPres[ip] - phaseVolFrac.value[ip] * phaseMolarDens.dP) / phaseMolarDens.value;

    phaseVolFrac.dTemp[ip] =
      (phaseFrac.dTemp[ip] - phaseVolFrac.value[ip] * phaseMolarDens.dT) / phaseMolarDens.value;

    totalMolDensInv += phaseVolFrac.value[ip];
    dTotalMolDensInv_dP += phaseVolFrac.dPres[ip];
    dTotalMolDensInv_dT += phaseVolFrac.dTemp[ip];

    for (localIndex jc = 0; jc < NC; ++jc)
    {
      phaseVolFrac.dComp[ip][jc] =
        (phaseFrac.dComp[ip][jc] - phaseVolFrac.value[ip] * phaseMolarDens.dz[jc]) / phaseMolarDens.value;

      dTotalMolDensInv_dC[jc] += phaseVolFrac.dComp[ip][jc];
    }
  }
  // normalization pass
  for (localIndex ip = 0; ip < NP; ++ip)
  {
    phaseVolFrac.value[ip] /= totalMolDensInv;
    phaseVolFrac.dPres[ip] = (phaseVolFrac.dPres[ip] - phaseVolFrac.value[ip] * dTotalMolDensInv_dP) / totalMolDensInv;
    phaseVolFrac.dTemp[ip] = (phaseVolFrac.dTemp[ip] - phaseVolFrac.value[ip] * dTotalMolDensInv_dT) / totalMolDensInv;

    for (localIndex ic = 0; ic < NC; ++ic)
    {
      phaseVolFrac.dComp[ip][ic] -= phaseVolFrac.value[ip] * dTotalMolDensInv_dC[ic];
      phaseVolFrac.dComp[ip][ic] /= totalMolDensInv;
    }
  }

  // 5. Compute total fluid mass/molar density
  totalDens.value = 0.0;
  totalDens.dPres = 0.0;
  totalDens.dTemp = 0.0;
  for (localIndex jc = 0; jc < NC; ++jc)
  {
    totalDens.dComp[jc] = 0.0;
  }
  for (localIndex ip = 0; ip < NP; ++ip)
  {
    totalDens.value += phaseDens.value[ip] * phaseVolFrac.value[ip];
    totalDens.dPres += phaseDens.dPres[ip] * phaseVolFrac.value[ip] + phaseDens.value[ip] * phaseVolFrac.dPres[ip];
    totalDens.dTemp += phaseDens.dTemp[ip] * phaseVolFrac.value[ip] + phaseDens.value[ip] * phaseVolFrac.dTemp[ip];

    for (localIndex jc = 0; jc < NC; ++jc)
    {
      totalDens.dComp[jc] += phaseDens.dComp[ip][jc] * phaseVolFrac.value[ip] + phaseDens.value[ip] * phaseVolFrac.dComp[ip][jc];
    }
  }

  if (m_useMassFractions)
  {
    // 6. Convert output mole fractions to mass fractions
    // At this point, phase and total densities are already mass densities

    // 6.1. Compute total molar weight and derivatives
    real64 dTotalMW_dC[MAX_NUM_COMPONENTS]{};

    real64 const totalMW = totalDens.value * totalMolDensInv;
    real64 const dTotalMW_dP = totalDens.dPres * totalMolDensInv + totalDens.value * dTotalMolDensInv_dP;
    real64 const dTotalMW_dT = totalDens.dTemp * totalMolDensInv + totalDens.value * dTotalMolDensInv_dT;

    for (localIndex jc = 0; jc < NC; ++jc)
    {
      dTotalMW_dC[jc] = totalDens.dComp[jc] * totalMolDensInv + totalDens.value * dTotalMolDensInv_dC[jc];
    }

    // 6.2. Finally, convert fractions
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      PHASE_TYPE phase = phaseNameDict.at(m_phases[ip]);
      PhaseProperties const & props = m_fluid->get_PhaseProperties(phase);

      auto const & phaseMW = props.MolecularWeight;

      phaseFrac.value[ip] *= phaseMW.value / totalMW;
      phaseFrac.dPres[ip] = phaseFrac.dPres[ip] * phaseMW.value / totalMW
                          + phaseFrac.value[ip] * (phaseMW.dP / phaseMW.value - dTotalMW_dP / totalMW);
      phaseFrac.dTemp[ip] = phaseFrac.dTemp[ip] * phaseMW.value / totalMW
                            + phaseFrac.value[ip] * (phaseMW.dT / phaseMW.value - dTotalMW_dT / totalMW);

      for (localIndex jc = 0; jc < NC; ++jc)
      {
        phaseFrac.dComp[ip][jc] = phaseFrac.dComp[ip][jc] * phaseMW.value / totalMW
                                + phaseFrac.value[ip] * (phaseMW.dz[jc] / phaseMW.value - dTotalMW_dC[jc] / totalMW);
      }

      for (localIndex ic = 0; ic < NC; ++ic)
      {

        real64 const compMW = m_componentMolarWeight[ic];

        phaseCompFrac.value[ip][ic] = phaseCompFrac.value[ip][ic] * compMW / phaseMW.value;
        phaseCompFrac.dPres[ip][ic] =
          (phaseCompFrac.dPres[ip][ic] * compMW - phaseCompFrac.value[ip][ic] * phaseMW.dP) / phaseMW.value;
        phaseCompFrac.dTemp[ip][ic] =
          (phaseCompFrac.dTemp[ip][ic] * compMW - phaseCompFrac.value[ip][ic] * phaseMW.dT) / phaseMW.value;

        for (localIndex jc = 0; jc < NC; ++jc)
        {
          phaseCompFrac.dComp[ip][ic][jc] =
            (phaseCompFrac.dComp[ip][ic][jc] * compMW - phaseCompFrac.value[ip][ic] * phaseMW.dz[jc]) / phaseMW.value;
        }
      }
    }

    // 7. Update derivatives w.r.t. mole fractions to derivatives w.r.t mass fractions

    array1d<real64> work( NC );
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      applyChainRuleInPlace( NC, dCompMoleFrac_dCompMassFrac, phaseFrac.dComp[ip], work );
      applyChainRuleInPlace( NC, dCompMoleFrac_dCompMassFrac, phaseDens.dComp[ip], work );
      applyChainRuleInPlace( NC, dCompMoleFrac_dCompMassFrac, phaseVolFrac.dComp[ip], work );

      for (localIndex ic = 0; ic < NC; ++ic)
      {
        applyChainRuleInPlace( NC, dCompMoleFrac_dCompMassFrac, phaseCompFrac.dComp[ip][ic], work );
      }
    }
    applyChainRuleInPlace( NC, dCompMoleFrac_dCompMassFrac, totalDens.dComp, work );
  }
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompositionalMultiphaseFluid, std::string const &, ManagedGroup * const )
} // namespace constitutive

} // namespace geosx
