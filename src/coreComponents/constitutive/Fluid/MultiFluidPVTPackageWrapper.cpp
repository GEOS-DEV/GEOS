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
  * @file MultiFluidPVTPackageWrapper.cpp
  */

#include "MultiFluidPVTPackageWrapper.hpp"

// PVTPackage includes
#include "MultiphaseSystem/MultiphaseSystem.hpp"

using namespace PVTPackage;

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

namespace
{

std::unordered_map<string, PHASE_TYPE> const PVTPackage_phaseDict =
{
  { "gas",   PHASE_TYPE::GAS },
  { "oil",   PHASE_TYPE::OIL },
  { "water", PHASE_TYPE::LIQUID_WATER_RICH }
};

}

MultiFluidPVTPackageWrapper::MultiFluidPVTPackageWrapper( std::string const & name, ManagedGroup * const parent )
  : MultiFluidBase( name, parent ),
    m_fluid( nullptr )
{

}

MultiFluidPVTPackageWrapper::~MultiFluidPVTPackageWrapper()
{
  delete m_fluid;
}

void MultiFluidPVTPackageWrapper::ProcessInputFile_PostProcess()
{
  MultiFluidBase::ProcessInputFile_PostProcess();

  localIndex const NP = numFluidPhases();

  m_pvtPackagePhaseTypes.resize( NP );

  for (localIndex ip = 0; ip < NP; ++ip)
  {
    auto it = PVTPackage_phaseDict.find( m_phaseNames[ip] );
    GEOS_ERROR_IF( it == PVTPackage_phaseDict.end(), "Phase not supported by PVTPackage: " << m_phaseNames[ip] );
    m_pvtPackagePhaseTypes[ip] = it->second;
  }
}

void MultiFluidPVTPackageWrapper::InitializePostSubGroups( ManagedGroup * const group )
{
  MultiFluidBase::InitializePostSubGroups(group);
  createFluid();
}

void MultiFluidPVTPackageWrapper::PointUpdate(real64 const & pres,
                                              real64 const & temp,
                                              arraySlice1d<real64 const> const & composition,
                                              localIndex const k,
                                              localIndex const q)
{
  // 0. set array views to the element/point data to avoid awkward quadruple indexing
  CompositionalVarContainer<1> phaseFrac {
    m_phaseFraction[k][q],
    m_dPhaseFraction_dPressure[k][q],
    m_dPhaseFraction_dTemperature[k][q],
    m_dPhaseFraction_dGlobalCompFraction[k][q]
  };

  CompositionalVarContainer<1> phaseDens {
    m_phaseDensity[k][q],
    m_dPhaseDensity_dPressure[k][q],
    m_dPhaseDensity_dTemperature[k][q],
    m_dPhaseDensity_dGlobalCompFraction[k][q]
  };

  CompositionalVarContainer<1> phaseVisc {
    m_phaseViscosity[k][q],
    m_dPhaseViscosity_dPressure[k][q],
    m_dPhaseViscosity_dTemperature[k][q],
    m_dPhaseViscosity_dGlobalCompFraction[k][q]
  };

  CompositionalVarContainer<2> phaseCompFrac {
    m_phaseCompFraction[k][q],
    m_dPhaseCompFraction_dPressure[k][q],
    m_dPhaseCompFraction_dTemperature[k][q],
    m_dPhaseCompFraction_dGlobalCompFraction[k][q]
  };

  CompositionalVarContainer<0> totalDens {
    m_totalDensity[k][q],
    m_dTotalDensity_dPressure[k][q],
    m_dTotalDensity_dTemperature[k][q],
    m_dTotalDensity_dGlobalCompFraction[k][q]
  };

  localIndex constexpr maxNumComp = MAX_NUM_COMPONENTS;

  localIndex const NP = numFluidPhases();
  localIndex const NC = numFluidComponents();

  // 1. Convert input mass fractions to mole fractions and keep derivatives

  std::vector<double> compMoleFrac( NC );
  stackArray2d<real64, maxNumComp * maxNumComp> dCompMoleFrac_dCompMassFrac( NC, NC );

  if (m_useMass)
  {
    dCompMoleFrac_dCompMassFrac.resize( NC, NC );
    dCompMoleFrac_dCompMassFrac = 0.0;

    real64 totalMolality = 0.0;
    for (localIndex ic = 0; ic < NC; ++ic)
    {
      real64 const mwInv = 1.0 / m_componentMolarWeight[ic];
      compMoleFrac[ic] = composition[ic] * mwInv; // this is molality (units of mole/mass)
      dCompMoleFrac_dCompMassFrac[ic][ic] = mwInv;
      totalMolality += compMoleFrac[ic];
    }

    real64 const totalMolalityInv = 1.0 / totalMolality;
    for (localIndex ic = 0; ic < NC; ++ic)
    {
      compMoleFrac[ic] *= totalMolalityInv;

      for (localIndex jc = 0; jc < NC; ++jc)
      {
        dCompMoleFrac_dCompMassFrac[ic][jc] -= compMoleFrac[ic] / m_componentMolarWeight[jc];
        dCompMoleFrac_dCompMassFrac[ic][jc] *= totalMolalityInv;
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

  GEOS_WARNING_IF( m_fluid->getState() != MultiphaseSystem::State::SUCCESS,
                   "Phase equilibrium calculations not converged at constitutive point (" << k << "," << q << ")");

  MultiphaseSystemProperties const & split = m_fluid->get_MultiphaseSystemProperties();

  // 3. Extract phase split, phase properties and derivatives from PVTPackage
  for (localIndex ip = 0; ip < NP; ++ip)
  {
    PHASE_TYPE phase = PVTPackage_phaseDict.at(m_phaseNames[ip]);
    PhaseProperties const & props = m_fluid->get_PhaseProperties(phase);

    auto const & frac = split.PhaseMoleFraction.at(phase);
    auto const & comp = props.MoleComposition;
    auto const & dens = m_useMass ? props.MassDensity : props.MoleDensity;

    phaseFrac.value[ip] = frac.value;
    phaseFrac.dPres[ip] = frac.dP;
    phaseFrac.dTemp[ip] = frac.dT;

    phaseDens.value[ip] = dens.value;
    phaseDens.dPres[ip] = dens.dP;
    phaseDens.dTemp[ip] = dens.dT;

    // TODO
    phaseVisc.value[ip] = 1.0;
    phaseVisc.dPres[ip] = 0.0;
    phaseVisc.dTemp[ip] = 0.0;

    for (localIndex jc = 0; jc < NC; ++jc)
    {
      phaseFrac.dComp[ip][jc] = frac.dz[jc];
      phaseDens.dComp[ip][jc] = dens.dz[jc];
      phaseVisc.dComp[ip][jc] = 0.0; // TODO

      phaseCompFrac.value[ip][jc] = comp.value[jc];
      phaseCompFrac.dPres[ip][jc] = comp.dP[jc];
      phaseCompFrac.dTemp[ip][jc] = comp.dT[jc];

      for (localIndex ic = 0; ic < NC; ++ic)
        phaseCompFrac.dComp[ip][ic][jc] = comp.dz[ic][jc];
    }
  }

  // 4. if mass variables used instead of molar, perform the conversion
  if (m_useMass)
  {
    // 4.1. Convert phase fractions (requires two passes)
    real64 totalMass = 0.0;
    real64 dTotalMass_dP = 0.0;
    real64 dTotalMass_dT = 0.0;
    real64 dTotalMass_dC[MAX_NUM_COMPONENTS]{};

    // 4.1.1. Compute mass of each phase and total mass (on a 1-mole basis)
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      PHASE_TYPE phase = PVTPackage_phaseDict.at(m_phaseNames[ip]);
      PhaseProperties const & props = m_fluid->get_PhaseProperties(phase);

      auto const & phaseMW = props.MolecularWeight;
      real64 const nu = phaseFrac.value[ip];

      phaseFrac.value[ip] *= phaseMW.value;
      phaseFrac.dPres[ip] = phaseFrac.dPres[ip] * phaseMW.value + nu * phaseMW.dP;
      phaseFrac.dTemp[ip] = phaseFrac.dTemp[ip] * phaseMW.value + nu * phaseMW.dT;

      totalMass += phaseFrac.value[ip];
      dTotalMass_dP += phaseFrac.dPres[ip];
      dTotalMass_dT += phaseFrac.dTemp[ip];

      for (localIndex jc = 0; jc < NC; ++jc)
      {
        phaseFrac.dComp[ip][jc] = phaseFrac.dComp[ip][jc] * phaseMW.value + nu * phaseMW.dz[jc];
        dTotalMass_dC[jc] += phaseFrac.dComp[ip][jc];
      }
    }

    // 4.1.2. Normalize to get mass fractions
    real64 const totalMassInv = 1.0 / totalMass;
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      phaseFrac.value[ip] *= totalMassInv;
      phaseFrac.dPres[ip] = (phaseFrac.dPres[ip] - phaseFrac.value[ip] * dTotalMass_dP) * totalMassInv;
      phaseFrac.dTemp[ip] = (phaseFrac.dTemp[ip] - phaseFrac.value[ip] * dTotalMass_dT) * totalMassInv;

      for (localIndex jc = 0; jc < NC; ++jc)
      {
        phaseFrac.dComp[ip][jc] = (phaseFrac.dComp[ip][jc] - phaseFrac.value[ip] * dTotalMass_dC[jc]) * totalMassInv;
      }
    }

    // 4.2. Convert phase compositions
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      PHASE_TYPE phase = PVTPackage_phaseDict.at(m_phaseNames[ip]);
      PhaseProperties const & props = m_fluid->get_PhaseProperties(phase);

      auto const & phaseMW = props.MolecularWeight;
      real64 const phaseMWInv = 1.0 / phaseMW.value;

      for (localIndex ic = 0; ic < NC; ++ic)
      {

        real64 const compMW = m_componentMolarWeight[ic];

        phaseCompFrac.value[ip][ic] = phaseCompFrac.value[ip][ic] * compMW * phaseMWInv;
        phaseCompFrac.dPres[ip][ic] =
          (phaseCompFrac.dPres[ip][ic] * compMW - phaseCompFrac.value[ip][ic] * phaseMW.dP) * phaseMWInv;
        phaseCompFrac.dTemp[ip][ic] =
          (phaseCompFrac.dTemp[ip][ic] * compMW - phaseCompFrac.value[ip][ic] * phaseMW.dT) * phaseMWInv;

        for (localIndex jc = 0; jc < NC; ++jc)
        {
          phaseCompFrac.dComp[ip][ic][jc] =
            (phaseCompFrac.dComp[ip][ic][jc] * compMW - phaseCompFrac.value[ip][ic] * phaseMW.dz[jc]) * phaseMWInv;
        }
      }
    }

    // 4.3. Update derivatives w.r.t. mole fractions to derivatives w.r.t mass fractions
    array1d<real64> work( NC );
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      applyChainRuleInPlace( NC, dCompMoleFrac_dCompMassFrac, phaseFrac.dComp[ip], work );
      applyChainRuleInPlace( NC, dCompMoleFrac_dCompMassFrac, phaseDens.dComp[ip], work );

      for (localIndex ic = 0; ic < NC; ++ic)
      {
        applyChainRuleInPlace( NC, dCompMoleFrac_dCompMassFrac, phaseCompFrac.dComp[ip][ic], work );
      }
    }
  }

  // 5. Compute total fluid mass/molar density and derivatives
  {
    totalDens.value = 0.0;
    totalDens.dPres = 0.0;
    totalDens.dTemp = 0.0;
    for (localIndex jc = 0; jc < NC; ++jc)
    {
      totalDens.dComp[jc] = 0.0;
    }

    // 5.1. Sum mass/molar fraction/density ratio over all phases to get the inverse of density
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      real64 const densInv = 1.0 / phaseDens.value[ip];
      real64 const value = phaseFrac.value[ip] * densInv;

      totalDens.value += value;
      totalDens.dPres += (phaseFrac.dPres[ip] - value * phaseDens.dPres[ip]) * densInv;
      totalDens.dTemp += (phaseFrac.dTemp[ip] - value * phaseDens.dTemp[ip]) * densInv;

      for (localIndex jc = 0; jc < NC; ++jc)
      {
        totalDens.dComp[jc] += (phaseFrac.dComp[ip][jc] - value * phaseDens.dComp[ip][jc]) * densInv;
      }
    }

    // 5.2. Invert the previous quantity to get actual density
    totalDens.value = 1.0 / totalDens.value;
    real64 const minusDens2 = -totalDens.value * totalDens.value;
    totalDens.dPres *= minusDens2;
    totalDens.dTemp *= minusDens2;

    for (localIndex jc = 0; jc < NC; ++jc)
    {
      totalDens.dComp[jc] *= minusDens2;
    }
  }
}

} //namespace constitutive

} //namespace geosx
