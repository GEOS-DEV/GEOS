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

std::unordered_map<string, EOS_TYPE> const PVTPackage_eosDict =
{
  { "PR",   EOS_TYPE::PENG_ROBINSON },
  { "SRK",  EOS_TYPE::REDLICH_KWONG_SOAVE }
};

}

CompositionalMultiphaseFluid::CompositionalMultiphaseFluid( std::string const & name, ManagedGroup * const parent )
  : MultiFluidPVTPackageWrapper( name, parent )
{
  RegisterViewWrapper( viewKeyStruct::equationsOfStateString, &m_equationsOfState, false );

  RegisterViewWrapper( viewKeyStruct::componentCriticalPressureString, &m_componentCriticalPressure, false );
  RegisterViewWrapper( viewKeyStruct::componentCriticalTemperatureString, &m_componentCriticalTemperature, false );
  RegisterViewWrapper( viewKeyStruct::componentAcentricFactorString, &m_componentAcentricFactor, false );
  RegisterViewWrapper( viewKeyStruct::componentVolumeShiftString, &m_componentVolumeShift, false );
  RegisterViewWrapper( viewKeyStruct::componentBinaryCoeffString, &m_componentBinaryCoeff, false );
}

CompositionalMultiphaseFluid::~CompositionalMultiphaseFluid()
{

}

std::unique_ptr<ConstitutiveBase>
CompositionalMultiphaseFluid::DeliverClone( string const & name, ManagedGroup * const parent ) const
{
  std::unique_ptr< CompositionalMultiphaseFluid > clone = std::make_unique<CompositionalMultiphaseFluid>( name, parent );

  clone->m_useMass = this->m_useMass;

  clone->m_componentNames   = this->m_componentNames;
  clone->m_componentMolarWeight = this->m_componentMolarWeight;

  clone->m_phaseNames           = this->m_phaseNames;
  clone->m_pvtPackagePhaseTypes = this->m_pvtPackagePhaseTypes;
  clone->m_equationsOfState     = this->m_equationsOfState;

  clone->m_componentCriticalPressure    = this->m_componentCriticalPressure;
  clone->m_componentCriticalTemperature = this->m_componentCriticalTemperature;
  clone->m_componentAcentricFactor      = this->m_componentAcentricFactor;
  clone->m_componentVolumeShift         = this->m_componentVolumeShift;
  clone->m_componentBinaryCoeff         = this->m_componentBinaryCoeff;

  clone->createFluid();

  std::unique_ptr<ConstitutiveBase> rval = std::move( clone );
  return rval;
}

void CompositionalMultiphaseFluid::FillDocumentationNode()
{
  MultiFluidPVTPackageWrapper::FillDocumentationNode();

  DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setShortDescription( "Compositional multiphase fluid model based on PVTPackage" );

  docNode->AllocateChildNode( viewKeyStruct::equationsOfStateString,
                              viewKeyStruct::equationsOfStateString,
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

  docNode->AllocateChildNode( viewKeyStruct::componentCriticalPressureString,
                              viewKeyStruct::componentCriticalPressureString,
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

  docNode->AllocateChildNode( viewKeyStruct::componentCriticalTemperatureString,
                              viewKeyStruct::componentCriticalTemperatureString,
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

  docNode->AllocateChildNode( viewKeyStruct::componentAcentricFactorString,
                              viewKeyStruct::componentAcentricFactorString,
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

  docNode->AllocateChildNode( viewKeyStruct::componentVolumeShiftString,
                              viewKeyStruct::componentVolumeShiftString,
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
  MultiFluidPVTPackageWrapper::ReadXML_PostProcess();

  localIndex const NC = numFluidComponents();
  localIndex const NP = numFluidPhases();

#define COMPFLUID_CHECK_INPUT_LENGTH( data, expected, attr ) \
  if (integer_conversion<localIndex>((data).size()) != integer_conversion<localIndex>(expected)) \
  { \
    GEOS_ERROR( "CompositionalMultiphaseFluid: invalid number of entries in " \
                << (attr) << " attribute (" \
                << (data).size() << "given, " \
                << (expected) << " expected)"); \
  }

  COMPFLUID_CHECK_INPUT_LENGTH( m_equationsOfState, NP, viewKeyStruct::equationsOfStateString )
  COMPFLUID_CHECK_INPUT_LENGTH( m_componentCriticalPressure, NC, viewKeyStruct::componentCriticalPressureString )
  COMPFLUID_CHECK_INPUT_LENGTH( m_componentCriticalTemperature, NC, viewKeyStruct::componentCriticalTemperatureString )
  COMPFLUID_CHECK_INPUT_LENGTH( m_componentAcentricFactor, NC, viewKeyStruct::componentAcentricFactorString )

  if (m_componentVolumeShift.empty())
  {
    m_componentVolumeShift.resize( NC );
    m_componentVolumeShift = 0.0;
  }

  COMPFLUID_CHECK_INPUT_LENGTH( m_componentVolumeShift, NC, viewKeyStruct::componentVolumeShiftString )

  //if (m_componentBinaryCoeff.empty()) TODO
  {
    m_componentBinaryCoeff.resize( NC, NC );
    m_componentBinaryCoeff = 0.0;
  }

  COMPFLUID_CHECK_INPUT_LENGTH( m_componentBinaryCoeff, NC * NC, viewKeyStruct::componentBinaryCoeffString )

#undef COMPFLUID_CHECK_INPUT_LENGTH
}

void CompositionalMultiphaseFluid::createFluid()
{
  localIndex const NC = numFluidComponents();
  localIndex const NP = numFluidPhases();

  std::vector<EOS_TYPE> eos( NP );

  for (localIndex ip = 0; ip < NP; ++ip)
  {
    auto it = PVTPackage_eosDict.find( m_equationsOfState[ip] );
    GEOS_ERROR_IF( it == PVTPackage_eosDict.end(), "Invalid eos name: " << m_equationsOfState[ip] );
    eos[ip] = it->second;
  }

  std::vector<PHASE_TYPE> phases( m_pvtPackagePhaseTypes.begin(), m_pvtPackagePhaseTypes.end() );
  std::vector<std::string> components( m_componentNames.begin(), m_componentNames.end() );
  std::vector<double> Pc( m_componentCriticalPressure.begin(), m_componentCriticalPressure.end() );
  std::vector<double> Tc( m_componentCriticalTemperature.begin(), m_componentCriticalTemperature.end() );
  std::vector<double> Mw( m_componentMolarWeight.begin(), m_componentMolarWeight.end() );
  std::vector<double> Omega( m_componentAcentricFactor.begin(), m_componentAcentricFactor.end() );

  const ComponentProperties CompProps( NC, components, Mw, Tc, Pc, Omega );
  // TODO choose flash type
  m_fluid = new CompositionalMultiphaseSystem( phases, eos, COMPOSITIONAL_FLASH_TYPE::NEGATIVE_OIL_GAS, CompProps );

}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompositionalMultiphaseFluid, std::string const &, ManagedGroup * const )
} // namespace constitutive

} // namespace geosx
