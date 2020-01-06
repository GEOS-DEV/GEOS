/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
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

CompositionalMultiphaseFluid::CompositionalMultiphaseFluid( std::string const & name, Group * const parent )
  : MultiFluidPVTPackageWrapper( name, parent )
{
  getWrapperBase( viewKeyStruct::componentNamesString )->setInputFlag(InputFlags::REQUIRED);
  getWrapperBase( viewKeyStruct::componentMolarWeightString )->setInputFlag(InputFlags::REQUIRED);
  getWrapperBase( viewKeyStruct::phaseNamesString )->setInputFlag(InputFlags::REQUIRED);

  registerWrapper( viewKeyStruct::equationsOfStateString, &m_equationsOfState, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("List of equation of state types for each phase");

  registerWrapper( viewKeyStruct::componentCriticalPressureString, &m_componentCriticalPressure, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Component critical pressures");

  registerWrapper( viewKeyStruct::componentCriticalTemperatureString, &m_componentCriticalTemperature, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Component critical temperatures");

  registerWrapper( viewKeyStruct::componentAcentricFactorString, &m_componentAcentricFactor, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Component acentric factors");

  registerWrapper( viewKeyStruct::componentVolumeShiftString, &m_componentVolumeShift, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Component volume shifts");

  registerWrapper( viewKeyStruct::componentBinaryCoeffString, &m_componentBinaryCoeff, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Table of binary interaction coefficients");
}

CompositionalMultiphaseFluid::~CompositionalMultiphaseFluid()
{

}

void
CompositionalMultiphaseFluid::DeliverClone( string const & name,
                                            Group * const parent,
                                            std::unique_ptr<ConstitutiveBase> & clone ) const
{
  std::unique_ptr< CompositionalMultiphaseFluid > newModel = std::make_unique<CompositionalMultiphaseFluid>( name, parent );

  newModel->m_useMass = this->m_useMass;

  newModel->m_componentNames   = this->m_componentNames;
  newModel->m_componentMolarWeight = this->m_componentMolarWeight;

  newModel->m_phaseNames           = this->m_phaseNames;
  newModel->m_pvtPackagePhaseTypes = this->m_pvtPackagePhaseTypes;
  newModel->m_equationsOfState     = this->m_equationsOfState;

  newModel->m_componentCriticalPressure    = this->m_componentCriticalPressure;
  newModel->m_componentCriticalTemperature = this->m_componentCriticalTemperature;
  newModel->m_componentAcentricFactor      = this->m_componentAcentricFactor;
  newModel->m_componentVolumeShift         = this->m_componentVolumeShift;
  newModel->m_componentBinaryCoeff         = this->m_componentBinaryCoeff;

  newModel->createFluid();

  clone = std::move( newModel );
}

void CompositionalMultiphaseFluid::PostProcessInput()
{
  MultiFluidPVTPackageWrapper::PostProcessInput();

  localIndex const NC = numFluidComponents();
  localIndex const NP = numFluidPhases();

#define COMPFLUID_CHECK_INPUT_LENGTH( data, expected, attr ) \
  if (integer_conversion<localIndex>((data).size()) != integer_conversion<localIndex>(expected)) \
  { \
    GEOSX_ERROR( "CompositionalMultiphaseFluid: invalid number of entries in " \
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

  //if (m_componentBinaryCoeff.empty()) TODO needs reading of 2D arrays
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
    GEOSX_ERROR_IF( it == PVTPackage_eosDict.end(), "Invalid eos name: " << m_equationsOfState[ip] );
    eos[ip] = it->second;
  }

  std::vector<PHASE_TYPE> const phases( m_pvtPackagePhaseTypes.begin(), m_pvtPackagePhaseTypes.end() );
  std::vector<std::string> const components( m_componentNames.begin(), m_componentNames.end() );
  std::vector<double> const Pc( m_componentCriticalPressure.begin(), m_componentCriticalPressure.end() );
  std::vector<double> const Tc( m_componentCriticalTemperature.begin(), m_componentCriticalTemperature.end() );
  std::vector<double> const Mw( m_componentMolarWeight.begin(), m_componentMolarWeight.end() );
  std::vector<double> const Omega( m_componentAcentricFactor.begin(), m_componentAcentricFactor.end() );

  ComponentProperties const compProps( NC, components, Mw, Tc, Pc, Omega );

  m_fluid = std::make_unique<CompositionalMultiphaseSystem>( phases, eos,
                                                             COMPOSITIONAL_FLASH_TYPE::NEGATIVE_OIL_GAS,
                                                             compProps );

}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompositionalMultiphaseFluid, std::string const &, Group * const )
} // namespace constitutive

} // namespace geosx
