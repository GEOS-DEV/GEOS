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
  * @file BlackOilFluid.cpp
  */

#include "BlackOilFluid.hpp"

#include "codingUtilities/Utilities.hpp"
#include "managers/ProblemManager.hpp"
#include "common/Path.hpp"

// PVTPackage includes
#include "MultiphaseSystem/BlackOilMultiphaseSystem.hpp"
#include "MultiphaseSystem/DeadOilMultiphaseSystem.hpp"


using namespace PVTPackage;

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

namespace constitutive
{

BlackOilFluid::FluidType BlackOilFluid::stringToFluidType( string const & str )
{
  if (str == "LiveOil")
  {
    return BlackOilFluid::FluidType::LiveOil;
  }
  else if (str == "DeadOil")
  {
    return BlackOilFluid::FluidType::DeadOil;
  }
  else
  {
    GEOSX_ERROR("Unrecognized black-oil fluid type: " << str);
  }
  return BlackOilFluid::FluidType::LiveOil; // keep compilers happy
}

BlackOilFluid::BlackOilFluid( std::string const & name, Group * const parent )
  : MultiFluidPVTPackageWrapper( name, parent )
{
  getWrapperBase( viewKeyStruct::componentMolarWeightString )->setInputFlag(InputFlags::REQUIRED);
  getWrapperBase( viewKeyStruct::phaseNamesString )->setInputFlag(InputFlags::REQUIRED);

  registerWrapper( viewKeyStruct::surfaceDensitiesString, &m_surfaceDensities, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("List of surface densities for each phase");

  registerWrapper( viewKeyStruct::tableFilesString, &m_tableFiles, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setRestartFlags(RestartFlags::NO_WRITE)->
    setDescription("List of filenames with input PVT tables");

  registerWrapper( viewKeyStruct::fluidTypeString, &m_fluidTypeString, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Type of black-oil fluid (LiveOil/DeadOil)");
}

BlackOilFluid::~BlackOilFluid()
{

}

void
BlackOilFluid::DeliverClone( string const & name,
                             Group * const parent,
                             std::unique_ptr<ConstitutiveBase> & clone ) const
{
  std::unique_ptr< BlackOilFluid > newModel = std::make_unique<BlackOilFluid>( name, parent );

  newModel->m_useMass = this->m_useMass;

  newModel->m_componentNames       = this->m_componentNames;
  newModel->m_componentMolarWeight = this->m_componentMolarWeight;

  newModel->m_phaseNames           = this->m_phaseNames;
  newModel->m_pvtPackagePhaseTypes = this->m_pvtPackagePhaseTypes;

  newModel->m_surfaceDensities = this->m_surfaceDensities;
  newModel->m_tableFiles       = this->m_tableFiles;

  newModel->m_fluidTypeString = this->m_fluidTypeString;
  newModel->m_fluidType       = this->m_fluidType;

  newModel->createFluid();

  clone = std::move( newModel );
}

void BlackOilFluid::PostProcessInput()
{
  // TODO maybe use different names?
  m_componentNames = m_phaseNames;

  MultiFluidPVTPackageWrapper::PostProcessInput();

  localIndex const NP = numFluidPhases();

#define BOFLUID_CHECK_INPUT_LENGTH( data, expected, attr ) \
  if (integer_conversion<localIndex>((data).size()) != integer_conversion<localIndex>(expected)) \
  { \
    GEOSX_ERROR( "BlackOilFluid: invalid number of entries in " \
                << (attr) << " attribute (" \
                << (data).size() << "given, " \
                << (expected) << " expected)" ); \
  }

  BOFLUID_CHECK_INPUT_LENGTH( m_surfaceDensities, NP, viewKeyStruct::surfaceDensitiesString )
  BOFLUID_CHECK_INPUT_LENGTH( m_tableFiles, NP, viewKeyStruct::surfaceDensitiesString )

#undef BOFLUID_CHECK_INPUT_LENGTH

  m_fluidType = stringToFluidType(m_fluidTypeString);
}

void BlackOilFluid::createFluid()
{
  std::vector<PHASE_TYPE> phases( m_pvtPackagePhaseTypes.begin(), m_pvtPackagePhaseTypes.end() );
  std::vector<std::string> tableFiles( m_tableFiles.begin(), m_tableFiles.end() );
  std::vector<double> densities( m_surfaceDensities.begin(), m_surfaceDensities.end() );
  std::vector<double> molarWeights( m_componentMolarWeight.begin(), m_componentMolarWeight.end() );

  switch (m_fluidType)
  {
    case FluidType::LiveOil:
      m_fluid = std::make_unique<BlackOilMultiphaseSystem>( phases, tableFiles, densities, molarWeights );
      break;
    case FluidType::DeadOil:
      m_fluid = std::make_unique<DeadOilMultiphaseSystem>( phases, tableFiles, densities, molarWeights );
      break;
    default:
      GEOSX_ERROR("Unknown fluid type");
  }
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, BlackOilFluid, std::string const &, Group * const )
} // namespace constitutive

} // namespace geosx
