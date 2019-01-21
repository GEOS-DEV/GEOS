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
  * @file BlackOilFluid.cpp
  */

#include "BlackOilFluid.hpp"

#include "codingUtilities/Utilities.hpp"
#include "managers/ProblemManager.hpp"
#include "fileIO/utils/utils.hpp"

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

namespace
{

BlackOilFluid::FluidType stringToFluidType(string const & str)
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
    GEOS_ERROR("Unrecognized black-oil fluid type: " << str);
  }
  return BlackOilFluid::FluidType::LiveOil; // keep compilers happy
}

}

BlackOilFluid::BlackOilFluid( std::string const & name, ManagedGroup * const parent )
  : MultiFluidPVTPackageWrapper( name, parent )
{
  getWrapperBase( viewKeyStruct::componentMolarWeightString )->setInputFlag(InputFlags::REQUIRED);
  getWrapperBase( viewKeyStruct::phaseNamesString )->setInputFlag(InputFlags::REQUIRED);

  RegisterViewWrapper( viewKeyStruct::surfaceDensitiesString, &m_surfaceDensities, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("List of surface densities for each phase");

  RegisterViewWrapper( viewKeyStruct::tableFilesString, &m_tableFiles, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("List of filenames with input PVT tables");

  RegisterViewWrapper( viewKeyStruct::fluidTypeString, &m_fluidTypeString, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Type of black-oil fluid (LiveOil/DeadOil)");
}

BlackOilFluid::~BlackOilFluid()
{

}

std::unique_ptr<ConstitutiveBase>
BlackOilFluid::DeliverClone( string const & name, ManagedGroup * const parent ) const
{
  std::unique_ptr< BlackOilFluid > clone = std::make_unique<BlackOilFluid>( name, parent );

  clone->m_useMass = this->m_useMass;

  clone->m_componentNames       = this->m_componentNames;
  clone->m_componentMolarWeight = this->m_componentMolarWeight;

  clone->m_phaseNames           = this->m_phaseNames;
  clone->m_pvtPackagePhaseTypes = this->m_pvtPackagePhaseTypes;

  clone->m_surfaceDensities = this->m_surfaceDensities;
  clone->m_tableFiles       = this->m_tableFiles;

  clone->m_fluidTypeString = this->m_fluidTypeString;
  clone->m_fluidType       = this->m_fluidType;

  clone->createFluid();

  return std::move( clone );
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
    GEOS_ERROR( "BlackOilFluid: invalid number of entries in " \
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

  // if table file names are not absolute paths, convert them to such, based on path to main input/restart file
  ProblemManager const * const problemManager = this->GetGroupByPath<ProblemManager>("/");
  if (problemManager != nullptr)
  {
    // hopefully at least one of input or restart file names is provided, otherwise '.' will be used
    string inputFileName = problemManager->getInputFileName();
    if (inputFileName.empty())
    {
      inputFileName = problemManager->getRestartFileName();
    }
    string inputFileDir;
    splitPath( inputFileName, inputFileDir, inputFileName );

    // if table file names are not full paths, convert them to such
    for (std::string & filename : tableFiles)
    {
      if (!isAbsolutePath(filename))
      {
        getAbsolutePath( inputFileDir + '/' + filename, filename );
      }
    }
  }

  switch (m_fluidType)
  {
    case FluidType::LiveOil:
      m_fluid = new BlackOilMultiphaseSystem(phases, tableFiles, densities, molarWeights);
      break;
    case FluidType::DeadOil:
      m_fluid = new DeadOilMultiphaseSystem(phases, tableFiles, densities, molarWeights);
      break;
  }
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, BlackOilFluid, std::string const &, ManagedGroup * const )
} // namespace constitutive

} // namespace geosx
