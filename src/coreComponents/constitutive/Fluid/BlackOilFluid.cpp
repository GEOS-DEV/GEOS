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
  * @file BlackOilFluid.cpp
  */

#include "BlackOilFluid.hpp"

#include "codingUtilities/Utilities.hpp"

// PVTPackage includes
#include "MultiphaseSystem/BlackOilMultiphaseSystem.hpp"


using namespace PVTPackage;

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

namespace constitutive
{

BlackOilFluid::BlackOilFluid( std::string const & name, ManagedGroup * const parent )
  : MultiFluidPVTPackageWrapper( name, parent )
{
  RegisterViewWrapper( viewKeyStruct::surfaceDensitiesString, &m_surfaceDensities, false );
  RegisterViewWrapper( viewKeyStruct::tableFilesString, &m_tableFiles, false );
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

  clone->createFluid();

  return std::move( clone );
}

void BlackOilFluid::FillDocumentationNode()
{
  MultiFluidPVTPackageWrapper::FillDocumentationNode();

  DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setShortDescription( "Black-oil fluid model based on PVTPackage" );

  docNode->AllocateChildNode( viewKeyStruct::surfaceDensitiesString,
                              viewKeyStruct::surfaceDensitiesString,
                              -1,
                              "real64_array",
                              "real64_array",
                              "List of surface densities for each phase",
                              "List of surface densities for each phase",
                              "REQUIRED",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::tableFilesString,
                              viewKeyStruct::tableFilesString,
                              -1,
                              "string_array",
                              "string_array",
                              "List of filenames with input PVT tables",
                              "List of filenames with input PVT tables",
                              "REQUIRED",
                              "",
                              1,
                              1,
                              0 );
}

void BlackOilFluid::ReadXML_PostProcess()
{
  // TODO maybe use different names?
  m_componentNames = m_phaseNames;

  MultiFluidPVTPackageWrapper::ReadXML_PostProcess();

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
}

void BlackOilFluid::createFluid()
{
  std::vector<PHASE_TYPE> phases( m_pvtPackagePhaseTypes.begin(), m_pvtPackagePhaseTypes.end() );
  std::vector<std::string> tableFiles( m_tableFiles.begin(), m_tableFiles.end() );
  std::vector<double> densities( m_surfaceDensities.begin(), m_surfaceDensities.end() );
  std::vector<double> molarWeights( m_componentMolarWeight.begin(), m_componentMolarWeight.end() );

  m_fluid = new BlackOilMultiphaseSystem( phases, tableFiles, densities, molarWeights );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, BlackOilFluid, std::string const &, ManagedGroup * const )
} // namespace constitutive

} // namespace geosx
