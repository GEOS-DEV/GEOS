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
 * @file ChomboIO.cpp
 */
#include "ChomboIO.hpp"
#include "mesh/MeshLevel.hpp"
#include "managers/DomainPartition.hpp"
#include "fileIO/coupling/ChomboCoupler.hpp"
#include <string>
#include <fstream>
#include <chrono>

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

ChomboIO::ChomboIO(std::string const & name, ManagedGroup * const parent):
  OutputBase(name, parent),
  m_coupler(nullptr),
  m_outputPath(),
  m_inputPath(),
  m_waitForInput()
{
  RegisterViewWrapper(viewKeyStruct::outputPathString, &m_outputPath, false)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Path at which the geosx to chombo file will be written.");

  RegisterViewWrapper(viewKeyStruct::waitForInputString, &m_waitForInput, false)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("True iff geosx should wait for chombo to write out a file. When true the inputPath must be set.");

  RegisterViewWrapper(viewKeyStruct::inputPathString, &m_inputPath, false)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Path at which the chombo to geosx file will be written.");
}

ChomboIO::~ChomboIO()
{
  delete m_coupler;
  m_coupler = nullptr;
}

void ChomboIO::Execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        dataRepository::ManagedGroup * const domain )
{
  GEOS_LOG("m_outputPath = " << m_outputPath);
  GEOS_LOG("m_inputPath = " << m_inputPath);
  GEOS_LOG("m_waitForInput = " << m_waitForInput);

  if (m_coupler == nullptr)
  {
    DomainPartition * const domainPartition = ManagedGroup::group_cast<DomainPartition*>(domain);
    MeshLevel * const meshLevel = domainPartition->getMeshBody(0)->getMeshLevel(0);
    m_coupler = new ChomboCoupler(MPI_COMM_GEOSX, m_outputPath, m_inputPath, *meshLevel);
  }
  
  m_coupler->write(dt);

  if (m_waitForInput)
  {
    m_coupler->read();
  }
}

REGISTER_CATALOG_ENTRY(OutputBase, ChomboIO, std::string const &, ManagedGroup * const)
} /* namespace geosx */
