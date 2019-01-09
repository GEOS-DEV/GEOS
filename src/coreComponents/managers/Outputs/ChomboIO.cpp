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
 * @file ChomboIO.cpp
 */
#include "ChomboIO.hpp"
#include "mesh/MeshLevel.hpp"
#include "managers/DomainPartition.hpp"
#include "fileIO/coupling/ChomboCoupler.hpp"
#include <string>

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

constexpr auto FILE_PATH = "geosxToChombo";

ChomboIO::ChomboIO(std::string const & name, ManagedGroup * const parent):
  OutputBase(name, parent),
  m_coupler(nullptr)
{}

ChomboIO::~ChomboIO()
{
  delete m_coupler;
  m_coupler = nullptr;
}

void ChomboIO::Execute( real64 const & time_n,
                        real64 const & dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const & eventProgress,
                        dataRepository::ManagedGroup * domain )
{
  if (m_coupler == nullptr)
  {
    DomainPartition* domainPartition = ManagedGroup::group_cast<DomainPartition*>(domain);
    MeshLevel * meshLevel = domainPartition->getMeshBody(0)->getMeshLevel(0);
    m_coupler = new ChomboCoupler(MPI_COMM_GEOSX, "", *meshLevel);
  }
  
  m_coupler->setPath(std::string(FILE_PATH) + "_" + std::to_string(cycleNumber) + ".hdf5");
  m_coupler->write(dt);
}

REGISTER_CATALOG_ENTRY(OutputBase, ChomboIO, std::string const &, ManagedGroup * const)
} /* namespace geosx */
