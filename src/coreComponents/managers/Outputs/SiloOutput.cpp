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
 * @file SiloOutput.cpp
 */

#include "SiloOutput.hpp"
#include "fileIO/silo/SiloFile.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/Functions/NewFunctionManager.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

SiloOutput::SiloOutput( std::string const & name,
                        ManagedGroup * const parent ):
  OutputBase( name, parent),
  m_plotFileRoot(),
  m_writeFaceMesh(),
  m_plotLevel()
{
  RegisterViewWrapper(viewKeysStruct::plotFileRoot, &m_plotFileRoot, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("");

  RegisterViewWrapper(viewKeysStruct::writeFEMFaces, &m_writeFaceMesh, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("");

  RegisterViewWrapper(viewKeysStruct::plotLevel, &m_plotLevel, false )->
    setApplyDefaultValue(1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("");

}

SiloOutput::~SiloOutput()
{}



void SiloOutput::Execute(real64 const time_n,
                         real64 const dt,
                         integer const cycleNumber,
                         integer const eventCounter,
                         real64 const eventProgress,
                         ManagedGroup * domain)
{
  DomainPartition* domainPartition = ManagedGroup::group_cast<DomainPartition*>(domain);
  SiloFile silo;

  integer rank;
  MPI_Comm_rank(MPI_COMM_GEOSX, &rank);
  MPI_Barrier( MPI_COMM_GEOSX );

  integer numFiles = this->parallelThreads();

  silo.setPlotLevel( getReference<integer>( viewKeysStruct::plotLevel ) );

  silo.Initialize(PMPIO_WRITE , numFiles );
  silo.WaitForBatonWrite(rank, cycleNumber, eventCounter, false );
  silo.WriteDomainPartition( *domainPartition, cycleNumber,  time_n + dt * eventProgress, 0);
  silo.HandOffBaton();
  silo.ClearEmptiesFromMultiObjects(cycleNumber);
  silo.Finish();

}


REGISTER_CATALOG_ENTRY( OutputBase, SiloOutput, std::string const &, ManagedGroup * const )
} /* namespace geosx */
