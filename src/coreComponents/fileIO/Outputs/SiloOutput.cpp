/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SiloOutput.cpp
 */

#include "SiloOutput.hpp"

#include "common/TimingMacros.hpp"
#include "fileIO/silo/SiloFile.hpp"
#include "mesh/DomainPartition.hpp"

namespace geosx
{

using namespace dataRepository;

SiloOutput::SiloOutput( string const & name,
                        Group * const parent ):
  OutputBase( name, parent ),
  m_plotFileRoot( "plot" ),
  m_writeEdgeMesh( 0 ),
  m_writeFaceMesh( 0 ),
  m_writeCellElementMesh( 1 ),
  m_writeFaceElementMesh( 1 ),
  m_plotLevel()
{
  registerWrapper( viewKeysStruct::plotFileRoot, &m_plotFileRoot ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "plot" ).
    setDescription( "" );

  registerWrapper( viewKeysStruct::writeEdgeMesh, &m_writeEdgeMesh ).
    setDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "" );

  registerWrapper( viewKeysStruct::writeFaceMesh, &m_writeFaceMesh ).
    setDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "" );

  registerWrapper( viewKeysStruct::writeCellElementMesh, &m_writeCellElementMesh ).
    setDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "" );

  registerWrapper( viewKeysStruct::writeFaceElementMesh, &m_writeFaceElementMesh ).
    setDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "" );

  registerWrapper( viewKeysStruct::plotLevel, &m_plotLevel ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "" );

}

SiloOutput::~SiloOutput()
{}



bool SiloOutput::execute( real64 const time_n,
                          real64 const dt,
                          integer const cycleNumber,
                          integer const eventCounter,
                          real64 const eventProgress,
                          DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  SiloFile silo;

  int const size = MpiWrapper::commSize( MPI_COMM_GEOSX );
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  MpiWrapper::barrier( MPI_COMM_GEOSX );

  integer const numFiles = parallelThreads() == 0 ? size : parallelThreads();

  // TODO set this during initialization
  //  silo.setOutputDirectory( getGlobalState().getCommandLineOptions().outputDirectory ),
  silo.setOutputDirectory( getOutputDirectory() ),
  silo.setPlotLevel( m_plotLevel );
  silo.setWriteEdgeMesh( m_writeEdgeMesh );
  silo.setWriteFaceMesh( m_writeFaceMesh );
  silo.setWriteCellElementMesh( m_writeCellElementMesh );
  silo.setWriteFaceElementMesh( m_writeFaceElementMesh );
  silo.setPlotFileRoot( m_plotFileRoot );
  silo.initialize( numFiles );
  silo.waitForBatonWrite( rank, cycleNumber, eventCounter, false );
  silo.writeDomainPartition( domain, cycleNumber, time_n + dt * eventProgress, 0 );
  silo.handOffBaton();
  silo.clearEmptiesFromMultiObjects( cycleNumber );
  silo.finish();

  return false;
}


REGISTER_CATALOG_ENTRY( OutputBase, SiloOutput, string const &, Group * const )
} /* namespace geosx */
