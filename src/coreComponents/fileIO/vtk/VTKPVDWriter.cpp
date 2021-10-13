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

#include "VTKPVDWriter.hpp"

#include "common/MpiWrapper.hpp"

namespace geosx
{
namespace vtk
{
VTKPVDWriter::VTKPVDWriter( string fileName ):
  m_fileName( std::move( fileName ) )
{
  // Declaration of XML version
  auto declarationNode = m_pvdFile.append_child( pugi::node_declaration );
  declarationNode.append_attribute( "version" ) = "1.0";

  // Declaration of the node VTKFile
  auto vtkFileNode = m_pvdFile.append_child( "VTKFile" );
  vtkFileNode.append_attribute( "type" ) = "Collection";
  vtkFileNode.append_attribute( "version" ) = "0.1";

  vtkFileNode.append_child( "Collection" );
}

void VTKPVDWriter::setFileName( string fileName )
{
  m_fileName = std::move( fileName );
}

void VTKPVDWriter::save() const
{
  int const mpiRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  if( mpiRank == 0 )
  {
    m_pvdFile.save_file( m_fileName.c_str() );
  }
}

void VTKPVDWriter::addData( real64 time, string const & filePath ) const
{
  auto collectionNode = m_pvdFile.child( "VTKFile" ).child( "Collection" );
  auto dataSetNode = collectionNode.append_child( "DataSet" );
  dataSetNode.append_attribute( "timestep" ) = time;
  dataSetNode.append_attribute( "file" ) = filePath.c_str();
}
}
}
