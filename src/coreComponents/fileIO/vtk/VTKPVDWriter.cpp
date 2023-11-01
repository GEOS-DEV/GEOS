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

namespace geos
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

// If running job from restart the pvd file should be read in to append to, not overwritten
void VTKPVDWriter::read()
{
  //CC: do I need to clear the document first?
  m_pvdFile.reset();

  // If restarting job and pvd already exists read in and append to that
  xmlWrapper::xmlResult const xmlResult = m_pvdFile.load_file( m_fileName.c_str() );
  GEOS_THROW_IF( !xmlResult, GEOS_FMT( "Errors found while parsing XML file {}\nDescription: {}\nOffset: {}",
                                       m_fileName, xmlResult.description(), xmlResult.offset ), InputError );
}

void VTKPVDWriter::save() const
{
  m_pvdFile.save_file( m_fileName.c_str() );
}

void VTKPVDWriter::addData( real64 time, string const & filePath ) const
{
  auto collectionNode = m_pvdFile.child( "VTKFile" ).child( "Collection" );
  auto dataSetNode = collectionNode.append_child( "DataSet" );
  dataSetNode.append_attribute( "timestep" ) = time;
  dataSetNode.append_attribute( "file" ) = filePath.c_str();
}

void VTKPVDWriter::reinitData() const
{
  auto collectionNode = m_pvdFile.child( "VTKFile" ).child( "Collection" );

  while( collectionNode.remove_child( "DataSet" ) )
  {}
}
}
}
