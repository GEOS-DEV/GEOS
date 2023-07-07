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
  auto declarationNode = m_pvdFile.appendChild( pugi::node_declaration );
  declarationNode.appendAttribute( "version" ) = "1.0";

  // Declaration of the node VTKFile
  auto vtkFileNode = m_pvdFile.appendChild( "VTKFile" );
  vtkFileNode.appendAttribute( "type" ) = "Collection";
  vtkFileNode.appendAttribute( "version" ) = "0.1";

  vtkFileNode.appendChild( "Collection" );
}

void VTKPVDWriter::setFileName( string fileName )
{
  m_fileName = std::move( fileName );
}

void VTKPVDWriter::save() const
{
  m_pvdFile.saveFile( m_fileName );
}

void VTKPVDWriter::addData( real64 time, string const & filePath ) const
{
  auto collectionNode = m_pvdFile.getChild( "VTKFile" ).getChild( "Collection" );
  auto dataSetNode = collectionNode.appendChild( "DataSet" );
  dataSetNode.appendAttribute( "timestep" ) = time;
  dataSetNode.appendAttribute( "file" ) = filePath.c_str();
}

void VTKPVDWriter::reinitData() const
{
  auto collectionNode = m_pvdFile.getChild( "VTKFile" ).getChild( "Collection" );

  while( collectionNode.remove_child( "DataSet" ) )
  {}
}
}
}
