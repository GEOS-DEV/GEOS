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


#include "VTKVTMWriter.hpp"

#include "common/MpiWrapper.hpp"

namespace geos
{
namespace vtk
{

VTKVTMWriter::VTKVTMWriter( string filePath )
  : m_filePath( std::move( filePath ) )
{
  // Declaration of XML version
  auto declarationNode = m_document.appendChild( pugi::node_declaration );
  declarationNode.appendAttribute( "version" ) = "1.0";

  // Declaration of the node VTKFile
  auto vtkFileNode = m_document.appendChild( "VTKFile" );
  vtkFileNode.appendAttribute( "type" ) = "vtkMultiBlockDataSet";
  vtkFileNode.appendAttribute( "version" ) = "1.0";

  m_blockRoot = vtkFileNode.appendChild( "vtkMultiBlockDataSet" );
}

void VTKVTMWriter::write() const
{
  m_document.saveFile( m_filePath );
}

void VTKVTMWriter::addDataSet( std::vector< string > const & blockPath,
                               string const & dataSetName,
                               string const & filePath ) const
{
  auto node = m_blockRoot;
  for( string const & blockName : blockPath )
  {
    auto const n = node.find_child_by_attribute( "Block", "name", blockName.c_str() );
    if( n )
    {
      node = n;
    }
    else
    {
      node = node.appendChild( "Block" );
      node.appendAttribute( "name" ) = blockName.c_str();
    }
  }
  node = node.appendChild( "DataSet" );
  node.appendAttribute( "name" ) = dataSetName.c_str();
  node.appendAttribute( "file" ) = filePath.c_str();
}

} // namespace vtk
} // namespace geos
