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
  auto declarationNode = m_document.append_child( pugi::node_declaration );
  declarationNode.append_attribute( "version" ) = "1.0";

  // Declaration of the node VTKFile
  auto vtkFileNode = m_document.append_child( "VTKFile" );
  vtkFileNode.append_attribute( "type" ) = "vtkMultiBlockDataSet";
  vtkFileNode.append_attribute( "version" ) = "1.0";

  m_blockRoot = vtkFileNode.append_child( "vtkMultiBlockDataSet" );
}

void VTKVTMWriter::write() const
{
  m_document.save_file( m_filePath.c_str() );
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
      node = node.append_child( "Block" );
      node.append_attribute( "name" ) = blockName.c_str();
    }
  }
  node = node.append_child( "DataSet" );
  node.append_attribute( "name" ) = dataSetName.c_str();
  node.append_attribute( "file" ) = filePath.c_str();
}

} // namespace vtk
} // namespace geos
