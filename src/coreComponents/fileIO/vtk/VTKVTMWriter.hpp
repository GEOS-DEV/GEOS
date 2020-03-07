/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_FILEIO_VTK_VTKVTMWRITER_HPP_
#define GEOSX_FILEIO_VTK_VTKVTMWRITER_HPP_

#include "dataRepository/xmlWrapper.hpp"

namespace geosx
{
namespace vtk
{
class VTKVTMWriter
{
  public:
  VTKVTMWriter( string const & name) :
    m_fileName( name )
  {
    // Declaration of XML version
    auto declarationNode = m_vtmFile.append_child(pugi::node_declaration);
    declarationNode.append_attribute("version") = "1.0";

    // Declaration of the node VTKFile
    auto vtkFileNode = m_vtmFile.append_child("VTKFile");
    vtkFileNode.append_attribute("type") = "vtkMultiBlockDataSet";
    vtkFileNode.append_attribute("version") = "1.0";

    vtkFileNode.append_child("vtkMultiBlockDataSet");
  }

  void Save() const
  {
    int const mpiRank = MpiWrapper::Comm_rank(MPI_COMM_GEOSX);
    if( mpiRank == 0 )
    {
       m_vtmFile.save_file( m_fileName.c_str() );
    }
  }

  void AddBlock( string const& blockName ) const
  {
    auto vtkMultiBlockNode = m_vtmFile.child("VTKFile").child("vtkMultiBlockDataSet");
    auto blockNode = vtkMultiBlockNode.append_child("Block");
    blockNode.append_attribute("name") = blockName.c_str();
  }

  void AddSubBlock( string const& blockName, string const & subBlockName ) const
  {
    auto blockNode = m_vtmFile.child("VTKFile").child("vtkMultiBlockDataSet").find_child_by_attribute("Block", "name", blockName.c_str() );
    auto subBlockNode = blockNode.append_child( "Block" );
    subBlockNode.append_attribute("name") = subBlockName.c_str();
  }

  void AddDataToSubBlock( string const& blockName, string const& subBlockName, string const& filePath, int mpiRank ) const
  {
    auto blockNode = m_vtmFile.child("VTKFile").child("vtkMultiBlockDataSet").find_child_by_attribute("Block", "name", blockName.c_str() );
    auto subBlockNode = blockNode.find_child_by_attribute("Block", "name", subBlockName.c_str() );
    auto dataNode = subBlockNode.append_child("DataSet");
    string name = "rank_" + std::to_string( mpiRank );
    dataNode.append_attribute("name") = name.c_str();
    dataNode.append_attribute("file") = filePath.c_str();
  }

  void AddDataToBlock( string const& blockName, string const& filePath ) const
  {
    int const mpiRank = MpiWrapper::Comm_rank(MPI_COMM_GEOSX);
    auto blockNode = m_vtmFile.child("VTKFile").child("vtkMultiBlockDataSet").find_child_by_attribute("Block", "name", blockName.c_str() );
    auto dataNode = blockNode.append_child("DataSet");
    string name = "rank_" + std::to_string( mpiRank );
    dataNode.append_attribute("name") = name.c_str();
    dataNode.append_attribute("file") = filePath.c_str();
  }
  private:
  
  /// VTM XML file
  xmlWrapper::xmlDocument m_vtmFile;

  /// Name of the XML File
  string const m_fileName;
};
}
};

#endif
