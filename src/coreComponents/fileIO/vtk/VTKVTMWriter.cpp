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


#include "VTKVTMWriter.hpp"

#include "mpiCommunications/MpiWrapper.hpp"

namespace geosx
{
namespace vtk
{
VTKVTMWriter::VTKVTMWriter( string const & filePath ):
  m_filePath( filePath )
{
  // Declaration of XML version
  auto declarationNode = m_vtmFile.append_child( pugi::node_declaration );
  declarationNode.append_attribute( "version" ) = "1.0";

  // Declaration of the node VTKFile
  auto vtkFileNode = m_vtmFile.append_child( "VTKFile" );
  vtkFileNode.append_attribute( "type" ) = "vtkMultiBlockDataSet";
  vtkFileNode.append_attribute( "version" ) = "1.0";

  vtkFileNode.append_child( "vtkMultiBlockDataSet" );
}

/*!
 * @brief Triggers the file output
 */
void VTKVTMWriter::Save() const
{
  int const mpiRank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  if( mpiRank == 0 )
  {
    m_vtmFile.save_file( m_filePath.c_str() );
  }
}

/*!
 * @brief Add a block to the VTM file
 * @details The first level of block is for the ElementRegion (\p blockName can
 * be CellElementRegion, FaceElementRegion or WellElementREgion)
 * @param[in] blockName Name of the block
 */
void VTKVTMWriter::AddBlock( string const & blockName ) const
{
  auto vtkMultiBlockNode = m_vtmFile.child( "VTKFile" ).child( "vtkMultiBlockDataSet" );
  auto blockNode = vtkMultiBlockNode.append_child( "Block" );
  blockNode.append_attribute( "name" ) = blockName.c_str();
}

/*!
 * @brief Add a subblock to the VTM file
 * @details The second level of block is for the different Regions
 * @param[in] blockName Name of the parent block
 * @param[in] subBlockName Name of the subBlock (usually the name of the Region)
 */
void VTKVTMWriter::AddSubBlock( string const & blockName, string const & subBlockName ) const
{
  auto blockNode = m_vtmFile.child( "VTKFile" ).child( "vtkMultiBlockDataSet" ).find_child_by_attribute( "Block", "name", blockName.c_str() );
  auto subBlockNode = blockNode.append_child( "Block" );
  subBlockNode.append_attribute( "name" ) = subBlockName.c_str();
}

/*!
 * @brief Add data to the subblock \p subBlockName
 * @details The final level : paths to the vtu file per rank
 * @param[in] blockName Name of the parent block
 * @param[in] subBlockName Name of the subBlock (usually the name of the Region)
 * @param[in] filePath path to the vtu file containing the unstructured mesh
 * @param[in] mpiRank the mpi rank.
 */
void VTKVTMWriter::AddDataToSubBlock( string const & blockName, string const & subBlockName, string const & filePath, int mpiRank ) const
{
  auto blockNode = m_vtmFile.child( "VTKFile" ).child( "vtkMultiBlockDataSet" ).find_child_by_attribute( "Block", "name", blockName.c_str() );
  auto subBlockNode = blockNode.find_child_by_attribute( "Block", "name", subBlockName.c_str() );
  auto dataNode = subBlockNode.append_child( "DataSet" );
  string name = "rank_" + std::to_string( mpiRank );
  dataNode.append_attribute( "name" ) = name.c_str();
  dataNode.append_attribute( "file" ) = filePath.c_str();
}
}
}

