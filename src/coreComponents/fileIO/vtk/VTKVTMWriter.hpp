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

#ifndef GEOSX_FILEIO_VTK_VTKVTMWRITER_HPP_
#define GEOSX_FILEIO_VTK_VTKVTMWRITER_HPP_

#include "dataRepository/xmlWrapper.hpp"

namespace geosx
{
namespace vtk
{
/*!
 * @brief VTM Writer class.
 * @details a VTM file is the root file for one time step. It will contain
 * path to all the component of the mesh (Surfaces, Volumes, Wells etc.).
 */
class VTKVTMWriter
{
public:
  /*!
   * @brief Build the VTM Writer
   * @param[in] filePath path to the file
   */
  explicit VTKVTMWriter( string filePath );

  /*!
   * @brief Triggers the file output
   */
  void save() const;

  /*!
   * @brief Add a block to the VTM file
   * @details The first level of block is for the MeshBody (\p blockName will
   * be the name of the MeshBody)
   * @param[in] blockName Name of the block (will be MeshBody name)
   */
  void addBlock( string const & blockName ) const;

  /*!
   * @brief Add a block to the VTM file
   * @details The first level of block is for the ElementRegion (\p blockName can
   * be CellElementRegion, FaceElementRegion or WellElementREgion)
   * @param[in] blockName Name of the block (will be MeshBody name)
   * @param[in] subBlockName Name of the subBlock (will be one of the ElementRegion types)
   */
  void addSubBlock( string const & blockName, string const & subBlockName ) const;

  /*!
   * @brief Add a sub-sub-block to the VTM file
   * @details The third level of block is for the different user-defined Regions
   * @param[in] blockName Name of the parent subblock (will be one of the ElementRegion types)
   * @param[in] subSubBlockName Name of the subSubBlock (usually the name of the Region)
   */
  void addSubSubBlock( string const & blockName, string const & subBlockName, string const & subSubBlockName ) const;

  /*!
   * @brief Add data to the subblock \p subBlockName
   * @details The final level : paths to the vtu file per rank
   * @param[in] subBlockName Name of the parent block (will be one of the ElementRegion types)
   * @param[in] subSubBlockName Name of the subSubBlock (usually the name of the Region)
   * @param[in] filePath path to the vtu file containing the unstructured mesh
   * @param[in] mpiRank the mpi rank.
   */
  void addDataToSubSubBlock( string const & blockName, string const & subBlockName, string const & subSubBlockName, string const & filePath, int mpiRank ) const;



private:

  /// VTM XML file
  xmlWrapper::xmlDocument m_vtmFile;

  /// Path to the XML File
  string const m_filePath;
};
}
}

#endif
