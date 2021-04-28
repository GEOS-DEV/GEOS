/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file VTKMeshGenerator.hpp
 */

#ifndef GEOSX_MESHUTILITIES_VTKMESHGENERATOR_HPP
#define GEOSX_MESHUTILITIES_VTKMESHGENERATOR_HPP

#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"
#include "codingUtilities/StringUtilities.hpp"

#include "MeshGeneratorBase.hpp"
#include "mesh/CellBlockManager.hpp"

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

namespace geosx
{


/**
 *  @class VTKMeshGenerator
 *  @brief The VTKMeshGenerator class provides a class implementation of VTK generated meshes.
 */
class VTKMeshGenerator : public MeshGeneratorBase
{
public:
/**
 * @brief Main constructor for MeshGenerator base class.
 * @param[in] name of the VTKMeshGenerator object
 * @param[in] parent the parent Group pointer for the MeshGenerator object
 */
  VTKMeshGenerator( const string & name,
                       Group * const parent );

  virtual ~VTKMeshGenerator() override;

/**
 * @brief Return the name of the VTKMeshGenerator in object Catalog.
 * @return string that contains the key name to VTKMeshGenerator in the Catalog
 */
  static string catalogName() { return "VTKMeshGenerator"; }

///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    constexpr static char const * filePathString() { return "file"; }
  };
/// @endcond

  /**
   * @brief Create a new geometric object (box, plane, etc) as a child of this group.
   * @param childKey the catalog key of the new geometric object to create
   * @param childName the name of the new geometric object in the repository
   * @return the group child
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  virtual void generateMesh( DomainPartition & domain ) override;

protected:

  /**
   * @brief This function provides capability to post process input values prior to
   * any other initialization operations.
   */
  void postProcessInput() override final;

private:

  /**
   * @brief Get the number of points of a cell knowing its vtk cell type
   * @param[in] cellType the vtk cell type
   * @return the number of points contained in the cell
   */
  localIndex GetNumberOfPoints(int cellType );

  /**
   * @brief Write a CellBlock of a given cell type
   * @param[in] name the name of the cellBlock to be written
   * @param[in] nbcells number of cells the CellBlock will contain
   * @param[in] region_id the id of the region
   * @param[in] cellType the vtk cell type for cells of the CellBlock being written
   * @param[in] cellBlockManager the CellBlockManager
   * @param[in] arraysTobeImported the list of arrays to be imported
   */
  void WriteCellBlock( string const & name, localIndex nbCells, int region_id, int cellType,
                       CellBlockManager & cellBlockManager,
                       std::vector< vtkDataArray * > const & arraysTobeImported ); 

  /**
   * @brief Write the hexahedron vertices
   * @details The node ordering from VTK differs from the node ordering in GEOSX
   * @param[in,out] cellToVertex list of nodes organized per cells
   */
  void WriteHexahedronVertices( CellBlock::NodeMapType & cellToVertex, int region_id );

  /**
   * @brief Write the pyramid vertices
   * @details The node ordering from VTK differs from the node ordering in GEOSX
   * @param[in,out] cellToVertex list of nodes organized per cells
   */
  void WritePyramidVertices( CellBlock::NodeMapType & cellToVertex );

private:

  /// Smart Pointer to the Mesh in the data structure of VTK.
  vtkSmartPointer< vtkUnstructuredGrid >  m_vtkMesh;

  /// Path to the mesh file
  Path m_filePath;
};

}

#endif /* GEOSX_MESHUTILITIES_VTKMESHGENERATOR_HPP */
