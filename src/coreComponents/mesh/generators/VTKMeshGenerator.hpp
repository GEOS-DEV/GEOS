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

#include <vtkDataArray.h>

#include <map>

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

  ~VTKMeshGenerator() override = default;

/**
 * @brief Return the name of the VTKMeshGenerator in object Catalog.
 * @return string that contains the key name to VTKMeshGenerator in the Catalog
 */
  static string catalogName() { return "VTKMeshGenerator"; }

protected:
  /**
   * @brief This function provides capability to post process input values prior to
   * any other initialization operations.
   */
  void postProcessInput() override final;

private:
///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    constexpr static char const * filePathString() { return "file"; }
  };
/// @endcond

  /**
   * @brief Generate the mesh using the VTK library.
   * @param[in] domain the DomainPartition to be written
   * @details This method leverages the VTK library to load the meshes.
   * The supported formats are the official VTK ones dedicated to
   * unstructured grids (.vtu, .pvtu and .vtk).\n\n
   *
   * Please note that this mesh generator works only with a number of MPI processes than
   * can be decomposed into a power of 2.\n\n
   *
   * - If a .vtu of .vtk file is used, the root MPI process will load it.
   *   The mesh will be then redistribute among all the avaible MPI processes
   * - If a .pvtu file is used, it means that the mesh is pre-partionned in the file system.
   *   The available MPI processes will load the pre-partionned mesh. The mesh will be then
   *   redistributed among ALL the available MPI processes.\n\n
   *
   * The properties on the mesh will be also and redistributed. The only compatible typs are double and float.
   * The properties can be multi-dimensional.\n
   * The name of the properties has to have the right name in order to be used by GEOSX. For instance,
   * the property that stored the input porosity in GEOSX is named "referencePorosity", so the mesh has to have
   * a property names "referencePorosity".\n\n
   *
   * The regions are defined using a property called "attribute" that can be defined in the input mesh. This property
   * will be held by each volume elements. This method will created several CellBlocks, named using the combination
   * of the attribute index and the type of the element.\n
   * For instance, the cells of a mesh with two regions will hold the attribute "1", or "2". The CellBlocks will
   * be instantiated according to the attribute and the type of the cells. If the region "1" has wedges, tetrahedron
   * and hexahedron, three CellBlocks will be created names 1_tetrahedron, 1_wedges and 1_hexahedron.
   * The ElementRegions have to be be defined in the XML file.\n\n
   *
   * The pointsets of surface are defined in the same way, using the same property names "attribute" defined in the
   * input mesh. The pointsets will hold a name that is just the attribute index. For instance, if a mesh has three
   * surfaces of interest, with triangles and/or quads holding an attribute value of 1, 2 or 3, three pointsets named
   * "1", "2" and "3" will be instantiated by this method
   */
  virtual void generateMesh( DomainPartition & domain ) override;

  virtual void importFields( DomainPartition & domain ) const override;

  /**
   * @brief Create a new geometric object (box, plane, etc) as a child of this group.
   * @param childKey the catalog key of the new geometric object to create
   * @param childName the name of the new geometric object in the repository
   * @return the group child
   */
  Group * createChild( string const & childKey, string const & childName ) override;

private:

  /// Path to the mesh file
  Path m_filePath;

  std::map< int, std::vector< vtkIdType > > m_regionsHex;
  std::map< int, std::vector< vtkIdType > > m_regionsTetra;
  std::map< int, std::vector< vtkIdType > > m_regionsWedges;
  std::map< int, std::vector< vtkIdType > > m_regionsPyramids;

  std::vector< vtkDataArray * > m_arraysToBeImported;
};
}
#endif /* GEOSX_MESHUTILITIES_VTKMESHGENERATOR_HPP */
