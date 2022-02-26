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

#ifndef GEOSX_MESH_GENERATORS_VTKMESHGENERATOR_HPP
#define GEOSX_MESH_GENERATORS_VTKMESHGENERATOR_HPP

#include "codingUtilities/StringUtilities.hpp"
#include "codingUtilities/Utilities.hpp"
#include "mesh/ElementType.hpp"
#include "mesh/generators/ExternalMeshGeneratorBase.hpp"

// TODO can we remove this and use unique_ptr to hold mesh?
#include <vtkSmartPointer.h>

#include <map>

class vtkUnstructuredGrid;
class vtkDataArray;

namespace geosx
{

class CellBlockManager;
class ElementRegionManager;

/**
 *  @class VTKMeshGenerator
 *  @brief The VTKMeshGenerator class provides a class implementation of VTK generated meshes.
 */
class VTKMeshGenerator : public ExternalMeshGeneratorBase
{
public:

  /**
   * @brief Main constructor for MeshGenerator base class.
   * @param[in] name of the VTKMeshGenerator object
   * @param[in] parent the parent Group pointer for the MeshGenerator object
   */
  VTKMeshGenerator( const string & name,
                    Group * const parent );

/**
 * @brief Return the name of the VTKMeshGenerator in object Catalog.
 * @return string that contains the key name to VTKMeshGenerator in the Catalog
 */
  static string catalogName() { return "VTKMesh"; }

  /**
   * @brief Generate the mesh using the VTK library.
   * @param[in] domain the DomainPartition to be written
   * @details This method leverages the VTK library to load the meshes.
   * The supported formats are the official VTK ones dedicated to
   * unstructured grids (.vtu, .pvtu and .vtk).
   *
   * Please note that this mesh generator works only with a number of MPI processes than
   * can be decomposed into a power of 2.
   *
   * - If a .vtu of .vtk file is used, the root MPI process will load it.
   *   The mesh will be then redistribute among all the available MPI processes
   * - If a .pvtu file is used, it means that the mesh is pre-partionned in the file system.
   *   The available MPI processes will load the pre-partionned mesh. The mesh will be then
   *   redistributed among ALL the available MPI processes.
   *
   * The properties on the mesh will be also and redistributed. The only compatible types are double and float.
   * The properties can be multi-dimensional.
   * The name of the properties has to have the right name in order to be used by GEOSX. For instance,
   * the property that stored the input porosity in GEOSX is named "referencePorosity", so the mesh has to have
   * a property names "referencePorosity".
   *
   * The regions are defined using a property called "attribute" that can be defined in the input mesh. This property
   * will be held by each volume elements. This method will created several CellBlocks, named using the combination
   * of the attribute index and the type of the element.
   * For instance, the cells of a mesh with two regions will hold the attribute "1", or "2". The CellBlocks will
   * be instantiated according to the attribute and the type of the cells. If the region "1" has wedges, tetrahedron
   * and hexahedron, three CellBlocks will be created names 1_tetrahedron, 1_wedges and 1_hexahedron.
   * The ElementRegions have to be be defined in the XML file.
   *
   * The node sets of surface are defined in the same way, using the same
   * property names "attribute" defined in the input mesh. The node sets will
   * hold a name that is just the attribute index. For instance, if a mesh has three
   * surfaces of interest, with triangles and/or quads holding an attribute value
   * of 1, 2 or 3, three node sets named "1", "2" and "3" will be instantiated by this method
   */
  virtual void generateMesh( DomainPartition & domain ) override;

  virtual void importFields( DomainPartition & domain ) const override;

  virtual void freeResources() override;

  /// Type of map used to store cell lists
  using CellMapType = std::map< ElementType, std::unordered_map< int, std::vector< vtkIdType > > >;

private:

  void buildCellBlocks( CellBlockManager & cellBlockManager ) const;

  void buildSurfaces( CellBlockManager & cellBlockManager ) const;

  void importFieldOnCellElementSubRegion( int const regionId,
                                          ElementType const elemType,
                                          std::vector< vtkIdType > const & cellIds,
                                          ElementRegionManager & elemManager,
                                          arrayView1d< string const > const & fieldNames,
                                          std::vector< vtkDataArray * > const & srcArrays ) const;

  ///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    constexpr static char const * regionAttributeString() { return "regionAttribute"; }
  };
  /// @endcond

  /**
   * @brief The VTK mesh to be imported into GEOSX.
   * @note We keep this smart pointer as a member for use in @p importFields().
   */
  vtkSmartPointer< vtkUnstructuredGrid > m_vtkMesh;

  /// Name of VTK dataset attribute used to mark regions
  string m_attributeName;

  /// Lists of VTK cell ids, organized by element type, then by region
  CellMapType m_cellMap;
};

} // namespace geosx

#endif /* GEOSX_MESH_GENERATORS_VTKMESHGENERATOR_HPP */
