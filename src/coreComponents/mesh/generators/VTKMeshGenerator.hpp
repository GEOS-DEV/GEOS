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
#ifdef GEOSX_USE_MPI
#include <vtkMPIController.h>
#include <vtkMPI.h>
#else
#include <vtkDummyController.h>
#endif

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
    /*
    constexpr static char const * regionsToImportString() { return "regionsToImport"; }
    constexpr static char const * surfacesToImportString() { return "surfacesToImport"; }
    */
  };
/// @endcond

  /**
   * @brief Create a new geometric object (box, plane, etc) as a child of this group.
   * @param childKey the catalog key of the new geometric object to create
   * @param childName the name of the new geometric object in the repository
   * @return the group child
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /**
   * @brief Generate the mesh using the VTK library.
   * @details This method leverages the VTK library to load the meshes.
   * The supported formats are the official VTK ones dedicated to
   * unstructured grids (.vtu, .pvtu and .vtk).
   *
   * Please note that this mesh generator works only with a number of MPI processes than
   * can be decomposed into a power of 2.
   * 
   * - If a .vtu of .vtk file is used, the root MPI process will load it.
   *   The mesh will be then redistribute among all the avaible MPI processes
   * - If a .pvtu file is used, it means that the mesh is pre-partionned in the file system.
   *   The available MPI processes will load the pre-partionned mesh. The mesh will be then
   *   redistributed among ALL the available MPI processes.
   *
   * The properties on the mesh will be also and redistributed. The only compatible typs are double and float. 
   * The properties can be multi-dimensional.$
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
   * The pointsets of surface are defined in the same way, using the same property names "attribute" defined in the 
   * input mesh. The pointsets will hold a name that is just the attribute index. For instance, if a mesh has three
   * surfaces of interest, with triangles and/or quads holding an attribute value of 1, 2 or 3, three pointsets named
   * "1", "2" and "3" will be instantiated by this method
   * @param[in] domain the DomainPartition to be written
   */
  virtual void generateMesh( DomainPartition & domain ) override;

protected:

  /**
   * @brief This function provides capability to post process input values prior to
   * any other initialization operations.
   */
  void postProcessInput() override final;

private:

  /**
   * @brief Return a VTK controller for multiprocessing.
   */
  vtkSmartPointer<vtkMultiProcessController> GetVTKController();

  /**
   * @brief Load the VTK file into the VTK data structure
   */
  vtkSmartPointer<vtkUnstructuredGrid> LoadVTKMesh();

  /**
   * @brief Redistribute the mesh among the available MPI ranks
   * @details this method will also generate global ids for points and cells in the VTK Mesh
   * @param[in] loadedMesh the mesh that was loaded on one or several MPI ranks
   */
  void RedistributeMesh(vtkUnstructuredGrid * loadedMesh);

  /**
   * @brief Copy the VTK mesh nodes into the nodeManager of GEOSX
   * @param[in] nodeManager the NodeManager of the domain in which the poiints will be copied.
   * @return the global length of the mesh
   */
  double WriteMeshNodes(NodeManager & nodeManager);

  /**
   * @brief Compute the potential rank neighbor list
   * @details Fills the metisNeighbor list in \p domain.
   * @param[in] domain the DomainPartition in which the neighbhor list will be computed
   */
  void ComputePotentialNeighborLists( DomainPartition & domain);

  /** 
   * @brief Get the attribute data array from the VTK mesh
   * @return a pointer to the vtkIntArray containing the attributes if it exists, nullptr otherwise.
   */
  vtkIntArray * GetAttributeDataArray();

  /**
   * @brief This methos is used to preprocess the the VTK mesh and count the number of cells, facets, regions
   * and surfaces.
   * @param[out] nbHex number of hexahedra
   * @param[out] nbTet number of tetra
   * @param[out] nbWedge number of wedges
   * @param[out] nbPyr number of pyramids
   * @param[out] regions_hex map from region index to the number of hexahedron in this region
   * @param[out] regions_tetra map from region index to the number of tetra in this region
   * @param[out] regions_wedges map from region index to the number of wedges in this region
   * @param[out] regions_pyramids map from region index to the number of pyramids in this region
   * @param[out] regions a set containing all the region indexes
   * @param[out] surfaces a set containing all the surface indexes, from this MPI rank
   * @param[out] allSurfaces a vector containing all the surfaces among all the MPI rank
   */
  void CountCellsAndFaces( localIndex & nbHex, localIndex & nbTet, localIndex & nbWedge, localIndex & nbPyr,
                           std::map<int,localIndex> & regions_hex, std::map<int,localIndex> & regions_tetra,
                           std::map<int,localIndex> & regions_wedges, std::map<int,localIndex> & regions_pyramids,
                           std::set< int > & regions, std::set< int > & surfaces, std::vector<int> & allSurfaces);

  /**
   * @brief Find the properties to be imported
   * @details all the float and double vtkArrays will be imported
   * @return a vector containing all the vtkDataArray that can be imported
   */
  std::vector< vtkDataArray * > FindArrayToBeImported();

  /**
   * @brief Write all the cell blocks
   * @param[in] domain the domain in which the cell blocks will be written
   * @param[in] nbHex number of hexahedra
   * @param[in] nbTet number of tetra
   * @param[in] nbWedge number of wedges
   * @param[in] nbPyr number of pyramids
   * @param[in] regions_hex map from region index to the number of hexahedron in this region
   * @param[in] regions_tetra map from region index to the number of tetra in this region
   * @param[in] regions_wedges map from region index to the number of wedges in this region
   * @param[in] regions_pyramids map from region index to the number of pyramids in this region
   * @param[in] regions a set containing all the region indexes
   * @param[in] arraysToBeImported a vector containing all the vtkDataArray that can be imported
   */
  void WriteCellBlocks( DomainPartition & domain, localIndex nbHex, localIndex nbTet, localIndex nbWedge, localIndex nbPyr,
                        std::map<int,localIndex> & regions_hex, std::map<int,localIndex> & regions_tetra,
                        std::map<int,localIndex> & regions_wedges, std::map<int,localIndex> & regions_pyramids,
                        std::set< int > & regions, std::vector< vtkDataArray * > & arraysToBeImported);


  /**
   * @param[in] nodeManager the NodeManager of the domain in which the poiints will be copied.
   * @param[in] allSurfaces the surfaces id to be imported
   */
  void WriteSurfaces( NodeManager & nodeManager, std::vector<int> const & allSurfaces );

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
  void WriteHexahedronVertices( CellBlock::NodeMapType & cellToVertex, int region_id,  arrayView1d< globalIndex > const & localToGlobal  );

  /**
   * @brief Write the wedge vertices
   * @details The node ordering from VTK differs from the node ordering in GEOSX
   * @param[in,out] cellToVertex list of nodes organized per cells
   */
  void WriteWedgeVertices( CellBlock::NodeMapType & cellToVertex, int region_id,  arrayView1d< globalIndex > const & localToGlobal );

private:

  /// Smart Pointer to the Mesh in the data structure of VTK.
  vtkSmartPointer< vtkUnstructuredGrid >  m_vtkMesh;

  /// Path to the mesh file
  Path m_filePath;

  /*
  /// Names of the surfaces to be imported
  string_array m_surfacesToImport;

  /// Names of the regions to be imported
  string_array m_regionsToImport;
  */
};

}

#endif /* GEOSX_MESHUTILITIES_VTKMESHGENERATOR_HPP */
