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
 * @file VTKCompositeMeshGenerator.hpp
 */

#ifndef GEOS_MESH_GENERATORS_VTKCOMPOSITEMESHGENERATOR_HPP
#define GEOS_MESH_GENERATORS_VTKCOMPOSITEMESHGENERATOR_HPP

#include "mesh/generators/ExternalMeshGeneratorBase.hpp"
#include "mesh/generators/VTKUtilities.hpp"

#include <vtkDataSet.h>

namespace geos
{

/**
 *  @class VTKCompositeMeshGenerator
 *  @brief The VTKCompositeMeshGenerator class provides a class implementation of VTK generated meshes.
 */
class VTKCompositeMeshGenerator : public ExternalMeshGeneratorBase
{
public:

  /**
   * @brief Main constructor for MeshGenerator base class.
   * @param[in] name of the VTKCompositeMeshGenerator object
   * @param[in] parent the parent Group pointer for the MeshGenerator object
   */
  VTKCompositeMeshGenerator( const string & name,
                             Group * const parent );

/**
 * @brief Return the name of the VTKCompositeMeshGenerator in object Catalog.
 * @return string that contains the key name to VTKCompositeMeshGenerator in the Catalog
 */
  static string catalogName() { return "VTKCompositeMesh"; }

  /**
   * @brief Generate meshes using the VTK library from composite data.
   * @param[inout] cellBlockManager the CellBlockManager that will receive the meshing information
   * @param[in] partition the number of domain in each direction (x,y,z) for InternalMesh only, not used here
   * @details This method leverages the VTK library to load the meshes from composite data sets.
   * The supported formats are the official VTK ones dedicated to
   * composite data sets (.vtm, .vtpc).
   */
  virtual void fillCellBlockManager( CellBlockManager & cellBlockManager, array1d< int >  const & partition ) override;

  void importFieldOnArray( Block block,
                           string const & blockName,
                           string const & meshFieldName,
                           bool isMaterialField,
                           dataRepository::WrapperBase & wrapper ) const override;

  virtual void freeResources() override;

private:

  ///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    constexpr static char const * regionAttributeString() { return "regionAttribute"; }
    constexpr static char const * mainBlockNameString() { return "mainBlockName"; }
    constexpr static char const * faceBlockNamesString() { return "faceBlocks"; }
    constexpr static char const * nodesetNamesString() { return "nodesetNames"; }
    constexpr static char const * useGlobalIdsString() { return "useGlobalIds"; }
  };
  /// @endcond


  /// Name of VTK dataset attribute used to mark regions
  string m_attributeName;

  /// Name of the main block to be imported (for multi-block files).
  string m_mainBlockName;

  /// Name of the face blocks to be imported (for multi-block files).
  array1d< string > m_faceBlockNames;

  /// Maps the face block name to its vtk mesh instance.
  std::map< string, vtkSmartPointer< vtkDataSet > > m_faceBlockMeshes;

  /// Names of VTK nodesets to import
  string_array m_nodesetNames;

  /// Whether global id arrays should be used, if available
  integer m_useGlobalIds = 0;

  /// Lists of VTK cell ids, organized by element type, then by region
  vtk::CellMapType m_cellMap;
};

} // namespace geos

#endif /* GEOS_MESH_GENERATORS_VTKCOMPOSITEMESHGENERATOR_HPP */
