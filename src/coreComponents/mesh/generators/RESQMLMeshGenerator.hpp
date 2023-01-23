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
 * @file RESQMLMeshGenerator.hpp
 */

#ifndef GEOSX_MESH_GENERATORS_RESQMLMESHGENERATOR_HPP
#define GEOSX_MESH_GENERATORS_RESQMLMESHGENERATOR_HPP

#include "codingUtilities/StringUtilities.hpp"
#include "codingUtilities/Utilities.hpp"
#include "mesh/ElementType.hpp"
#include "mesh/generators/ExternalMeshGeneratorBase.hpp"
#include "mesh/generators/VTKUtilities.hpp"
#include "mesh/FieldIdentifiers.hpp"

#include "fesapi/common/EpcDocument.h"

// TODO can we remove this and use unique_ptr to hold mesh?
#include <vtkSmartPointer.h>

#include <map>
#include <unordered_map>

class vtkDataSet;

namespace geosx
{

class CellBlockManager;
class ElementRegionManager;

/**
 *  @class RESQMLMeshGenerator
 *  @brief The RESQMLMeshGenerator class provides a class implementation of RESQML generated meshes from the Fesapi library.
 */
class RESQMLMeshGenerator : public ExternalMeshGeneratorBase
{
public:

  /**
   * @brief Main constructor for MeshGenerator base class.
   * @param[in] name of the RESQMLMeshGenerator object
   * @param[in] parent the parent Group pointer for the MeshGenerator object
   */
  RESQMLMeshGenerator( const string & name,
                       Group * const parent );

  /**
   * @brief Return the name of the RESQMLMeshGenerator in object Catalog.
   * @return string that contains the key name to RESQMLMeshGenerator in the Catalog
   */
  static string catalogName() { return "RESQMLMesh"; }

  /**
   * @brief Deserializes the RESQML epc file into a Data Repository.
   * This repository is then used to populate the vtk data structures.
   */
  virtual void postProcessInput() override;

  /**
   * @brief Generate the mesh using fesapi library for reading RESQML data
   * @param[in] domain in the DomainPartition to be written
   * @details This method leverages the fesapi library to load meshes into
   * vtk data structures.
   */
  virtual void generateMesh( DomainPartition & domain ) override;

  /**
   * @brief Imports field data from RESQML data
   * @param[in] domain in the DomainPartition to be written
   */
  virtual void importFields( DomainPartition & domain ) const override;

  /**
   * @brief Free the memory of the temporary objects used to load the file.
   */
  virtual void freeResources() override;

  /**
   * @brief Get the Parent Representation object
   * @return a tuple (uuid, name) of the loaded representation
   */
  std::tuple< string, string > getParentRepresentation() const;

  /**
   * @brief Type of map used to store cell lists.
   *
   * This should be an unordered_map, but some outdated standard libraries on some systems
   * do not provide std::hash specialization for enums. This is not performance critical though.
   */
  using CellMapType = std::map< ElementType, std::unordered_map< int, std::vector< vtkIdType > > >;

private:

  /**
   * @brief Load a list of surfaces from fesapi into CellData of
   * a vtkDataSet
   * @param[in] mesh The dataset in which load the surfaces
   * @return the dataset with the surfaces
   */
  vtkSmartPointer< vtkDataSet > loadSurfaces( vtkSmartPointer< vtkDataSet > mesh );

  /**
   * @brief Load a list of regions from fesapi into CellData of
   * a vtkDataSet
   * @param[in] mesh The dataset in which load the regions
   * @return the dataset with the regions
   */
  vtkSmartPointer< vtkDataSet > loadRegions( vtkSmartPointer< vtkDataSet > mesh );

  /**
   * @brief Load a list of fields from fesapi into CellData of
   * a vtkDataSet
   * @param[in] mesh The dataset in which load the fields
   * @return the dataset with the fields
   */
  vtkSmartPointer< vtkDataSet > loadProperties( vtkSmartPointer< vtkDataSet > mesh );

  /**
   * @brief Looking for the UnstructuredGrid with fesapi and
   * Load DataObject into a vtkDataSet
   * @return the loaded object into a dataset
   */
  vtkSmartPointer< vtkDataSet > retrieveUnstructuredGrid();

  /**
   * @brief Load the RESQML data into the VTK data structure
   * @return a vtk mesh
   */
  vtkSmartPointer< vtkDataSet > loadMesh();

  ///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    constexpr static char const * uuidString() { return "UUID"; }
    constexpr static char const * titleInFileString() { return "titleInFile"; }
    constexpr static char const * uuidsRegionsToImportString() { return "UUIDsRegionsToImport"; }
    constexpr static char const * regionAttributeString() { return "regionAttribute"; }
    constexpr static char const * uuidsFieldsToImportString() { return "UUIDsFieldsToImport"; }
    constexpr static char const * uuidsSurfacesToImportString() { return "UUIDsSurfacesToImport"; }
    constexpr static char const * nodesetNamesString() { return "nodesetNames"; }
    constexpr static char const * partitionRefinementString() { return "partitionRefinement"; }
    constexpr static char const * partitionMethodString() { return "partitionMethod"; }
    constexpr static char const * useGlobalIdsString() { return "useGlobalIds"; }
  };
  /// @endcond


  ///Repository of RESQML objects
  common::DataObjectRepository * m_repository;

  /// UUID and title of the mesh
  string m_parent_uuid;
  string m_parent_title;

  /**
   * @brief The VTK mesh to be imported into GEOSX.
   * @note We keep this smart pointer as a member for use in @p importFields().
   */
  vtkSmartPointer< vtkDataSet > m_vtkMesh;

  /// UUIDs of the subrepresentation to import as regions
  string_array m_uuidsRegionsToImport;

  /// Name of VTK dataset attribute used to mark regions
  string m_attributeName;

  /// UUIDs of the fields to import
  string_array m_uuidsFieldsToImport;

  /// UUIDs of the surfaces to import
  string_array m_uuidsSurfacesToImport;

  /// Names of VTK nodesets to import
  string_array m_nodesetNames;

  /// Number of graph partitioning refinement iterations
  integer m_partitionRefinement = 0;

  /// Whether global id arrays should be used, if available
  integer m_useGlobalIds = 0;

  /// Method (library) used to partition the mesh
  vtk::PartitionMethod m_partitionMethod = vtk::PartitionMethod::parmetis;

  /// Lists of VTK cell ids, organized by element type, then by region
  CellMapType m_cellMap;
};

} // namespace geosx

#endif /* GEOSX_MESH_GENERATORS_RESQMLMESHGENERATOR_HPP */
