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

/**
 * @file MeshGeneratorBase.hpp
 */

#ifndef GEOSX_MESH_GENERATORS_MESHGENERATORBASE_HPP
#define GEOSX_MESH_GENERATORS_MESHGENERATORBASE_HPP

#include <string>

#include "dataRepository/Group.hpp"
#include "dataRepository/WrapperBase.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"


namespace geosx
{

namespace dataRepository
{}

class DomainPartition;

/**
 *  @class MeshGeneratorBase
 *  @brief The MeshGeneratorBase class provides an abstract base class implementation for different mesh types.
 *	   The MeshGeneratorBase is the Group specialization for different type of mesh handling.
 */
class MeshGeneratorBase : public dataRepository::Group
{
public:

  /**
   * @brief Main constructor for MeshGenerator base class.
   * @param[in] name of the MeshGenerator object
   * @param[in] parent the parent Group pointer for the MeshGenerator object
   */
  explicit MeshGeneratorBase( string const & name,
                              Group * const parent );

  /**
   * @brief Return the name of the MeshGenerator in object catalog.
   * @return string that contains the catalog name of the MeshGenerator
   */
  static string catalogName() { return "MeshGeneratorBase"; }

  /// using alias for templated Catalog meshGenerator type
  using CatalogInterface = dataRepository::CatalogInterface< MeshGeneratorBase, string const &, Group * const >;

  /**
   * @brief Create a new geometric object (box, plane, etc) as a child of this group.
   * @param childKey the catalog key of the new geometric object to create
   * @param childName the name of the new geometric object in the repository
   * @return the group child
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /**
   * @brief Accessor for the singleton Catalog object
   * @return a static reference to the Catalog object
   */
  static CatalogInterface::CatalogType & getCatalog();

  /**
   * @brief Generate the mesh object the input mesh object.
   * @param[in] domain the domain partition from which to construct the mesh object
   */
  virtual void generateMesh( DomainPartition & domain ) = 0;

  /**
   * @brief import fields from the mesh  on the array accessible via the given wrapper.
   * @param cellBlockName name of the cell block to copy data from.
   * @param meshFieldName name of the field in the meshd
   * @param wrapper Wrapper to access the array
   * @param importMaterial Indicate if we want to import material or regular fields
   */
  virtual void importFieldsOnArray( string const cellBlockName, string const meshFieldName, dataRepository::WrapperBase & wrapper, bool importMaterial ) const = 0;

  /**
   * @brief Free internal resources associated with mesh/data import.
   *
   * This is relevant for mesh generators that load external mesh files.
   * Once this method is called, they can release any memory allocated.
   */
  virtual void freeResources() {}

  virtual string getSourceName( localIndex index ) const { return ""; }//= 0;


 /**
  * @brief Get the name mapping between mesh field names and Internal GEOSX field names.
  * @return The string to string mapping of field names.
  */
std::map< string, string > getFieldsMapping( ) { return m_fieldsMapping; }

protected:
  /// Mesh to GEOSX field names mapping
  std::map< string, string > m_fieldsMapping;
};
}

#endif /* GEOSX_MESH_GENERATORS_MESHGENERATORBASE_HPP */
