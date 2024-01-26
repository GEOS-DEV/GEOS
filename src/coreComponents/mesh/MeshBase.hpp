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
 * @file MeshBase.hpp
 */

#ifndef GEOS_MESH_MESHBASE_HPP
#define GEOS_MESH_MESHBASE_HPP

#include "dataRepository/Group.hpp"
#include "dataRepository/WrapperBase.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"


namespace geos
{

/**
 *  @class MeshBase
 *  @brief The MeshBase class provides an abstract base class implementation for different mesh types.
 *	   The MeshBase is the Group specialization for different type of mesh handling.
 */
class MeshBase : public dataRepository::Group
{
public:

  /**
   * @brief Main constructor for MeshBase base class.
   * @param[in] name of the MeshBase object
   * @param[in] parent the parent Group pointer for the MeshBase object
   */
  explicit MeshBase( string const & name,
                     Group * const parent );

  /**
   * @brief Return the name of the MeshBase in object catalog.
   * @return string that contains the catalog name of the MeshBase
   */
  static string catalogName() { return "MeshBase"; }

  /// This function is used to expand any catalogs in the data structure
//   virtual void expandObjectCatalogs() override;

  /// using alias for templated Catalog MeshBase type
  using CatalogInterface = dataRepository::CatalogInterface< MeshBase, string const &, Group * const >;

  /**
   * @brief Create a new geometric object (box, plane, etc) as a child of this group.
   * @param childKey the catalog key of the new geometric object to create
   * @param childName the name of the new geometric object in the repository
   * @return the group child
   */
//   virtual Group * createChild( string const & childKey, string const & childName ) override;

  /**
   * @brief Accessor for the singleton Catalog object
   * @return a static reference to the Catalog object
   */
  static CatalogInterface::CatalogType & getCatalog();

};

}

#endif /* GEOS_MESH_MESHBASE_HPP */
