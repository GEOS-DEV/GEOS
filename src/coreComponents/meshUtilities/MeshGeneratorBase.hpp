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
 * @file MeshGeneratorBase.hpp
 */

#ifndef GEOSX_MESHUTILITIES_MESHGENERATORBASE_HPP
#define GEOSX_MESHUTILITIES_MESHGENERATORBASE_HPP

#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"

namespace geosx
{

namespace dataRepository
{}

class NodeManager;
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
   * @brief Destructor for MeshGenerator
   */
  virtual ~MeshGeneratorBase();

  /**
   * @brief Return the name of the MeshGenerator in object catalog.
   * @return string that contains the catalog name of the MeshGenerator
   */
  static string catalogName() { return "MeshGeneratorBase"; }

  /// using alias for templated Catalog meshGenerator type
  using CatalogInterface = dataRepository::CatalogInterface< MeshGeneratorBase, string const &, Group * const >;

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
};
}

#endif /* GEOSX_MESHUTILITIES_MESHGENERATORBASE_HPP */
