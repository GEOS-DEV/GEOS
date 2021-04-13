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

/*
 * @file SimpleGeometricObjectBase.hpp
 */

#ifndef GEOSX_MESH_SIMPLEGEOMETRICOBJECTS_SIMPLEGEOMETRICOBJECTBASE_HPP_
#define GEOSX_MESH_SIMPLEGEOMETRICOBJECTS_SIMPLEGEOMETRICOBJECTBASE_HPP_

#include "dataRepository/Group.hpp"
#include "codingUtilities/StringUtilities.hpp"
#include "dataRepository/ObjectCatalog.hpp"

class Function;

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const geometricObjects( "GeometricObjects" );
}
}

/**
 * @class SimpleGeometricObjectBase
 * @brief Base class for the geometric objects (box, plane, cylinder).
 */
class SimpleGeometricObjectBase : public dataRepository::Group
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor.
   * @param name name of the object in the data hierarchy.
   * @param parent pointer to the parent group in the data hierarchy.
   */
  explicit SimpleGeometricObjectBase( string const & name,
                                      Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~SimpleGeometricObjectBase();

  ///@}

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Get the catalog name.
   * @return the name of this class in the catalog
   */
  static string catalogName() { return "SimpleGeometricObjectBase"; }

  /**
   * @brief Type alias for catalog interface used by this class.
   */
  using CatalogInterface = dataRepository::CatalogInterface< SimpleGeometricObjectBase, string const &, Group * const >;

  /**
   * @copydoc catalogName()
   */
  static CatalogInterface::CatalogType & getCatalog();

  ///@}

  /**
   * @brief Check if the input coordinates are in the object.
   * @param[in] coord the coordinates to test
   * @return true if the coordinates are in the object, false otherwise
   */
  virtual bool isCoordInObject( real64 const ( &coord ) [3] ) const = 0;

};


}
#endif /* GEOSX_MESH_SIMPLEGEOMETRICOBJECTS_SIMPLEGEOMETRICOBJECTBASE_HPP_ */
