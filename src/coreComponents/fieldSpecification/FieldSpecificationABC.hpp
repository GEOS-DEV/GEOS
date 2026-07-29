/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FieldSpecificationABC.hpp
 */

#ifndef GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONABC_HPP
#define GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONABC_HPP


#include "common/DataTypes.hpp"
#include "dataRepository/Group.hpp"
#include "functions/FunctionBase.hpp"
#include "functions/FunctionManager.hpp"

namespace geos
{
class Function;


/**
 * @class FieldSpecificationABC
 *
 * Abstract Base Class for concrete field modifiers (`FieldSpecification`) and high-level user-defined field specification.
 */
class FieldSpecificationABC : public dataRepository::Group
{
public:

  /**
   * @defgroup alias and functions to define statically initialized catalog
   * @{
   */

  /**
   * alias to define the catalog type for this abstract type
   */
  using CatalogInterface = dataRepository::CatalogInterface< FieldSpecificationABC,
                                                             string const &,
                                                             dataRepository::Group * const >;

  /**
   * @brief static function to return static catalog.
   * @return the static catalog to create derived types through the static factory methods.
   */
  static CatalogInterface::CatalogType & getCatalog();

  /**
   * @brief return the catalog name
   * @return the catalog name
   */
  virtual const string getCatalogName() const = 0;

  /**
   * @}
   */


  /**
   * @brief constructor
   * @param name the name of the FieldSpecificationABC in the data repository
   * @param parent the parent group of this group.
   */
  FieldSpecificationABC( string const & name, dataRepository::Group * parent );

  /**
   * destructor
   */
  virtual ~FieldSpecificationABC() override;


  /// Deleted copy constructor
  FieldSpecificationABC( FieldSpecificationABC const & ) = delete;

  /// Defaulted move constructor
  FieldSpecificationABC( FieldSpecificationABC && ) = default;

  /// deleted copy assignment
  FieldSpecificationABC & operator=( FieldSpecificationABC const & ) = delete;

  /// deleted move assignement
  FieldSpecificationABC & operator=( FieldSpecificationABC && ) = delete;

  /**
   * @brief View keys
   */
  struct viewKeyStruct
  {};

};

}

#endif //GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONABC_HPP
