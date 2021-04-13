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
 * @file FunctionManager.hpp
 */

#ifndef GEOSX_FUNCTIONS_FUNCTIONMANAGER_HPP_
#define GEOSX_FUNCTIONS_FUNCTIONMANAGER_HPP_

#include "FunctionBase.hpp"

#include "dataRepository/Group.hpp"

namespace geosx
{

/**
 * @class FunctionManager
 *
 * A class to manage arbitrary N-dimensional functions.
 */
class FunctionManager : public dataRepository::Group
{
public:
  /// @copydoc geosx::dataRepository::Group::Group( string const & name, Group * const parent )
  FunctionManager( const string & name,
                   dataRepository::Group * const parent );

  /**
   * @brief destructor
   */
  virtual ~FunctionManager() override;

  /**
   * @brief @return A pointer to the FunctionManager.
   */
  static FunctionManager & getInstance();

  /**
   * @brief Static Factory Catalog Functions
   * @return the catalog name
   */
  static string catalogName() { return "FunctionManager"; }

  /**
   * @brief Create a new FunctionManager object as a child of this group.
   * @param functionCatalogKey the catalog key of the new FunctionManager derived type to create
   * @param functionName the name of the new FunctionManager object in the repository
   * @return the group child
   */
  virtual Group * createChild( string const & functionCatalogKey, string const & functionName ) override;

  /**
   * @brief This function is used to expand any catalogs in the data structure
   */
  virtual void expandObjectCatalogs() override;

private:
  static FunctionManager * m_instance;
};

} /* namespace geosx */

#endif /* GEOSX_FUNCTIONS_FUNCTIONMANAGER_HPP_ */
