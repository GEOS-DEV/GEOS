/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TasksBase.hpp
 */

#ifndef TASKBASE_HPP_
#define TASKBASE_HPP_

#include <string>
#include <limits>

#include "dataRepository/ExecutableGroup.hpp"
#include "common/DataTypes.hpp"
namespace geosx
{

/**
 * @class TaskBase
 * @brief A class to execute a task when directed by the TasksManager
 */
class TaskBase : public ExecutableGroup
{
public:
  /// @copydoc geosx::dataRepository::Group::Group(std::string const & name, Group * const parent)
  explicit TaskBase( std::string const & name,
                     Group * const parent );
  /**
   * @brief Universal reference copy constructor
   * @param The other TaskBase
   */
  TaskBase( TaskBase && ) = default;

  /// Destructor
  virtual ~TaskBase() override;

  /// Delete default constructor
  TaskBase() = delete;
  /// Delete copy constructor
  TaskBase( TaskBase const & ) = delete;
  /// Delete assignment constructor
  TaskBase & operator=( TaskBase const & ) = delete;
  /// Delete universal reference assignment constructor
  TaskBase & operator=( TaskBase && ) = delete;

  /**
   * @brief Catalog name interface
   * @return This type's catalog name
   */
  static string CatalogName() { return "TaskBase"; }

  using CatalogInterface = dataRepository::CatalogInterface< TaskBase, std::string const &, Group * const >;
  /// Get the catalog interface for the TaskBase
  static CatalogInterface::CatalogType & GetCatalog();

  /// @copydoc geosx::dataRepository::Group::PostProcessInput( )
  void PostProcessInput() override;
};

} /* namespace */

#endif /* TASKBASE_HPP_ */
