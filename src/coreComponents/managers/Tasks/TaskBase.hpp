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
 * @file TaskBase.hpp
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
 * A class to execute a task when directed by the TasksManager
 */
class TaskBase : public ExecutableGroup
{
public:
  /// @copydoc geosx::dataRepository::Group::Group(std::string const & name, Group * const parent)
  explicit TaskBase( std::string const & name,
                     Group * const parent );
  virtual ~TaskBase( ) override;

  /**
   * @brief Catalog name interface
   * @return This type's catalog name
   */
  static string CatalogName() { return "TaskBase"; }

  /// The catalog interface type for TaskBase
  using CatalogInterface = dataRepository::CatalogInterface< TaskBase, std::string const &, Group * const >;
  /**
   * @brief Get the catalog interface for the TaskBase
   * @return the Catalog for TaskBase.
   */
  static CatalogInterface::CatalogType & GetCatalog();

  /// @copydoc geosx::ExecutableGroup::Execute
  virtual void Execute( real64 const timeN,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        dataRepository::Group * domain ) override
  {
    GEOSX_UNUSED_VAR( timeN );
    GEOSX_UNUSED_VAR( dt );
    GEOSX_UNUSED_VAR( cycleNumber );
    GEOSX_UNUSED_VAR( eventCounter );
    GEOSX_UNUSED_VAR( eventProgress );
    GEOSX_UNUSED_VAR( domain );
    GEOSX_ERROR( "NOT IMPLEMENTED" );
  }

  /// @copydoc geosx::dataRepository::Group::PostProcessInput( )
  void PostProcessInput() override;
};

} /* namespace */

#endif /* TASKBASE_HPP_ */
