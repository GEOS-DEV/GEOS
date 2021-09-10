/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TaskBase.hpp
 */

#ifndef GEOSX_EVENTS_TASKS_TASKBASE_HPP_
#define GEOSX_EVENTS_TASKS_TASKBASE_HPP_


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
  /// @copydoc geosx::dataRepository::Group::Group(string const & name, Group * const parent)
  explicit TaskBase( string const & name,
                     Group * const parent );
  virtual ~TaskBase( ) override;

  /**
   * @brief Catalog name interface
   * @return This type's catalog name
   */
  static string catalogName() { return "TaskBase"; }

  /// The catalog interface type for TaskBase
  using CatalogInterface = dataRepository::CatalogInterface< TaskBase, string const &, Group * const >;
  /**
   * @brief Get the catalog interface for the TaskBase
   * @return the Catalog for TaskBase.
   */
  static CatalogInterface::CatalogType & getCatalog();

  /// @copydoc geosx::ExecutableGroup::execute
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override
  {
    GEOSX_UNUSED_VAR( time_n );
    GEOSX_UNUSED_VAR( dt );
    GEOSX_UNUSED_VAR( cycleNumber );
    GEOSX_UNUSED_VAR( eventCounter );
    GEOSX_UNUSED_VAR( eventProgress );
    GEOSX_UNUSED_VAR( domain );
    GEOSX_ERROR( "NOT IMPLEMENTED" );
    return false;
  }

  /// @copydoc geosx::dataRepository::Group::postProcessInput( )
  void postProcessInput() override;
};

} /* namespace */

#endif /* GEOSX_EVENTS_TASKS_TASKBASE_HPP_ */
