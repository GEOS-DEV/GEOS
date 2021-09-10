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
 * @file TasksManager.hpp
 */

#ifndef GEOSX_EVENTS_TASKS_TASKSMANAGER_HPP_
#define GEOSX_EVENTS_TASKS_TASKSMANAGER_HPP_

#include "dataRepository/Group.hpp"
#include "TaskBase.hpp"

namespace geosx
{

/**
 * @class TasksManager
 * @brief A class to manage and execute tasks.
 */
class TasksManager : public dataRepository::Group
{
public:
  /// @copydoc geosx::dataRepository::Group::Group
  TasksManager( string const & name, Group * const parent );
  /// Destructor
  virtual ~TasksManager() override;

  /// @copydoc geosx::dataRepository::Group::createChild
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /// This function is used to expand any catalogs in the data structure
  virtual void expandObjectCatalogs() override;

private:
  TasksManager() = delete;

};

} /* namespace geosx */

#endif
