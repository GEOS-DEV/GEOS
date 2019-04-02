/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef SRC_COMPONENTS_CORE_SRC_TOOLSMANAGER_TOOLSMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_TOOLSMANAGER_TOOLSMANAGER_HPP_

#include "dataRepository/ManagedGroup.hpp"

namespace pugi
{
class xml_node;
}

namespace geosx
{

class TasksManager : public dataRepository::ManagedGroup
{
public:
  TasksManager( std::string const & name,
                ManagedGroup * const parent );

  virtual ~TasksManager() override;

  virtual ManagedGroup * CreateChild( string const & childKey, string const & childName ) override;

  /// This function is used to expand any catalogs in the data structure
  virtual void ExpandObjectCatalogs() override;

  /*
  struct viewKeyStruct
  {
      constexpr static auto gravityVectorString = "gravityVector";
      constexpr static auto blockSystemRepositoryString = "blockSystemRepository";
  } viewKeys;
  */

private:
  TasksManager() = delete;

};

} /* namespace geosx */

#endif
