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

/**
 * @file OutputManager.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_OUTPUTMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_OUTPUTMANAGER_HPP_

#include "dataRepository/Group.hpp"


namespace geosx
{
namespace dataRepository
{
namespace keys
{
}
}

/**
 * @class OutputManager
 *
 * An class for managing output types
 */
class OutputManager : public dataRepository::Group
{
public:
  /// Main constructor
  OutputManager( std::string const & name,
                 Group * const parent );

  /// Destructor
  virtual ~OutputManager() override;

  /// Catalog name interface
  virtual Group * CreateChild( string const & childKey, string const & childName ) override;

  /// This function is used to expand any catalogs in the data structure
  virtual void ExpandObjectCatalogs() override;

  struct viewKeyStruct
  {
    dataRepository::ViewKey time = { "time" };
  } viewKeys;

};


} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_OUTPUTMANAGER_HPP_ */
