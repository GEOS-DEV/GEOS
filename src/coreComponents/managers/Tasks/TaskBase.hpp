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

#ifndef TOOLBASE_HPP_
#define TOOLBASE_HPP_



#include <string>
#include <limits>

#include "dataRepository/ManagedGroup.hpp"
#include "dataRepository/ExecutableGroup.hpp"
#include "common/DataTypes.hpp"
namespace geosx
{

class TaskBase : public ExecutableGroup
{
public:

  explicit TaskBase( std::string const & name,
                       ManagedGroup * const parent );

  TaskBase( TaskBase && ) = default;

  virtual ~TaskBase() override;

  TaskBase() = delete;
  TaskBase( TaskBase const & ) = delete;
  TaskBase& operator=( TaskBase const & ) = delete;
  TaskBase& operator=( TaskBase&& ) = delete;

  static string CatalogName() { return "TaskBase"; }

  using CatalogInterface = cxx_utilities::CatalogInterface< TaskBase, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

  void PostProcessInput() override;


};

} /* namespace */


#endif /* TOOLBASE_HPP_ */
