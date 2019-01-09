/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
 * @file OutputBase.hpp
 */
#ifndef SRC_COMPONENTS_CORE_SRC_OUTPUTBASE_HPP_
#define SRC_COMPONENTS_CORE_SRC_OUTPUTBASE_HPP_

#include "dataRepository/ManagedGroup.hpp"
#include "dataRepository/ExecutableGroup.hpp"


namespace geosx
{

/**
 * @class OutputBase
 *
 * A base class for output types
 */
class OutputBase : public ExecutableGroup
{
public:
  /// Main constructor
  explicit OutputBase( std::string const & name,
                       ManagedGroup * const parent );

  /// Destructor
  virtual ~OutputBase() override;

  /// Catalog name interface
  static string CatalogName() { return "OutputBase"; }

  OutputBase() = default;
  OutputBase( OutputBase const & ) = default;
  OutputBase( OutputBase &&) = default;
  OutputBase& operator=( OutputBase const & ) = default;
  OutputBase& operator=( OutputBase&& ) = default;

  /// Method for setting up output directories
  virtual void SetupDirectoryStructure();

  /// Catalog interface
  using CatalogInterface = cxx_utilities::CatalogInterface< OutputBase, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

  struct viewKeysStruct
  {
    static constexpr auto slaveDirectoryString = "slaveDirectory";
    static constexpr auto parallelThreadsString = "parallelThreads";
  } outputBaseViewKeys;

  string slaveDirectory() const { return m_slaveDirectory; }
  integer parallelThreads() const { return m_parallelThreads; }

protected:
  virtual void InitializePreSubGroups( ManagedGroup * const group ) override;

private:
  string m_slaveDirectory;
  integer m_parallelThreads;

};


} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_OUTPUTBASE_HPP_ */
