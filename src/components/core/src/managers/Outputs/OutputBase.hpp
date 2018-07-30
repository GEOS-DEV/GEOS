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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
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

class OutputBase : public ExecutableGroup
{
public:
  explicit OutputBase( std::string const & name,
                       ManagedGroup * const parent );

  virtual ~OutputBase() override;

  static string CatalogName() { return "OutputBase"; }

  OutputBase() = default;
  OutputBase( OutputBase const & ) = default;
  OutputBase( OutputBase &&) = default;
  OutputBase& operator=( OutputBase const & ) = default;
  OutputBase& operator=( OutputBase&& ) = default;
    
  virtual void FillDocumentationNode() override;

  virtual void Initialize( ManagedGroup * const group ) override;

  virtual void SetupDirectoryStructure();

  using CatalogInterface = cxx_utilities::CatalogInterface< OutputBase, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

  struct viewKeyStruct
  {
    dataRepository::ViewKey slaveDirectory = { "slaveDirectory" };
    dataRepository::ViewKey parallelThreads = { "parallelThreads" };
  } viewKeys;

};


} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_OUTPUTBASE_HPP_ */
