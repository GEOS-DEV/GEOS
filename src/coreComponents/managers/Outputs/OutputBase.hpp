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
 * @file OutputBase.hpp
 */
#ifndef GEOSX_MANAGERS_OUTPUTS_OUTPUTBASE_HPP_
#define GEOSX_MANAGERS_OUTPUTS_OUTPUTBASE_HPP_

#include "dataRepository/Group.hpp"
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
                       Group * const parent );

  /// Destructor
  virtual ~OutputBase() override;

  /// Catalog name interface
  static string CatalogName() { return "OutputBase"; }

  /// Method for setting up output directories
  virtual void SetupDirectoryStructure();

  /// Catalog interface
  using CatalogInterface = dataRepository::CatalogInterface< OutputBase, std::string const &, Group * const >;
  static CatalogInterface::CatalogType & GetCatalog();

  struct viewKeysStruct
  {
    static constexpr auto slaveDirectoryString = "slaveDirectory";
    static constexpr auto parallelThreadsString = "parallelThreads";
  } outputBaseViewKeys;

  string slaveDirectory() const { return m_slaveDirectory; }
  integer parallelThreads() const { return m_parallelThreads; }

protected:
  virtual void InitializePreSubGroups( Group * const group ) override;

private:
  string m_slaveDirectory;
  integer m_parallelThreads;

};


} /* namespace geosx */

#endif /* GEOSX_MANAGERS_OUTPUTS_OUTPUTBASE_HPP_ */
