/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MemoryStatsOutput.hpp
 */

#ifndef GEOS_FILEIO_OUTPUTS_UMPIRESTATSOUTPUT_HPP_
#define GEOS_FILEIO_OUTPUTS_UMPIRESTATSOUTPUT_HPP_

#include "OutputBase.hpp"


namespace geos
{

/**
 * @brief An output Task to control the logging of memory statistics.
 */
class MemoryStatsOutput final : public OutputBase
{
public:

  /// @copydoc geos::dataRepository::Group::Group(string const & name, Group * const parent)
  MemoryStatsOutput( string const & name, Group * const parent );

  /// Destructor
  virtual ~MemoryStatsOutput() override
  {}

  /**
   * @brief Catalog name interface
   * @return This type's catalog name
   */
  static string catalogName() { return "MemoryStats"; }

  /**
   * @brief Writes out a set of vtk files.
   * @copydoc EventBase::execute()
   */
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /// @cond DO_NOT_DOCUMENT
  struct viewKeysStruct : OutputBase::viewKeysStruct
  {
    static constexpr auto writeCSV = "writeCSV";
  };
  /// @endcond

protected:

  void postInputInitialization() override;

  /**
   * @copydoc OutputBase::getTimerCategory
   */
  logInfo::OutputTimerBase const & getTimerCategory() const override;

private:

  // Flag to enable writing CSV output
  integer m_writeCSV;

};

} /* namespace geos */

#endif /* GEOS_FILEIO_OUTPUTS_UMPIRESTATSOUTPUT_HPP_ */
