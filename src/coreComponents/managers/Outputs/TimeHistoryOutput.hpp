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
 * @file TimeHistoryOutput.hpp
 */

#ifndef GEOSX_TIME_HISTORY_OUTPUT_HPP_
#define GEOSX_TIME_HISTORY_OUTPUT_HPP_

#include "OutputBase.hpp"
#include "managers/TimeHistory/TimeHistoryCollection.hpp"
#include "managers/TimeHistory/HistoryIO.hpp"
#include "fileIO/hdf/HDFFile.hpp"

#include "LvArray/src/Array.hpp" // just for collector

namespace geosx
{

/**
 * @class TimeHistoryOutput
 *
 * A class for creating time history output
 */
  class TimeHistoryOutput : public OutputBase
  {
  public:
    /// @copydoc geosx::dataRepository::Group::Group(std::string const & name, Group * const parent)
    TimeHistoryOutput( string const & name,
                       Group * const parent );

    /// Destructor
    virtual ~TimeHistoryOutput() override
    { }

    /**
     * @brief Catalog name interface
     * @return This type's catalog name
     */
    static string CatalogName() { return "TimeHistory"; }

    /**
     * @brief Perform initalization after all subgroups have been initialized.
     *        Check for existing files and data spaces/sets on restart, else create the
     *        file and data spaces/sets and output any set index metatata for the time
     *        history.
     * @param group The problem manager cast to a group.
     * @note There are operations in this function that are collective on the GEOSX comm.
     */
    virtual void InitializePostSubGroups( Group * const group ) override;

    /**
     * @brief Writes out a time history file.
     * @copydoc EventBase::Execute()
     */
    virtual void Execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                          real64 const GEOSX_UNUSED_PARAM( dt ),
                          integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                          integer const GEOSX_UNUSED_PARAM( eventCounter ),
                          real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                          dataRepository::Group * GEOSX_UNUSED_PARAM( domain ) ) override;
    /**
     * @brief Writes out a time history file at the end of the simulation.
     * @copydoc ExecutableGroup::Cleanup()
     */
    virtual void Cleanup( real64 const time_n,
                          integer const cycleNumber,
                          integer const eventCounter,
                          real64 const eventProgress,
                          dataRepository::Group * domain ) override;

    /// @cond DO_NOT_DOCUMENT
    struct viewKeys
    {
      static constexpr auto timeHistoryOutputTarget = "sources";
      static constexpr auto timeHistoryOutputFilename = "filename";
      static constexpr auto timeHistoryOutputFormat = "format";
      static constexpr auto timeHistoryRestart = "restart";
    } timeHistoryOutputViewKeys;
    /// @endcond

  private:
    string_array m_collector_paths;
    string m_format;
    string m_filename;
    integer m_record_count;
    std::vector< std::pair< std::unique_ptr< BufferedHistoryIO >, std::unique_ptr< BufferedHistoryIO > > > m_io;
  };
}

#endif
