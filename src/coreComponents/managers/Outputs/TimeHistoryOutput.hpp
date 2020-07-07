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
#include "TimeHistoryCollector.hpp"
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
                       Group * const parent ):
      OutputBase(name,parent),
      m_collector_paths( ),
      m_format( ),
      m_filename( ),
      m_record_count(0),
      m_io( )
    {
      registerWrapper(viewKeys::timeHistoryOutputTarget, &m_collector_paths)->
        setInputFlag(InputFlags::REQUIRED)->
        setDescription("A list of collectors from which to collect and output time history information.");

      registerWrapper(viewKeys::timeHistoryOutputFilename, &m_filename)->
        setApplyDefaultValue("TimeHistory")->
        setInputFlag(InputFlags::OPTIONAL)->
        setDescription("The filename to which to write time history output.");

      registerWrapper(viewKeys::timeHistoryOutputFormat, &m_format)->
        setApplyDefaultValue("hdf")->
        setInputFlag(InputFlags::OPTIONAL)->
        setDescription("The output file format for time history output.");

      registerWrapper(viewKeys::timeHistoryRestart, &m_record_count)->
        setApplyDefaultValue(0)->
        setInputFlag(InputFlags::FALSE)->
        setRestartFlags(RestartFlags::WRITE_AND_READ)->
        setDescription("The current history record to be written, on restart from an earlier time allows use to remove invalid future history.");

    }

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
    virtual void InitializePostSubGroups( Group * const group ) override
    {
      // handle recreating the file up front to avoid issues with creating/adding multiple datasets to the file on seperate file accesses
      {
        // if the record count is zero, this isn't a restart, so delete any existing data
        // this is supposed to be filetype agnostic, so using system interfaces to delete the file would be preferable
        HDFFile( m_filename, (m_record_count == 0) );
      }
      for( auto collector_path : m_collector_paths )
      {
        Group * tmp = this->GetGroupByPath( collector_path );
        HistoryCollection * collector = Group::group_cast< HistoryCollection* >( tmp );
        GEOSX_ERROR_IF( collector == nullptr, "The target of a time history output event must be a collector! " << collector_path );
        // todo: switch based on m_format, always hdf for now
        HistoryMetadata metadata = collector->GetMetadata( group );
        int rank = MpiWrapper::Comm_rank();
        if( rank == 0 )
        {
          HistoryMetadata time_metadata = collector->GetTimeMetadata( );
          time_metadata.setName(metadata.getName() + string(" ") + time_metadata.getName());
          m_io.emplace_back( std::make_pair( std::make_unique<HDFHistIO>( m_filename, metadata, m_record_count ),
                                             std::make_unique<HDFHistIO>( m_filename, time_metadata, m_record_count, 4, MPI_COMM_SELF ) ) );
          collector->RegisterTimeBufferCall([this]() { return this->m_io.back().second->GetBufferHead( ); });
          m_io.back().second->Init( ( m_record_count > 0 ) );
        }
        else
        {
          m_io.emplace_back( std::make_pair( std::make_unique<HDFHistIO>( m_filename, metadata, m_record_count ), std::unique_ptr<HDFHistIO>(nullptr) ) );
        }
        collector->RegisterBufferCall([this]() { return this->m_io.back().first->GetBufferHead( ); });
        m_io.back().first->Init( ( m_record_count > 0 ) );

        if ( m_record_count == 0 )
        {
          // do any 1-time metadata output
          globalIndex global_rank_offset = m_io.back().first->GetRankOffset( );
          localIndex meta_collector_count = collector->GetNumMetaCollectors( );
          Group * domain_group = dynamicCast< Group * >( dynamicCast< ProblemManager * >( group )->getDomainPartition( ) );
          for( localIndex meta_idx = 0; meta_idx < meta_collector_count; ++meta_idx )
          {
            std::unique_ptr<HistoryCollection> meta_collector = collector->GetMetaCollector( group, meta_idx, global_rank_offset );
            HistoryMetadata meta_metadata = meta_collector->GetMetadata( group );
            meta_metadata.setName(metadata.getName() + " " +meta_metadata.getName());
            std::unique_ptr<HDFHistIO> meta_io = std::make_unique<HDFHistIO>( m_filename, meta_metadata, 0, 1 );
            meta_collector->RegisterBufferCall([&meta_io] () { return meta_io->GetBufferHead( ); });
            meta_io->Init( false );
            meta_collector->Execute( 0.0, 0.0, 0, 0, 0, domain_group );
            meta_io->Write( );
          }
        }
      }
    }

    /**
     * @brief Writes out a time history file.
     * @copydoc EventBase::Execute()
     */
    virtual void Execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                          real64 const GEOSX_UNUSED_PARAM( dt ),
                          integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                          integer const GEOSX_UNUSED_PARAM( eventCounter ),
                          real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                          dataRepository::Group * GEOSX_UNUSED_PARAM( domain ) ) override
    {
      for ( auto & io_pair : m_io )
      {
        m_record_count += io_pair.first->GetBufferedCount( );
        io_pair.first->Write( );
        if( io_pair.second )
        {
          io_pair.second->Write( );
        }
      }
    }

    /**
     * @brief Writes out a time history file at the end of the simulation.
     * @copydoc ExecutableGroup::Cleanup()
     */
    virtual void Cleanup( real64 const time_n,
                          integer const cycleNumber,
                          integer const eventCounter,
                          real64 const eventProgress,
                          dataRepository::Group * domain ) override
    {
      Execute(time_n,0.0,cycleNumber,eventCounter,eventProgress,domain);
      // remove any unused trailing space reserved to write additional histories
      for( auto & io_pair : m_io )
      {
        io_pair.first->CompressInFile();
        if( io_pair.second )
        {
          io_pair.second->CompressInFile();
        }
      }
    }

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
