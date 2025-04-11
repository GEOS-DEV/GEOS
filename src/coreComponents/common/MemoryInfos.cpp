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

#include "MemoryInfos.hpp"

#include "MpiWrapper.hpp"
#if defined( GEOS_USE_CALIPER ) and defined( GEOS_USE_ADIAK )
#include <adiak.hpp>
#endif

namespace geos
{

MemoryInfos::MemoryInfos( umpire::MemoryResourceTraits::resource_type resourceType ):
  m_totalMemory( 0 ),
  m_availableMemory( 0 ),
  m_physicalMemoryHandled( 1 )
{
  switch( resourceType )
  {
    case umpire::MemoryResourceTraits::resource_type::host:
    case umpire::MemoryResourceTraits::resource_type::pinned:
      #if defined( _SC_PHYS_PAGES ) && defined( _SC_PAGESIZE )
      m_totalMemory = sysconf( _SC_PHYS_PAGES ) * sysconf( _SC_PAGESIZE );
      #if defined(_SC_AVPHYS_PAGES)
      m_availableMemory = sysconf( _SC_AVPHYS_PAGES ) * sysconf( _SC_PAGESIZE );
      #else
      GEOS_WARNING( "Unknown device avaialable memory size getter for this system." );
      m_availableMemory = 0;
      #endif
      #else
      GEOS_WARNING( "Unknown device physical memory size getter for this compiler." );
      m_physicalMemoryHandled = 0;
      #endif
      break;
    case umpire::MemoryResourceTraits::resource_type::device:
    case umpire::MemoryResourceTraits::resource_type::device_const:
    case umpire::MemoryResourceTraits::resource_type::um:
      #if defined( GEOS_USE_CUDA )
      cudaMemGetInfo( &m_availableMemory, &m_totalMemory );
      #else
      GEOS_WARNING( "Unknown device physical memory size getter for this compiler." );
      m_physicalMemoryHandled = 0;
      #endif
      break;
    default:
      GEOS_WARNING( "Physical memory lookup not implemented" );
      m_physicalMemoryHandled = 0;
      break;
  }
}

size_t MemoryInfos::getTotalMemory() const
{
  return m_totalMemory;
}

size_t MemoryInfos::getAvailableMemory() const
{
  return m_availableMemory;
}

bool MemoryInfos::isPhysicalMemoryHandled() const
{
  return m_physicalMemoryHandled;
}

MemoryLogging::MemoryLogging():
  m_umpireStatsLogReport( true ),
  m_umpireStatsCsvReport( false ),
  m_umpireStatsCsvReportFilename( "./memoryStats.csv" ),
  m_currentCycle( 0 ),
  m_currentTime( 0.0 )
{
  TableLayout memoryStatLayout = { "Cycle",
                                   "Time",
                                   "Umpire Memory Pool\n(reserved / % over total)",
                                   "Min over ranks",
                                   "Max over ranks",
                                   "Avg over ranks",
                                   "Sum over ranks" };
  m_memoryStatLogFormatter = std::make_unique< TableTextFormatter >( memoryStatLayout );

  // We do not need to output the cycle & time in the log table
  memoryStatLayout.getColumns()[0].setVisibility( false );
  memoryStatLayout.getColumns()[1].setVisibility( false );
  m_memoryStatCsvFormatter = std::make_unique< TableCSVFormatter >( memoryStatLayout );
}

MemoryLogging & MemoryLogging::getInstance()
{
  static MemoryLogging instance;
  return instance;
}

void MemoryLogging::enableUmpireStatsCsvReport( bool enable )
{
  m_umpireStatsCsvReport = enable;
  if( enable )
  {
    std::ofstream csvFile{ m_umpireStatsCsvReportFilename };
    m_memoryStatCsvFormatter->headerToStream( csvFile );
  }
}

void MemoryLogging::memoryStatsReport() const
{
  if( !m_umpireStatsLogReport && !m_umpireStatsCsvReport )
    return;

  umpire::ResourceManager & rm = umpire::ResourceManager::getInstance();
  integer size;
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  size_t nbRank = (std::size_t)size;
  // Get a list of all the allocators and sort it so that it's in the same order on each rank.
  std::vector< string > allocatorNames = rm.getAllocatorNames();
  std::sort( allocatorNames.begin(), allocatorNames.end() );

  // If each rank doesn't have the same number of allocators you can't aggregate them.
  std::size_t const numAllocators = allocatorNames.size();
  std::size_t const minNumAllocators = MpiWrapper::min( numAllocators );

  if( numAllocators != minNumAllocators )
  {
    GEOS_WARNING( "Not all ranks have created the same number of umpire allocators, cannot compute high water marks." );
    return;
  }

  // Loop over the allocators to collect stats.
  TableData tableData;
  for( string const & allocatorName : allocatorNames )
  {
    // Skip umpire internal allocators.
    if( allocatorName.rfind( "__umpire_internal", 0 ) == 0 )
      continue;

    static constexpr int MAX_NAME_LENGTH = 100;
    GEOS_ERROR_IF_GT( allocatorName.size(), MAX_NAME_LENGTH );
    string allocatorNameFixedSize = allocatorName;
    allocatorNameFixedSize.resize( MAX_NAME_LENGTH, '\0' );
    string allocatorNameMinChars = string( MAX_NAME_LENGTH, '\0' );

    // Make sure that each rank is looking at the same allocator.
    MpiWrapper::allReduce( allocatorNameFixedSize, allocatorNameMinChars, MpiWrapper::Reduction::Min, MPI_COMM_GEOS );
    if( allocatorNameFixedSize != allocatorNameMinChars )
    {
      GEOS_WARNING( "Not all ranks have an allocator named " << allocatorNameFixedSize << ", cannot compute high water mark." );
      continue;
    }

    umpire::Allocator allocator = rm.getAllocator( allocatorName );
    umpire::strategy::AllocationStrategy const * allocationStrategy = allocator.getAllocationStrategy();
    umpire::MemoryResourceTraits const traits = allocationStrategy->getTraits();
    umpire::MemoryResourceTraits::resource_type resourceType = traits.resource;
    MemoryInfos const memInfos( resourceType );

    if( !memInfos.isPhysicalMemoryHandled() )
    {
      continue;
    }

    // Get the total number of bytes allocated with this allocator across ranks.
    // This is a little redundant since
    std::size_t const mark = allocator.getHighWatermark();
    std::size_t const minMark = MpiWrapper::min( mark );
    std::size_t const maxMark = MpiWrapper::max( mark );
    std::size_t const sumMark = MpiWrapper::sum( mark );

    string percentage;
    if( memInfos.getTotalMemory() == 0 )
    {
      percentage = 0.0;
      GEOS_WARNING( "umpire memory percentage could not be resolved" );
    }
    else
    {
      percentage = GEOS_FMT( "({:.1f}%)", ( 100.0f * (float)mark ) / (float)memInfos.getTotalMemory() );
    }

    string const minMarkValue = GEOS_FMT( "{} {:>8}",
                                          LvArray::system::calculateSize( minMark ), percentage );
    string const maxMarkValue = GEOS_FMT( "{} {:>8}",
                                          LvArray::system::calculateSize( maxMark ), percentage );
    string const avgMarkValue = GEOS_FMT( "{} {:>8}",
                                          LvArray::system::calculateSize( sumMark / nbRank ), percentage );
    string const sumMarkValue = GEOS_FMT( "{} {:>8}",
                                          LvArray::system::calculateSize( sumMark ), percentage );

    tableData.addRow( m_currentCycle,
                      m_currentTime,
                      allocatorName,
                      minMarkValue,
                      maxMarkValue,
                      avgMarkValue,
                      sumMarkValue );

#if defined( GEOS_USE_CALIPER )and defined( GEOS_USE_ADIAK )
    pushStatsIntoAdiak( allocatorName + " sum across ranks", mark );
    pushStatsIntoAdiak( allocatorName + " rank max", mark );
#endif
  }

  { // output statistics
    if( m_umpireStatsLogReport )
    {
      GEOS_LOG_RANK_0( m_memoryStatLogFormatter->toString( tableData ) );
    }

    if( m_umpireStatsCsvReport )
    {
      std::ofstream csvFile{ m_umpireStatsCsvReportFilename };
      m_memoryStatCsvFormatter->dataToStream( csvFile, tableData );
    }
  }
}

} /* namespace geos */
