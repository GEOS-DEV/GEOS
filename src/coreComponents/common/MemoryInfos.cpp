/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-LiCense-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "MemoryInfos.hpp"

namespace geos
{
MemoryInfos::MemoryInfos( umpire::MemoryResourceTraits::resource_type resourceType )
{
  switch( resourceType )
  {
    case umpire::MemoryResourceTraits::resource_type::host:
      GEOS_LOG_RANK_0( "host" );
      #if defined( _SC_PHYS_PAGES ) && defined( _SC_PAGESIZE )
      m_totalMemory = sysconf( _SC_PHYS_PAGES ) * sysconf( _SC_PAGESIZE );
      m_availableMemory = sysconf( _SC_AVPHYS_PAGES ) * sysconf( _SC_PAGESIZE );
      #else
      GEOS_ERROR( "Unsupported resource type : host" );
      #endif
      break;
    case umpire::MemoryResourceTraits::resource_type::device:
      GEOS_LOG_RANK_0( "device" );
      #if defined( GEOS_USE_CUDA )
      cudaMemGetInfo( &m_availableMemory, &m_totalMemory );
      #else
      GEOS_ERROR( "Unsupported resource type : device" );
      #endif
      break;
    case umpire::MemoryResourceTraits::resource_type::device_const:
      GEOS_LOG_RANK_0( "device_const" );
      m_totalMemory = 0;
      m_availableMemory = 0;
      //to implement
      break;
    case umpire::MemoryResourceTraits::resource_type::pinned:
      GEOS_LOG_RANK_0( "pinned" );
      m_totalMemory = 0;
      m_availableMemory = 0;
      //regarder a chaque appel de la pinned memory, maintenir un compteur
      break;
    default:
      GEOS_LOG_RANK_0( "unknown" );
      GEOS_ERROR( "Unsupported resource type : device" );
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

}
