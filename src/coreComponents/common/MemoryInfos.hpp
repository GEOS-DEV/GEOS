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

#ifndef GEOS_COMMON_MemoryInfos_HPP_
#define GEOS_COMMON_MemoryInfos_HPP_

#include "umpire/util/MemoryResourceTraits.hpp"
#include "common/Logger.hpp"
#include <unistd.h>
#include <iostream>
#if defined( GEOS_USE_CUDA )
#include <cuda.h>
#endif

namespace geos
{

class MemoryInfos
{
public:
  MemoryInfos( umpire::MemoryResourceTraits::resource_type resourceType );

  size_t getTotalMemory() const;
  size_t getAvailableMemory() const;

private:
  size_t m_totalMemory;
  size_t m_availableMemory;
};

}

#endif
