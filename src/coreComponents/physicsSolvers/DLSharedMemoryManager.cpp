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

#include "DLSharedMemoryManager.hpp"
#include "PhysicsSolverManager.hpp"

#include "physicsSolvers/LogLevelsInfo.hpp"
#include "common/format/LogPart.hpp"
#include "common/TimingMacros.hpp"
#include "linearAlgebra/solvers/KrylovSolver.hpp"
#include "mesh/DomainPartition.hpp"
#include "math/interpolation/Interpolation.hpp"
#include "common/Timer.hpp"
#include "common/Units.hpp"
#include "dataRepository/LogLevelsInfo.hpp"

#if defined(GEOS_USE_PYGEOSX)
#include "python/PySolverType.hpp"
#endif

namespace geos
{

using namespace dataRepository;

DLSharedMemoryManager::DLSharedMemoryManager( string const & name,
                                              Group * const parent )
  : dataRepository::Group( name, parent )
{}

DLSharedMemoryManager::~DLSharedMemoryManager() = default;


sem_t * DLSharedMemoryManager::openASemaphore( const std::string & mem_name )
{
  // TODO: check access permissions and handle errors
  sem_t * ptr = sem_open( mem_name.c_str(), O_CREAT | O_RDWR, 0666, 0 );

  return ptr;
}


#if defined(GEOS_USE_PYGEOSX)
PyTypeObject *DLSharedMemoryManager::getPythonType() const
{
  return python::getPySolverType();
}
#endif

} // namespace geos
