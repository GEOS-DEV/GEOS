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
 * @file DLSharedMemoryManager.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_DLSHAREDMEMORYMANAGER_HPP_
#define GEOS_PHYSICSSOLVERS_DLSHAREDMEMORYMANAGER_HPP_

#include "codingUtilities/traits.hpp"
#include "common/DataTypes.hpp"
#include "common/format/LogPart.hpp"
#include "dataRepository/ExecutableGroup.hpp"
#include "dataRepository/RestartFlags.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/utilities/LinearSolverResult.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "mesh/MeshBody.hpp"
#include "physicsSolvers/NonlinearSolverParameters.hpp"
#include "physicsSolvers/LinearSolverParameters.hpp"
#include "physicsSolvers/SolverStatistics.hpp"
#include "physicsSolvers/PhysicsSolverBase.hpp"

// Shared memory includes
#include <sys/mman.h>  // Provides functions for memory management (e.g., mmap, munmap).
#include <fcntl.h>     // Provides file control options (e.g., open, O_CREAT).
#include <unistd.h>    // Provides access to the POSIX operating system API (e.g., read, write, close, fork, exec).
#include <cstring>     // Provides functions for string manipulation (e.g., memcpy, memmove, memset, strlen).
#include <semaphore.h> // Provides POSIX semaphore functionality for inter-process synchronization (e.g., sem_open, sem_wait, sem_post).

#include <limits>

namespace geos
{

/**
 * @class DLSharedMemoryManager
 * @brief a class for intitalizing and managing shared memory in DL simulations
 *
 */
class DLSharedMemoryManager : public dataRepository::Group
{
public:
  /**
   * @brief Constructor for DLSharedMemoryManager
   * @param name the name of this instantiation of DLSharedMemoryManager
   * @param parent the parent group of this instantiation of DLSharedMemoryManager
   */
  explicit DLSharedMemoryManager( string const & name,
                                  Group * const parent );

  /**
   * @brief Move constructor for DLSharedMemoryManager
   */
  DLSharedMemoryManager( DLSharedMemoryManager && ) = default;

  /**
   * @brief Destructor for DLSharedMemoryManager
   */
  virtual ~DLSharedMemoryManager() override;

  /**
   * @brief Deleted constructor
   */
  DLSharedMemoryManager() = delete;

  /**
   * @brief Deleted copy constructor
   */
  DLSharedMemoryManager( DLSharedMemoryManager const & ) = delete;

  /**
   * @brief Deleted copy assignment operator
   */
  DLSharedMemoryManager & operator=( DLSharedMemoryManager const & ) = delete;

  /**
   * @brief Deleted move assignment operator
   */
  DLSharedMemoryManager & operator=( DLSharedMemoryManager && ) = delete;

  // Function to initialize a shared memory for a struct
  // As a template, the full definition will be in the header file
  template< typename T >
  T * initializeASharedMemory( const std::string & mem_name, size_t size )
  {
    // TODO: check access permissions and handling errors
    int shm_fd = shm_open( mem_name.c_str(), O_CREAT | O_RDWR, 0666 );
    if( shm_fd == -1 )
    {
      std::cerr << "Error opening shared memory object: " << mem_name << std::endl;
      return nullptr;
    }

    if( ftruncate( shm_fd, size ) == -1 )
    {
      std::cerr << "Error in ftruncate for shared memory: " << mem_name << std::endl;
      return nullptr;
    }

    void *ptr = mmap( 0, size, PROT_READ | PROT_WRITE, MAP_SHARED, shm_fd, 0 );
    if( ptr == MAP_FAILED )
    {
      std::cerr << "Error mapping shared memory: " << mem_name << std::endl;
      return nullptr;
    }

    return static_cast< T * >(ptr);
  }

  sem_t *openASemaphore( const std::string & mem_name );

#if defined(GEOS_USE_PYGEOSX)
  PyTypeObject * getPythonType() const;
#endif

protected:
private:
};

} // namespace geos

#endif /* GEOS_PHYSICSSOLVERS_DLSHAREDMEMORYMANAGER_HPP_ */
