/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
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
#ifndef TWOMUTEXESLOCK_HPP
#define TWOMUTEXESLOCK_HPP

#include <mutex>

namespace geosx
{

/**
 * @brief Class to handle locks using 2 mutexes
 */
class TwoMutexesLock
{
public:
  /**
   * @brief Construct a dual mutex lock and lock the mutexes.
   * @param mutex1 First mutex of the dual mutex
   * @param mutex2 Second mutex of the dual mutex
   */
  TwoMutexesLock( std::mutex & mutex1, std::mutex & mutex2 ):
    m_islocked( false ), m_mutex1( &mutex1 ), m_mutex2( &mutex2 )
  {
    lock();
  }

  /**
   * @brief Unlock the mutexes and destroy the locks.
   */
  ~TwoMutexesLock()
  {
    unlock();
  }

  /**
   * @brief Lock the two mutexes using std::lock is not already owning lock
   */
  void lock()
  {
    if( m_islocked ) return;
    std::lock( *m_mutex1, *m_mutex2 );
    m_islocked = true;
  }

  /**
   * @brief Unlock the two mutexes is owning them.
   */
  void unlock()
  {
    if( !m_islocked ) return;
    m_mutex1->unlock();
    m_mutex2->unlock();
    m_islocked = false;
  }
private:
  /// Indicate if the mutexes are owned by the lock
  bool m_islocked;
  /// Pointer to the first mutex
  std::mutex *m_mutex1;
  /// Pointer to the second mutex
  std::mutex *m_mutex2;
};
}
#endif // TWOMUTEXESLOCK_HPP
