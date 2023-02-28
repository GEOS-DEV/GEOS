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
#ifndef FOURMUTEXESLOCK_HPP
#define FOURMUTEXESLOCK_HPP

#include <mutex>

namespace geosx
{

/**
 * @brief Class to handle locks using 4 mutexes
 */
class FourMutexesLock
{
public:
  /**
   * @brief Construct a four mutex lock and lock the mutexes.
   * @param mutex1 First mutex of the dual mutex
   * @param mutex2 Second mutex of the dual mutex
   * @param mutex3 Third mutex of the dual mutex
   * @param mutex4 Fourth mutex of the dual mutex
   */
  FourMutexesLock( std::mutex & mutex1, std::mutex & mutex2, std::mutex & mutex3, std::mutex & mutex4 ):
    m_islocked( false ), m_mutex1( &mutex1 ), m_mutex2( &mutex2 ), m_mutex3( &mutex3 ), m_mutex4( &mutex4 )
  {
    lock();
  }

  /**
   * @brief Unlock the mutexes and destroy the locks.
   */
  ~FourMutexesLock()
  {
    unlock();
  }

  /**
   * @brief Lock the two mutexes using std::lock is not already owning lock
   */
  void lock()
  {
    if( m_islocked ) return;
    std::lock( *m_mutex1, *m_mutex2, *m_mutex3, *m_mutex4 );
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
    m_mutex3->unlock();
    m_mutex4->unlock();
    m_islocked = false;
  }
private:
  /// Indicate if the mutexes are owned by the lock
  bool m_islocked;
  /// Pointer to the first mutex
  std::mutex *m_mutex1;
  /// Pointer to the second mutex
  std::mutex *m_mutex2;
  /// Pointer to the third mutex
  std::mutex *m_mutex3;
  /// Pointer to the fourth mutex
  std::mutex *m_mutex4;
};
}
#endif // FOURMUTEXESLOCK
