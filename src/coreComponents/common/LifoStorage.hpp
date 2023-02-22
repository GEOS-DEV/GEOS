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
#ifndef LIFOSTORAGE_HPP
#define LIFOSTORAGE_HPP

#include <deque>
#include <future>
#include <mutex>
#include <condition_variable>
#include <camp/camp.hpp>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <algorithm>

#include "common/FixedSizeDeque.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "common/TimingMacros.hpp"

#ifdef LIFO_DISABLE_CALIPER
#define LIFO_MARK_FUNCTION
#define LIFO_MARK_SCOPE(a)
#define LIFO_LOG_RANK(a) std::cerr << a << std::endl;
#else
#define LIFO_MARK_FUNCTION GEOSX_MARK_FUNCTION
#define LIFO_MARK_SCOPE(a) GEOSX_MARK_SCOPE(a)
#define LIFO_LOG_RANK(a) GEOSX_LOG_RANK(a)
#endif

namespace geosx
{

/**
 * @brief Class to handle locks using 2 mutexes
 */
class TwoMutexLock
{
public:
  /**
   * @brief Construct a dual mutex lock and lock the mutexes.
   * @param mutex1 First mutex of the dual mutex
   * @param mutex2 Second mutex of the dual mutex
   */
  TwoMutexLock( std::mutex & mutex1, std::mutex & mutex2 ):
    m_islocked( false ), m_mutex1( &mutex1 ), m_mutex2( &mutex2 )
  {
    lock();
  }

  /**
   * @brief Unlock the mutexes and destroy the locks.
   */
  ~TwoMutexLock()
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

/**
 * @brief Class to handle locks using 4 mutexes
 */
class FourMutexLock
{
public:
  /**
   * @brief Construct a four mutex lock and lock the mutexes.
   * @param mutex1 First mutex of the dual mutex
   * @param mutex2 Second mutex of the dual mutex
   * @param mutex3 Third mutex of the dual mutex
   * @param mutex4 Fourth mutex of the dual mutex
   */
  FourMutexLock( std::mutex & mutex1, std::mutex & mutex2, std::mutex & mutex3, std::mutex & mutex4 ):
    m_islocked( false ), m_mutex1( &mutex1 ), m_mutex2( &mutex2 ), m_mutex3( &mutex3 ), m_mutex4( &mutex4 )
  {
    lock();
  }

  /**
   * @brief Unlock the mutexes and destroy the locks.
   */
  ~FourMutexLock()
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

/// Associate mutexes with the fixedSizeDeque
template< typename T >
class FixedSizeDequeAndMutexes : public FixedSizeDeque< T, int >
{
public:
  /// Mutex to protect access to the front
  std::mutex m_frontMutex;
  /// Mutex to protect access to the back
  std::mutex m_backMutex;
  /// Mutex to prevent two simulteaneous pop (can be an issue for last one)
  std::mutex m_popMutex;
  /// Mutex to prevent two simulteaneous emplace (can be an issue for last one)
  std::mutex m_emplaceMutex;
  /// Condition used to notify when queue is not full
  std::condition_variable_any m_notFullCond;
  /// Condition used to notify when queue is not empty
  std::condition_variable_any m_notEmptyCond;

  /**
   * Create a fixed size double ended queue with associated mutexes and condition variables.
   *
   * @param maxEntries     Maximum number of array to store in the queue.
   * @param valuesPerEntry Number of values in each array of the deque.
   * @param space          Space used to store que queue.
   */
  FixedSizeDequeAndMutexes( int maxEntries, int valuesPerEntry, LvArray::MemorySpace space ): FixedSizeDeque< T, int >( maxEntries, valuesPerEntry, space,
#ifdef GEOSX_USE_CUDA
                                                                                                                        camp::resources::Resource{ camp::resources::Cuda{} }
#else
                                                                                                                        camp::resources::Resource{ camp::resources::Host{} }
#endif
                                                                                                                        ) {}

  /**
   * Emplace on front from array with locks.
   *
   * @param array The array to emplace
   * @returns Event to sync with the memcpy.
   */
  camp::resources::Event emplaceFront( arrayView1d< T > array )
  {
    LIFO_MARK_FUNCTION;
    camp::resources::Event e;
    {
      TwoMutexLock lock( m_emplaceMutex, m_frontMutex );
      {
        LIFO_MARK_SCOPE( waitingForBuffer );
        m_notFullCond.wait( lock, [ this ]  { return !this->full(); } );
      }
      {
        LIFO_MARK_SCOPE( copy );
        e = FixedSizeDeque< T, int >::emplace_front( array.toSliceConst() );
      }
    }
    m_notEmptyCond.notify_all();
    return e;
  }

  /**
   * Pop from front to array with locks
   *
   * @param array The array to copy data into
   * @returns Event to sync with the memcpy.
   */
  camp::resources::Event popFront( arrayView1d< T > array )
  {
    LIFO_MARK_FUNCTION;
    camp::resources::Event e;
    {
      TwoMutexLock lock( m_popMutex, m_frontMutex );
      {
        LIFO_MARK_SCOPE( waitingForBuffer );
        m_notEmptyCond.wait( lock, [ this ]  { return !this->empty(); } );
      }
      // deadlock can occur if frontMutex is taken after an
      // emplaceMutex (inside pushAsync) but this is prevented by the
      // pushWait() in popAsync.
      {
        LIFO_MARK_SCOPE( copy );
        camp::resources::Resource r = this->getStream();
        e = LvArray::memcpy( r, array.toSlice(), this->front() );
        this->pop_front();
      }
    }
    m_notFullCond.notify_all();
    return e;
  }

  /**
   * Emplace front from back of given queue
   *
   * @param q2 The queue to copy data from.
   */
  void emplaceFrontFromBack( FixedSizeDequeAndMutexes< T > & q2 )
  {
    LIFO_MARK_FUNCTION;
    {
      FourMutexLock lock( m_emplaceMutex, q2.m_popMutex, m_frontMutex, q2.m_backMutex );
      while( this->full() || q2.empty() )
      {
        {
          LIFO_MARK_SCOPE( WaitForBufferToEmplace );
          m_notFullCond.wait( lock, [ this ]  { return !this->full(); } );
        }
        {
          LIFO_MARK_SCOPE( WaitForBufferToPop );
          q2.m_notEmptyCond.wait( lock, [ &q2 ] { return !q2.empty(); } );
        }
      }
      LIFO_MARK_SCOPE( Transfert );
      this->emplace_front( q2.back() ).wait();
      q2.pop_back();
    }
    q2.m_notFullCond.notify_all();
    m_notEmptyCond.notify_all();
  }

  /**
   * Emplace back from front of given queue
   *
   * @param q2 The queue to copy data from.
   */
  void emplaceBackFromFront( FixedSizeDequeAndMutexes< T > & q2 )
  {
    LIFO_MARK_FUNCTION;
    {
      FourMutexLock lock( m_emplaceMutex, q2.m_popMutex, m_backMutex, q2.m_frontMutex );
      while( this->full() || q2.empty() )
      {
        m_notFullCond.wait( lock, [ this ]  { return !this->full(); } );
        q2.m_notEmptyCond.wait( lock, [ &q2 ] { return !q2.empty(); } );
      }
      this->emplace_back( q2.front() ).wait();
      q2.pop_front();
    }
    m_notEmptyCond.notify_all();
    q2.m_notFullCond.notify_all();
  }
};

/**
 * This class is used to store in a LIFO way buffers, first on device, then on host, then on disk.
 */
template< typename T >
class LifoStorage
{

private:

  /// number of buffers to be inserted into the LIFO
  int m_maxNumberOfBuffers;
  /// size of one buffer in bytes
  size_t m_bufferSize;
  /// name used to store data on disk
  std::string m_name;
#ifdef GEOSX_USE_CUDA
  /// ueue of data stored on device
  std::unique_ptr< FixedSizeDequeAndMutexes< T > > m_deviceDeque;
#endif
  /// ueue of data stored on host memory
  std::unique_ptr< FixedSizeDequeAndMutexes< T > > m_hostDeque;

  /// counter of buffer stored in LIFO
  int m_bufferCount;
  /// counter of buffer pushed to host
  int m_bufferToHostCount;
  /// counter of buffer pushed to disk
  int m_bufferToDiskCount;

#ifdef GEOSX_USE_CUDA
  // Events associated to ith  copies to device buffer
  std::vector< camp::resources::Event > m_pushToDeviceEvents;
#endif
  // Futures associated to push to host in case we have no device buffers
  std::vector< std::future< void > > m_pushToHostFutures;
#ifdef GEOSX_USE_CUDA
  // Events associated to ith  copies from device buffer
  std::vector< camp::resources::Event > m_popFromDeviceEvents;
#endif
  // Futures associated to pop from host in case we have no device buffers
  std::vector< std::future< void > > m_popFromHostFutures;

  /// condition used to tell m_worker queue has been filled or processed is stopped.
  std::condition_variable m_task_queue_not_empty_cond[2];
  /// mutex to protect access to m_task_queue.
  std::mutex m_task_queue_mutex[2];
  /// queue of task to be executed by m_worker.
  std::deque< std::packaged_task< void() > > m_task_queue[2];
  /// thread to execute tasks.
  std::thread m_worker[2];
  /// boolean to keep m_worker alive.
  bool m_continue;
  /// marker to detect first pop
  bool m_hasPoppedBefore = false;
public:


  /**
   * A LIFO storage will store numberOfBuffersToStoreDevice buffer on
   * deevice, numberOfBuffersToStoreHost on host and the rest on disk.
   *
   * @param name                           Prefix of the files used to save the occurenncy of the saved buffer on disk.
   * @param elemCnt                        Number of elments in the LvArray we want to store in the LIFO storage.
   * @param numberOfBuffersToStoreOnDevice Maximum number of array to store on device memory ( -1 = use 80% of remaining memory ).
   * @param numberOfBuffersToStoreOnHost   Maximum number of array to store on host memory ( -1 = use 80% of remaining memory ).
   * @param maxNumberOfBuffers             Number of arrays expected to be stores in the LIFO.
   */
  LifoStorage( std::string name, size_t elemCnt, int numberOfBuffersToStoreOnDevice, int numberOfBuffersToStoreOnHost, int maxNumberOfBuffers ):
    m_maxNumberOfBuffers( maxNumberOfBuffers ),
    m_bufferSize( elemCnt*sizeof( T ) ),
    m_name( name ),
    m_bufferCount( 0 ), m_bufferToHostCount( 0 ), m_bufferToDiskCount( 0 ),
#ifndef GEOSX_USE_CUDA
    m_pushToHostFutures( maxNumberOfBuffers ),
    m_popFromHostFutures( maxNumberOfBuffers ),
#endif
    m_continue( true )
  {
    LIFO_LOG_RANK(" LIFO : maximum size "<< m_maxNumberOfBuffers << " buffers ");
    double bufferSize = ( ( double ) m_bufferSize ) / ( 1024.0 * 1024.0 );
    LIFO_LOG_RANK(" LIFO : buffer size "<< bufferSize << "MB");
#ifndef GEOSX_USE_CUDA
    numberOfBuffersToStoreOnDevice = 0;
#else
    if ( numberOfBuffersToStoreOnDevice == -1 )
    {
      size_t free, total;
      GEOSX_ERROR_IF( cudaSuccess != cudaMemGetInfo( &free, &total ), "Error getting CUDA device available memory" );
      double freeGB = ( ( double ) free ) / ( 1024.0 * 1024.0 * 1024.0 );
      LIFO_LOG_RANK(" LIFO : available memory on device " << freeGB << " GB");
      numberOfBuffersToStoreOnDevice = std::min( (int)( 0.8 * free / m_bufferSize ), m_maxNumberOfBuffers );
    }
    m_deviceDeque = std::unique_ptr< FixedSizeDequeAndMutexes< T > >( new FixedSizeDequeAndMutexes< T >( numberOfBuffersToStoreOnDevice, elemCnt, LvArray::MemorySpace::cuda ) );
    m_pushToDeviceEvents = std::vector< camp::resources::Event >( (numberOfBuffersToStoreOnDevice > 0)?maxNumberOfBuffers:0 );
    m_pushToHostFutures = std::vector< std::future< void > >( (numberOfBuffersToStoreOnDevice > 0)?0:maxNumberOfBuffers );
    m_popFromDeviceEvents = std::vector< camp::resources::Event >( (numberOfBuffersToStoreOnDevice > 0)?maxNumberOfBuffers:0 );
    m_popFromHostFutures = std::vector< std::future< void > >( (numberOfBuffersToStoreOnDevice > 0)?0:maxNumberOfBuffers );
    LIFO_LOG_RANK(" LIFO : allocating "<< numberOfBuffersToStoreOnDevice <<" buffers on device");
#endif
    if ( numberOfBuffersToStoreOnHost == -1 )
    {
      size_t free = sysconf(_SC_AVPHYS_PAGES) * sysconf(_SC_PAGESIZE);
      numberOfBuffersToStoreOnHost = std::max( 1 , std::min( (int)( 0.8 * free / m_bufferSize ), m_maxNumberOfBuffers - numberOfBuffersToStoreOnDevice ) );
      double freeGB = ( ( double ) free ) / ( 1024.0 * 1024.0 * 1024.0 );
      LIFO_LOG_RANK(" LIFO : available memory on host " << freeGB << " GB");
    }
    LIFO_LOG_RANK(" LIFO : allocating "<< numberOfBuffersToStoreOnHost <<" buffers on host");
    m_hostDeque = std::unique_ptr< FixedSizeDequeAndMutexes< T > >( new FixedSizeDequeAndMutexes< T >( numberOfBuffersToStoreOnHost, elemCnt, LvArray::MemorySpace::host ) );
    m_worker[0] = std::thread( &LifoStorage< T >::wait_and_consume_tasks, this, 0 );
    m_worker[1] = std::thread( &LifoStorage< T >::wait_and_consume_tasks, this, 1 );
  }

  /**
   * Build a LIFO storage for a given LvArray array.
   *
   * @param name                           Prefix of the files used to save the occurenncy of the saved buffer on disk.
   * @param array                          The LvArray that will be store in the LIFO.
   * @param numberOfBuffersToStoreOnDevice Maximum number of array to store on device memory.
   * @param numberOfBuffersToStoreOnHost   Maximum number of array to store on host memory.
   * @param maxNumberOfBuffers             Number of arrays expected to be stores in the LIFO.
   */
  LifoStorage( std::string name, arrayView1d< T > array, int numberOfBuffersToStoreOnDevice, int numberOfBuffersToStoreOnHost, int maxNumberOfBuffers ):
    LifoStorage( name, array.size(), numberOfBuffersToStoreOnDevice, numberOfBuffersToStoreOnHost, maxNumberOfBuffers ) {}

  ~LifoStorage()
  {
    m_continue = false;
    m_task_queue_not_empty_cond[0].notify_all();
    m_task_queue_not_empty_cond[1].notify_all();
    m_worker[0].join();
    m_worker[1].join();
  }

  /**
   * Asynchroneously push a copy of the given LvArray into the LIFO
   *
   * @param array The LvArray to store in the LIFO, should match the size of the data used in constructor.
   */
  void pushAsync( arrayView1d< T > array )
  {
    LIFO_MARK_FUNCTION;
    //To be sure 2 pushes are not mixed
    pushWait();
    int id = m_bufferCount++;
    GEOSX_ERROR_IF( m_hostDeque->capacity() == 0,
                    "Cannot save on a Lifo without host storage (please set lifoSize, lifoOnDevice and lifoOnHost in xml file)" );

#ifdef GEOSX_USE_CUDA
    if( m_deviceDeque->capacity() > 0 )
    {
      m_pushToDeviceEvents[id] = m_deviceDeque->emplaceFront( array );

      if( m_maxNumberOfBuffers - id > (int)m_deviceDeque->capacity() )
      {
        LIFO_MARK_SCOPE( geosx::lifoStorage< T >::pushAddTasks );
        // This buffer will go to host memory, and maybe on disk
        std::packaged_task< void() > task( std::bind( &LifoStorage< T >::deviceToHost, this, m_bufferToHostCount++ ) );
        {
          std::unique_lock< std::mutex > lock( m_task_queue_mutex[0] );
          m_task_queue[0].emplace_back( std::move( task ) );
        }
        m_task_queue_not_empty_cond[0].notify_all();
      }
    }
    else
#endif
    {
      std::packaged_task< void() > task( std::bind( [ this ] ( int pushId, arrayView1d< T > pushedArray ) {
        m_hostDeque->emplaceFront( pushedArray );

        if( m_maxNumberOfBuffers - pushId > (int)m_hostDeque->capacity() )
        {
          LIFO_MARK_SCOPE( geosx::lifoStorage< T >::pushAddTasks );
          // This buffer will go to host memory, and maybe on disk
          std::packaged_task< void() > t2( std::bind( &LifoStorage< T >::hostToDisk, this, m_bufferToDiskCount++ ) );
          {
            std::unique_lock< std::mutex > l2( m_task_queue_mutex[1] );
            m_task_queue[1].emplace_back( std::move( t2 ) );
          }
          m_task_queue_not_empty_cond[1].notify_all();
        }
      }, id, array ) );
      m_pushToHostFutures[id] = task.get_future();
      {
        std::unique_lock< std::mutex > lock( m_task_queue_mutex[0] );
        m_task_queue[0].emplace_back( std::move( task ) );
      }
      m_task_queue_not_empty_cond[0].notify_all();
    }

  }

  /**
   * Waits for last push to be terminated
   */
  void pushWait()
  {
    LIFO_MARK_FUNCTION;

    if( m_bufferCount > 0 )
    {
#ifdef GEOSX_USE_CUDA
      if( m_deviceDeque->capacity() > 0 )
      {
        m_pushToDeviceEvents[m_bufferCount-1].wait();
      }
      else
#endif
      {
        m_pushToHostFutures[m_bufferCount-1].wait();
      }
    }
  }

  /**
   * Push a copy of the given LvArray into the LIFO
   *
   * @param array The LvArray to store in the LIFO, should match the size of the data used in constructor.
   */
  void push( arrayView1d< T > array )
  {
    LIFO_MARK_FUNCTION;
    pushAsync( array );
    pushWait();
  }

  /**
   * Asynchroneously copy last data from the LIFO into the LvArray.
   *
   * @param array LvArray to store data from the LIFO into it.
   */
  void popAsync( arrayView1d< T > array )
  {
    LIFO_MARK_FUNCTION;
    //wait the last push to avoid race condition
    pushWait();
    if ( m_hasPoppedBefore )
    {
      // Ensure last pop is finished
      popWait();
    }
    else
    {
      if ( m_maxNumberOfBuffers != m_bufferCount ) LIFO_LOG_RANK(" LIFO : warning number of entered buffered (" << m_bufferCount
                                                                 << ") != max LIFO size (" << m_maxNumberOfBuffers << ") !" );
      // Ensure that all push step are ended
      for ( int queueId = 0; queueId < 2; queueId++ )
      {
        std::unique_lock< std::mutex > lock( m_task_queue_mutex[queueId] );
        m_task_queue_not_empty_cond[queueId].wait( lock, [ this, &queueId ] { return m_task_queue[queueId].empty(); } );
      }
    }
    m_hasPoppedBefore = true;
    int id = --m_bufferCount;

#ifdef GEOSX_USE_CUDA
    if( m_deviceDeque->capacity() > 0 )
    {

      if( m_bufferToHostCount > 0 )
      {
        LIFO_MARK_SCOPE( geosx::LifoStorage< T >::popAddTasks );
        // Trigger pull one buffer from host, and maybe from disk
        std::packaged_task< void() > task( std::bind( &LifoStorage< T >::hostToDevice, this, --m_bufferToHostCount, id ) );
        {
          std::unique_lock< std::mutex > lock( m_task_queue_mutex[0] );
          m_task_queue[0].emplace_back( std::move( task ) );
        }
        m_task_queue_not_empty_cond[0].notify_all();
      }

      m_popFromDeviceEvents[id] = m_deviceDeque->popFront( array );
    }
    else
#endif
    {
      std::packaged_task< void() > task( std::bind ( [ this ] ( int popId, arrayView1d< T > poppedArray ) {
        m_hostDeque->popFront( poppedArray );

        if( m_bufferToDiskCount > 0 )
        {
          LIFO_MARK_SCOPE( geosx::LifoStorage< T >::popAddTasks );
          // Trigger pull one buffer from host, and maybe from disk
          std::packaged_task< void() > task2( std::bind( &LifoStorage< T >::diskToHost, this, --m_bufferToDiskCount ) );
          {
            std::unique_lock< std::mutex > lock2( m_task_queue_mutex[1] );
            m_task_queue[1].emplace_back( std::move( task2 ) );
          }
          m_task_queue_not_empty_cond[1].notify_all();
        }
      }, id, array ) );
      m_popFromHostFutures[id] = task.get_future();
      {
        std::unique_lock< std::mutex > lock( m_task_queue_mutex[0] );
        m_task_queue[0].emplace_back( std::move( task ) );
      }
      m_task_queue_not_empty_cond[0].notify_all();
    }
  }

  /**
   * Waits for last pop to be terminated
   */
  void popWait()
  {
    LIFO_MARK_FUNCTION;
    if( m_bufferCount < m_maxNumberOfBuffers )
    {
#ifdef GEOSX_USE_CUDA
      if( m_deviceDeque->capacity() > 0 )
      {
        auto *cuda_event = m_popFromDeviceEvents[m_bufferCount].try_get< camp::resources::CudaEvent >();
        if ( cuda_event ) cudaEventSynchronize( cuda_event->getCudaEvent_t() );
      }
      else
#endif
      {
        m_popFromHostFutures[m_bufferCount].wait();
      }
    }
  }

  /**
   * Copy last data from the LIFO into the LvArray.
   *
   * @param array LvArray to store data from the LIFO into it.
   */
  void pop( arrayView1d< T > array )
  {
    LIFO_MARK_FUNCTION;
    popAsync( array );
    popWait();
  }

  /**
   * @brief Returns true if the LIFO does not contain a buffer.
   */
  bool isEmpty()
  {
    return m_bufferCount == 0;
  }
private:

  /**
   * Copy data from device memory to host memory
   *
   * @param ID of the buffer to copy on host.
   */
#ifdef GEOSX_USE_CUDA
  void deviceToHost( int id )
  {
    LIFO_MARK_FUNCTION;
    // The copy to host will only start when the data is copied on device buffer
    m_hostDeque->getStream().wait_for( const_cast< camp::resources::Event * >( &m_pushToDeviceEvents[id] ) );
    m_hostDeque->emplaceFrontFromBack( *m_deviceDeque );

    if( m_maxNumberOfBuffers - id > (int)(m_deviceDeque->capacity() + m_hostDeque->capacity()) )
    {
      // This buffer will go to host then maybe to disk
      std::packaged_task< void() > task( std::bind( &LifoStorage< T >::hostToDisk, this, m_bufferToDiskCount++ ) );
      {
        std::unique_lock< std::mutex > lock( m_task_queue_mutex[1] );
        m_task_queue[1].emplace_back( std::move( task ) );
      }
      m_task_queue_not_empty_cond[1].notify_all();
    }
  }
#endif

  /**
   * Copy data from host memory to disk
   *
   * @param id ID of the buffer to store on disk.
   */
  void hostToDisk( int id )
  {
    LIFO_MARK_FUNCTION;
    {
      TwoMutexLock lock( m_hostDeque->m_popMutex, m_hostDeque->m_backMutex );
      writeOnDisk( m_hostDeque->back().dataIfContiguous(), id );
      m_hostDeque->pop_back();
    }
    m_hostDeque->m_notFullCond.notify_all();
  }

  /**
   * Copy data from host memory to device memory
   *
   * @param id ID of the buffer to load from host memory.
   * @param id_pop ID of the last popped buffer from device
   */
#ifdef GEOSX_USE_CUDA
  void hostToDevice( int id, int id_pop )
  {
    LIFO_MARK_FUNCTION;
    // enqueue diskToHost on worker #2 if needed
    if( m_bufferToDiskCount > 0 )
    {
      // This buffer will go to host then to disk
      std::packaged_task< void() > task( std::bind( &LifoStorage< T >::diskToHost, this, --m_bufferToDiskCount ) );
      {
        std::unique_lock< std::mutex > lock( m_task_queue_mutex[1] );
        m_task_queue[1].emplace_back( std::move( task ) );
      }
      m_task_queue_not_empty_cond[1].notify_all();
    }

    m_deviceDeque->emplaceBackFromFront( *m_hostDeque );
  }
#endif

  /**
   * Copy data from disk to host memory
   *
   * @param id ID of the buffer to read on disk.
   */
  void diskToHost( int id )
  {
    LIFO_MARK_FUNCTION;
    {
      TwoMutexLock lock( m_hostDeque->m_emplaceMutex, m_hostDeque->m_backMutex );
      m_hostDeque->m_notFullCond.wait( lock, [ this ]  { return !( m_hostDeque->full() ); } );
      readOnDisk( const_cast< T * >(m_hostDeque->next_back().dataIfContiguous()), id );
      m_hostDeque->inc_back();
    }
    m_hostDeque->m_notEmptyCond.notify_all();
  }
  /**
   * Checks if a directory exists.
   *
   * @param dirName Directory name to check existence of.
   * @return true is dirName exists and is a directory.
   */
  bool dirExists( const std::string & dirName )
  {
    struct stat buffer;
    return stat( dirName.c_str(), &buffer ) == 0;
  }

  /**
   * Write data on disk
   *
   * @param d Data to store on disk.
   * @param id ID of the buffer to read on disk
   */
  void writeOnDisk( const T * d, int id )
  {
    LIFO_MARK_FUNCTION;
    std::string fileName = GEOSX_FMT( "{}_{:08}.dat", m_name, id );
    int lastDirSeparator = fileName.find_last_of( "/\\" );
    std::string dirName = fileName.substr( 0, lastDirSeparator );
    if( string::npos != (size_t)lastDirSeparator && !dirExists( dirName ))
      makeDirsForPath( dirName );

    std::ofstream wf( fileName, std::ios::out | std::ios::binary );
    GEOSX_ERROR_IF( !wf || wf.fail() || !wf.is_open(),
                    "Could not open file "<< fileName << " for writting" );
    wf.write( (const char *)d, m_bufferSize );
    GEOSX_ERROR_IF( wf.bad() || wf.fail(),
                    "An error occured while writting "<< fileName );
    wf.close();
  }

  /**
   * Read data from disk
   *
   * @param d  Buffer to store data read from disk.
   * @param id ID of the buffer on disk.
   */
  void readOnDisk( T * d, int id )
  {
    LIFO_MARK_FUNCTION;
    std::string fileName = GEOSX_FMT( "{}_{:08}.dat", m_name, id );
    std::ifstream wf( fileName, std::ios::in | std::ios::binary );
    GEOSX_ERROR_IF( !wf,
                    "Could not open file "<< fileName << " for reading" );
    wf.read( (char *)d, m_bufferSize );
    wf.close();
    remove( fileName.c_str() );
  }

  /**
   * Execute the tasks pushed in m_task_queue in the given order.
   *
   * @param queueId index of the queue (0: queue to handle device/host transfers; 1: host/disk transfers)
   */
  void wait_and_consume_tasks( int queueId )
  {
    LIFO_MARK_FUNCTION;
    while( m_continue )
    {
      std::unique_lock< std::mutex > lock( m_task_queue_mutex[queueId] );
      {
        LIFO_MARK_SCOPE( waitForTask );
        m_task_queue_not_empty_cond[queueId].wait( lock, [ this, &queueId ] { return !( m_task_queue[queueId].empty()  && m_continue ); } );
      }
      if( m_continue == false ) break;
      std::packaged_task< void() > task( std::move( m_task_queue[queueId].front() ) );
      m_task_queue[queueId].pop_front();
      lock.unlock();
      m_task_queue_not_empty_cond[queueId].notify_all();
      {
        LIFO_MARK_SCOPE( runningTask );
        task();
      }
    }
  }
};
}
#endif // LIFOSTORAGE_HPP
