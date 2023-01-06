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

#include "common/fixedSizeDeque.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "common/TimingMacros.hpp"

#define LIFO_MARK_FUNCTION if( std::getenv( "LIFO_TRACE_ON" ) != NULL ) GEOSX_MARK_FUNCTION;
#define LIFO_MARK_SCOPE( a )  if( std::getenv( "LIFO_TRACE_ON" ) != NULL ) GEOSX_MARK_SCOPE( a );


namespace geosx
{


class twoMutexLock
{
public:
  twoMutexLock( std::mutex & mutex1, std::mutex & mutex2 ):
    m_islocked( false ), m_mutex1( &mutex1 ), m_mutex2( &mutex2 )
  {
    lock();
  }

  ~twoMutexLock()
  {
    unlock();
  }

  void lock()
  {
    if( m_islocked ) return;
    std::lock( *m_mutex1, *m_mutex2 );
    m_islocked = true;
  }

  void unlock()
  {
    if( !m_islocked ) return;
    m_mutex1->unlock();
    m_mutex2->unlock();
    m_islocked = false;
  }
private:
  bool m_islocked;
  std::mutex *m_mutex1;
  std::mutex *m_mutex2;
};

class fourMutexLock
{
public:
  fourMutexLock( std::mutex & mutex1, std::mutex & mutex2, std::mutex & mutex3, std::mutex & mutex4 ):
    m_islocked( false ), m_mutex1( &mutex1 ), m_mutex2( &mutex2 ), m_mutex3( &mutex3 ), m_mutex4( &mutex4 )
  {
    lock();
  }

  ~fourMutexLock()
  {
    unlock();
  }

  void lock()
  {
    if( m_islocked ) return;
    std::lock( *m_mutex1, *m_mutex2, *m_mutex3, *m_mutex4 );
    m_islocked = true;
  }

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
  bool m_islocked;
  std::mutex *m_mutex1;
  std::mutex *m_mutex2;
  std::mutex *m_mutex3;
  std::mutex *m_mutex4;
};

/// Associate mutexes with the fixedSizeDeque
template< typename T >
class fixedSizeDequeAndMutexes : public fixedSizeDeque< T, int >
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
  fixedSizeDequeAndMutexes( int maxEntries, int valuesPerEntry, LvArray::MemorySpace space ): fixedSizeDeque< T, int >( maxEntries, valuesPerEntry, space,
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
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    camp::resources::Event e;
    {
      twoMutexLock lock( m_emplaceMutex, m_frontMutex );
      {
        LIFO_MARK_SCOPE( waitingForBuffer );
        m_notFullCond.wait( lock, [ this ]  { return !this->full(); } );
      }
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      {
        LIFO_MARK_SCOPE( copy );
        e = fixedSizeDeque< T, int >::emplace_front( array.toSliceConst() );
      }
    }
    m_notEmptyCond.notify_all();
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
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
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    LIFO_MARK_FUNCTION;
    camp::resources::Event e;
    {
      twoMutexLock lock( m_popMutex, m_frontMutex );
      {
        LIFO_MARK_SCOPE( waitingForBuffer );
        m_notEmptyCond.wait( lock, [ this ]  { return !this->empty(); } );
      }
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
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
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    return e;
  }

  /**
   * Emplace front from back of given queue
   *
   * @param q2 The queue to copy data from.
   */
  void emplaceFrontFromBack( fixedSizeDequeAndMutexes< T > & q2 )
  {
    LIFO_MARK_FUNCTION;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    {
      fourMutexLock lock( m_emplaceMutex, q2.m_popMutex, m_frontMutex, q2.m_backMutex );
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
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
  }

  /**
   * Emplace back from front of given queue
   *
   * @param q2 The queue to copy data from.
   */
  void emplaceBackFromFront( fixedSizeDequeAndMutexes< T > & q2 )
  {
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    LIFO_MARK_FUNCTION;
    {
      fourMutexLock lock( m_emplaceMutex, q2.m_popMutex, m_backMutex, q2.m_frontMutex );
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
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
  }
};

/**
 * This class is used to store in a LIFO way buffers, first on device, then on host, then on disk.
 */
template< typename T >
class lifoStorage
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
  fixedSizeDequeAndMutexes< T > m_deviceDeque;
#endif
  /// ueue of data stored on host memory
  fixedSizeDequeAndMutexes< T > m_hostDeque;

  /// counter of buffer stored in LIFO
  int m_bufferCount;
  /// counter of buffer stored on disk
  int m_bufferOnDiskCount;

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
  bool m_continue = true;

public:


  /**
   * A LIFO storage will store numberOfBuffersToStoreDevice buffer on
   * deevice, numberOfBuffersToStoreHost on host and the rest on disk.
   *
   * @param name                           Prefix of the files used to save the occurenncy of the saved buffer on disk.
   * @param elemCnt                        Number of elments in the LvArray we want to store in the LIFO storage.
   * @param numberOfBuffersToStoreOnDevice Maximum number of array to store on device memory.
   * @param numberOfBuffersToStoreOnHost   Maximum number of array to store on host memory.
   * @param maxNumberOfBuffers             Number of arrays expected to be stores in the LIFO.
   */
  lifoStorage( std::string name, size_t elemCnt, int numberOfBuffersToStoreOnDevice, int numberOfBuffersToStoreOnHost, int maxNumberOfBuffers ):
    m_maxNumberOfBuffers( maxNumberOfBuffers ),
    m_bufferSize( elemCnt*sizeof( T ) ),
    m_name( name ),
#ifdef GEOSX_USE_CUDA
    m_deviceDeque( numberOfBuffersToStoreOnDevice, elemCnt, LvArray::MemorySpace::cuda ),
#endif
    m_hostDeque( numberOfBuffersToStoreOnHost, elemCnt, LvArray::MemorySpace::host ),
    m_bufferCount( 0 ), m_bufferOnDiskCount( 0 ),
#ifdef GEOSX_USE_CUDA
    m_pushToDeviceEvents( (numberOfBuffersToStoreOnDevice > 0)?maxNumberOfBuffers:0 ),
    m_pushToHostFutures( (numberOfBuffersToStoreOnDevice > 0)?0:maxNumberOfBuffers ),
    m_popFromDeviceEvents( (numberOfBuffersToStoreOnDevice > 0)?maxNumberOfBuffers:0 ),
    m_popFromHostFutures( (numberOfBuffersToStoreOnDevice > 0)?0:maxNumberOfBuffers )
#else
    m_pushToHostFutures( maxNumberOfBuffers ),
    m_popFromHostFutures( maxNumberOfBuffers )
#endif
  {
#ifndef GEOSX_USE_CUDA
    GEOSX_UNUSED_VAR( numberOfBuffersToStoreOnDevice );
#endif
    m_worker[0] = std::thread( &lifoStorage< T >::wait_and_consume_tasks, this, 0 );
    m_worker[1] = std::thread( &lifoStorage< T >::wait_and_consume_tasks, this, 1 );
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
  lifoStorage( std::string name, arrayView1d< T > array, int numberOfBuffersToStoreOnDevice, int numberOfBuffersToStoreOnHost, int maxNumberOfBuffers ):
    lifoStorage( name, array.size(), numberOfBuffersToStoreOnDevice, numberOfBuffersToStoreOnHost, maxNumberOfBuffers ) {}

  ~lifoStorage()
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
    GEOSX_ERROR_IF( m_hostDeque.capacity() == 0,
                    "Cannot save on a Lifo without host storage (please set lifoSize, lifoOnDevice and lifoOnHost in xml file)" );

#ifdef GEOSX_USE_CUDA
    if( m_deviceDeque.capacity() > 0 )
    {
      m_pushToDeviceEvents[id] = m_deviceDeque.emplaceFront( array );

      if( m_maxNumberOfBuffers - id > (int)m_deviceDeque.capacity() )
      {
        LIFO_MARK_SCOPE( geosx::lifoStorage< T >::pushAddTasks );
        // This buffer will go to host memory, and maybe on disk
        std::packaged_task< void() > task( std::bind( &lifoStorage< T >::deviceToHost, this, id ) );
        {
          std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
          std::unique_lock< std::mutex > lock( m_task_queue_mutex[0] );
          std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
          m_task_queue[0].emplace_back( std::move( task ) );
        }
        m_task_queue_not_empty_cond[0].notify_all();
      }
    }
    else
#endif
    {
      std::packaged_task< void() > task( std::bind( [ this ] ( int pushId, arrayView1d< T > pushedArray ) {
        m_hostDeque.emplaceFront( pushedArray );

        if( m_maxNumberOfBuffers - pushId > (int)m_hostDeque.capacity() )
        {
          LIFO_MARK_SCOPE( geosx::lifoStorage< T >::pushAddTasks );
          // This buffer will go to host memory, and maybe on disk
          std::packaged_task< void() > t2( std::bind( &lifoStorage< T >::hostToDisk, this, pushId ) );
          {
            std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
            std::unique_lock< std::mutex > l2( m_task_queue_mutex[1] );
            std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
            m_task_queue[1].emplace_back( std::move( t2 ) );
          }
          m_task_queue_not_empty_cond[1].notify_all();
        }
      }, id, array ) );
      m_pushToHostFutures[id] = task.get_future();
      {
        std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
        std::unique_lock< std::mutex > lock( m_task_queue_mutex[0] );
        std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
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
      if( m_deviceDeque.capacity() > 0 )
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
    // Ensure last pop is finished
    popWait();
    int id = --m_bufferCount;
#ifdef GEOSX_USE_CUDA
    if( m_deviceDeque.capacity() > 0 )
    {
      m_popFromDeviceEvents[id] = m_deviceDeque.popFront( array );

      if( id >= (int)m_deviceDeque.capacity() )
      {
        LIFO_MARK_SCOPE( geosx::lifoStorage< T >::popAddTasks );
        // Trigger pull one buffer from host, and maybe from disk
        std::packaged_task< void() > task( std::bind( &lifoStorage< T >::hostToDevice, this, id - m_deviceDeque.capacity() ) );
        {
          std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
          std::unique_lock< std::mutex > lock( m_task_queue_mutex[0] );
          std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
          m_task_queue[0].emplace_back( std::move( task ) );
        }
        m_task_queue_not_empty_cond[0].notify_all();
      }
    }
    else
#endif
    {
      std::packaged_task< void() > task( std::bind ( [ this ] ( int popId, arrayView1d< T > poppedArray ) {
        m_hostDeque.popFront( poppedArray );

        if( popId >= (int)m_hostDeque.capacity() )
        {
          LIFO_MARK_SCOPE( geosx::lifoStorage< T >::popAddTasks );
          // Trigger pull one buffer from host, and maybe from disk
          std::packaged_task< void() > task2( std::bind( &lifoStorage< T >::diskToHost, this, popId  - m_hostDeque.capacity() ) );
          {
            std::cerr << __FILE__ << ":" << __LINE__ << " " << popId << std::endl;
            std::unique_lock< std::mutex > lock2( m_task_queue_mutex[1] );
            std::cerr << __FILE__ << ":" << __LINE__ << " " << popId << std::endl;
            m_task_queue[1].emplace_back( std::move( task2 ) );
          }
          m_task_queue_not_empty_cond[1].notify_all();
        }
      }, id, array ) );
      m_popFromHostFutures[id] = task.get_future();
      {
        std::cerr << __FILE__ << ":" << __LINE__ << " " << id << std::endl;
        std::unique_lock< std::mutex > lock( m_task_queue_mutex[0] );
        std::cerr << __FILE__ << ":" << __LINE__ << " " << id << std::endl;
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
      if( m_deviceDeque.capacity() > 0 )
      {
#ifdef GEOSX_USE_CUDA
        cudaEventSynchronize( m_popFromDeviceEvents[m_bufferCount].get< camp::resources::CudaEvent >().getCudaEvent_t() );
#endif
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
    m_hostDeque.getStream().wait_for( const_cast< camp::resources::Event * >( &m_pushToDeviceEvents[id] ) );
    m_hostDeque.emplaceFrontFromBack( m_deviceDeque );

    if( m_maxNumberOfBuffers - id > (int)(m_deviceDeque.capacity() + m_hostDeque.capacity()) )
    {
      // This buffer will go to host then maybe to disk
      std::packaged_task< void() > task( std::bind( &lifoStorage< T >::hostToDisk, this, id ) );
      {
        std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
        std::unique_lock< std::mutex > lock( m_task_queue_mutex[1] );
        std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
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
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      twoMutexLock lock( m_hostDeque.m_popMutex, m_hostDeque.m_backMutex );
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      writeOnDisk( m_hostDeque.back().dataIfContiguous(), id );
      m_hostDeque.pop_back();
    }
    m_hostDeque.m_notFullCond.notify_all();
  }

  /**
   * Copy data from host memory to device memory
   *
   * @param id ID of the buffer to load from host memory.
   */
#ifdef GEOSX_USE_CUDA
  void hostToDevice( int id )
  {
    LIFO_MARK_FUNCTION;
    m_hostDeque.getStream().wait_for( const_cast< camp::resources::Event * >( &m_popFromDeviceEvents[ id + m_deviceDeque.capacity() ] ) );
    m_deviceDeque.emplaceBackFromFront( m_hostDeque );

    // enqueue diskToHost on worker #2 if needed
    if( id >= (int)m_hostDeque.capacity() )
    {
      // This buffer will go to host then to disk
      std::packaged_task< void() > task( std::bind( &lifoStorage< T >::diskToHost, this, id - m_hostDeque.capacity() ) );
      {
        std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
        std::unique_lock< std::mutex > lock( m_task_queue_mutex[1] );
        std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
        m_task_queue[1].emplace_back( std::move( task ) );
      }
      m_task_queue_not_empty_cond[1].notify_all();
    }
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
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      twoMutexLock lock( m_hostDeque.m_emplaceMutex, m_hostDeque.m_backMutex );
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      m_hostDeque.m_notFullCond.wait( lock, [ this ]  { return !( m_hostDeque.full() ); } );
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      readOnDisk( const_cast< T * >(m_hostDeque.next_back().dataIfContiguous()), id );
      m_hostDeque.inc_back();
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    }
    m_hostDeque.m_notEmptyCond.notify_all();
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
      std::cerr << __FILE__ << ":" << __LINE__ << " " << queueId << std::endl;
      std::unique_lock< std::mutex > lock( m_task_queue_mutex[queueId] );
      std::cerr << __FILE__ << ":" << __LINE__ << " " << queueId << std::endl;
      {
        LIFO_MARK_SCOPE( waitForTask );
        m_task_queue_not_empty_cond[queueId].wait( lock, [ this, &queueId ] { return !( m_task_queue[queueId].empty()  && m_continue ); } );
        std::cerr << __FILE__ << ":" << __LINE__ << " " << queueId << std::endl;
      }
      if( m_continue == false ) break;
      std::packaged_task< void() > task( std::move( m_task_queue[queueId].front() ) );
      m_task_queue[queueId].pop_front();
      lock.unlock();
      std::cerr << __FILE__ << ":" << __LINE__ << " " << queueId << std::endl;
      {
        LIFO_MARK_SCOPE( runningTask );
        task();
      }
      std::cerr << __FILE__ << ":" << __LINE__ << " " << queueId << std::endl;
    }
    std::cerr << __FILE__ << ":" << __LINE__ << " " << queueId << std::endl;
  }
};
}
#endif // LIFOSTORAGE_HPP
