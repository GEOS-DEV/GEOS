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
#include <condition_variable>
#include <camp/camp.hpp>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "common/fixedSizeDeque.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "common/TimingMacros.hpp"

#define LIFO_MARK_FUNCTION if( std::getenv( "LIFO_TRACE_ON" ) != NULL ) GEOSX_MARK_FUNCTION;
#define LIFO_MARK_SCOPE( a )  if( std::getenv( "LIFO_TRACE_ON" ) != NULL ) GEOSX_MARK_SCOPE( a );
namespace geosx
{

/// Associate mutexes with the fixedSizeDeque
template< typename T >
class fixedSizeDequeAndMutexes : public fixedSizeDeque< T, int >
{
public:
  // Mutex to protect access to the front
  std::mutex m_frontMutex;
  // Mutex to protect access to the back
  std::mutex m_backMutex;
  // Mutex to prevent two simulteaneous pop (can be an issue for last one)
  std::mutex m_popMutex;
  // Mutex to prevent two simulteaneous emplace (can be an issue for last one)
  std::mutex m_emplaceMutex;
  /// Condition used to notify when device queue is not full
  std::condition_variable m_notFullCond;
  /// Condition used to notify when device queue is not empty
  std::condition_variable m_notEmptyCond;

  /**
   * Create a fixed size double ended queue with associated mutexes and condition variables.
   *
   * @param maxEntries     Maximum number of array to store in the queue.
   * @param valuesPerEntry Number of values in each array of the deque.
   * @param space          Space used to store que queue.
   */
  fixedSizeDequeAndMutexes( int maxEntries, int valuesPerEntry, LvArray::MemorySpace space ): fixedSizeDeque< T, int >( maxEntries, valuesPerEntry, space,
                                                                                                                        camp::resources::Resource{ camp::resources::Cuda{} } ) {}

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
      {
        LIFO_MARK_SCOPE( waitingForMutex );
        std::lock( m_emplaceMutex, m_frontMutex );
      }
      std::unique_lock< std::mutex > l1( m_emplaceMutex, std::adopt_lock );
      std::unique_lock< std::mutex > l2( m_frontMutex, std::adopt_lock );
      {
        LIFO_MARK_SCOPE( waitingForBuffer );
        m_notFullCond.wait( l1, [ this ]  { return !this->full(); } );
      }
      {
        LIFO_MARK_SCOPE( copy );
        e = fixedSizeDeque< T, int >::emplace_front( array.toSliceConst() );
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
      {
        LIFO_MARK_SCOPE( waitingForMutex );
        std::lock( m_popMutex, m_frontMutex );
      }
      std::unique_lock< std::mutex > l1( m_popMutex, std::adopt_lock );
      std::unique_lock< std::mutex > l2( m_frontMutex, std::adopt_lock );
      {
        LIFO_MARK_SCOPE( waitingForBuffer );
        m_notEmptyCond.wait( l1, [ this ]  { return !this->empty(); } );
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
  void emplaceFrontFromBack( fixedSizeDequeAndMutexes< T > & q2 )
  {
    LIFO_MARK_FUNCTION;
    {
      std::lock( m_emplaceMutex, m_frontMutex, q2.m_popMutex, q2.m_backMutex );
      std::unique_lock< std::mutex > lockQ1Emplace( m_emplaceMutex, std::adopt_lock );
      std::unique_lock< std::mutex > lockQ1Front( m_frontMutex, std::adopt_lock );
      std::unique_lock< std::mutex > lockQ2Pop( q2.m_popMutex, std::adopt_lock );
      std::unique_lock< std::mutex > lockQ2Back( q2.m_backMutex, std::adopt_lock );
      {
        LIFO_MARK_SCOPE( WaitForBufferToEmplace );
        m_notFullCond.wait( lockQ1Emplace, [ this ]  { return !this->full(); } );
      }
      {
        LIFO_MARK_SCOPE( WaitForBufferToPop );
        q2.m_notEmptyCond.wait( lockQ2Pop, [ &q2 ] { return !q2.empty(); } );
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
  void emplaceBackFromFront( fixedSizeDequeAndMutexes< T > & q2 )
  {
    LIFO_MARK_FUNCTION;
    {
      std::lock( q2.m_frontMutex, q2.m_popMutex, m_emplaceMutex, m_backMutex );
      std::unique_lock< std::mutex > lockQ1Emplace( m_emplaceMutex, std::adopt_lock );
      std::unique_lock< std::mutex > lockQ1Back( m_backMutex, std::adopt_lock );
      std::unique_lock< std::mutex > lockQ2Pop( q2.m_popMutex, std::adopt_lock );
      std::unique_lock< std::mutex > lockQ2Front( q2.m_frontMutex, std::adopt_lock );
      m_notFullCond.wait( lockQ1Emplace, [ this ]  { return !this->full(); } );
      q2.m_notEmptyCond.wait( lockQ2Pop, [ &q2 ] { return !q2.empty(); } );
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
class lifoStorage
{

private:

  /// number of buffers to be inserted into the LIFO
  int m_maxNumberOfBuffers;
  /// size of one buffer in bytes
  size_t m_bufferSize;
  /// name used to store data on disk
  std::string m_name;
  /// ueue of data stored on device
  fixedSizeDequeAndMutexes< T > m_deviceDeque;
  /// ueue of data stored on host memory
  fixedSizeDequeAndMutexes< T > m_hostDeque;

  /// counter of buffer stored in LIFO
  int m_bufferCount;
  /// counter of buffer stored on disk
  int m_bufferOnDiskCount;

  // Events associated to ith  copies to device buffer
  std::vector< camp::resources::Event > m_pushToDeviceEvents;
  // Futures associated to push to host in case we have no device buffers
  std::vector< std::future< void > > m_pushToHostFutures;
  // Events associated to ith  copies from device buffer
  std::vector< camp::resources::Event > m_popFromDeviceEvents;
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
    m_deviceDeque( numberOfBuffersToStoreOnDevice, elemCnt, LvArray::MemorySpace::cuda ),
    m_hostDeque( numberOfBuffersToStoreOnHost, elemCnt, LvArray::MemorySpace::host ),
    m_bufferCount( 0 ), m_bufferOnDiskCount( 0 ),
    m_pushToDeviceEvents( (numberOfBuffersToStoreOnDevice > 0)?maxNumberOfBuffers:0 ),
    m_pushToHostFutures( (numberOfBuffersToStoreOnDevice > 0)?0:maxNumberOfBuffers ),
    m_popFromDeviceEvents( (numberOfBuffersToStoreOnDevice > 0)?maxNumberOfBuffers:0 ),
    m_popFromHostFutures( (numberOfBuffersToStoreOnDevice > 0)?0:maxNumberOfBuffers )
  {
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
    GEOSX_ASSERT( m_deviceStorage.empty() && m_hostStorage.empty() && m_diskStorage.empty() && m_task_queue.empty() );
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

    if( m_deviceDeque.capacity() > 0 )
    {
      m_pushToDeviceEvents[id] = m_deviceDeque.emplaceFront( array );

      if( m_maxNumberOfBuffers - id > m_deviceDeque.capacity() )
      {
        LIFO_MARK_SCOPE( geosx::lifoStorage< T >::pushAddTasks );
        // This buffer will go to host memory, and maybe on disk
        std::packaged_task< void() > task( std::bind( &lifoStorage< T >::deviceToHost, this, id ) );
        std::unique_lock< std::mutex > lock( m_task_queue_mutex[0] );
        m_task_queue[0].emplace_back( std::move( task ) );
        lock.unlock();
        m_task_queue_not_empty_cond[0].notify_all();
      }
    }
    else
    {
      std::packaged_task< void() > task( std::bind( [ this ] ( int id, arrayView1d< T > array ) {
        m_hostDeque.emplaceFront( array );

        if( m_maxNumberOfBuffers - id > m_hostDeque.capacity() )
        {
          LIFO_MARK_SCOPE( geosx::lifoStorage< T >::pushAddTasks );
          // This buffer will go to host memory, and maybe on disk
          std::packaged_task< void() > task( std::bind( &lifoStorage< T >::hostToDisk, this, id ) );
          std::unique_lock< std::mutex > lock( m_task_queue_mutex[1] );
          m_task_queue[1].emplace_back( std::move( task ) );
          lock.unlock();
          m_task_queue_not_empty_cond[1].notify_all();
        }
      }, id, array ) );
      m_pushToHostFutures[id] = task.get_future();
      std::unique_lock< std::mutex > lock( m_task_queue_mutex[0] );
      m_task_queue[0].emplace_back( std::move( task ) );
      lock.unlock();
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
      if( m_deviceDeque.capacity() > 0 )
      {
        m_pushToDeviceEvents[m_bufferCount-1].wait();
      }
      else
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
    if( m_deviceDeque.capacity() > 0 )
    {
      m_popFromDeviceEvents[id] = m_deviceDeque.popFront( array );

      if( id >= m_deviceDeque.capacity() )
      {
        LIFO_MARK_SCOPE( geosx::lifoStorage< T >::popAddTasks );
        // Trigger pull one buffer from host, and maybe from disk
        std::packaged_task< void() > task( std::bind( &lifoStorage< T >::hostToDevice, this, id - m_deviceDeque.capacity() ) );
        std::unique_lock< std::mutex > lock( m_task_queue_mutex[0] );
        m_task_queue[0].emplace_back( std::move( task ) );
        lock.unlock();
        m_task_queue_not_empty_cond[0].notify_all();
      }
    }
    else
    {
      std::packaged_task< void() > task( std::bind ( [ this ] ( int id, arrayView1d< T > array ) {
        m_hostDeque.popFront( array );

        if( id >= m_hostDeque.capacity() )
        {
          LIFO_MARK_SCOPE( geosx::lifoStorage< T >::popAddTasks );
          // Trigger pull one buffer from host, and maybe from disk
          std::packaged_task< void() > task( std::bind( &lifoStorage< T >::diskToHost, this, id  - m_hostDeque.capacity() ) );
          std::unique_lock< std::mutex > lock( m_task_queue_mutex[1] );
          m_task_queue[1].emplace_back( std::move( task ) );
          lock.unlock();
          m_task_queue_not_empty_cond[1].notify_all();
        }
      }, id, array ) );
      m_popFromHostFutures[id] = task.get_future();
      std::unique_lock< std::mutex > lock( m_task_queue_mutex[0] );
      m_task_queue[0].emplace_back( std::move( task ) );
      lock.unlock();
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
      if( m_deviceDeque.capacity() > 0 )
      {
        cudaEventSynchronize( m_popFromDeviceEvents[m_bufferCount].get< camp::resources::CudaEvent >().getCudaEvent_t() );
      }
      else
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
  void deviceToHost( int id )
  {
    LIFO_MARK_FUNCTION;
    // The copy to host will only start when the data is copied on device buffer
    m_hostDeque.getStream().wait_for( const_cast< camp::resources::Event * >( &m_pushToDeviceEvents[id] ) );
    m_hostDeque.emplaceFrontFromBack( m_deviceDeque );

    if( m_maxNumberOfBuffers - id > m_deviceDeque.capacity() + m_hostDeque.capacity() )
    {
      // This buffer will go to host then maybe to disk
      std::packaged_task< void() > task( std::bind( &lifoStorage< T >::hostToDisk, this, id ) );
      std::unique_lock< std::mutex > lock( m_task_queue_mutex[1] );
      m_task_queue[1].emplace_back( std::move( task ) );
      lock.unlock();
      m_task_queue_not_empty_cond[1].notify_all();
    }
  }

  /**
   * Copy data from host memory to disk
   *
   * @param id ID of the buffer to store on disk.
   */
  void hostToDisk( int id )
  {
    LIFO_MARK_FUNCTION;
    {
      std::lock( m_hostDeque.m_popMutex, m_hostDeque.m_backMutex );
      std::unique_lock< std::mutex > l1( m_hostDeque.m_popMutex, std::adopt_lock );
      std::unique_lock< std::mutex > l2( m_hostDeque.m_backMutex, std::adopt_lock );
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
  void hostToDevice( int id )
  {
    LIFO_MARK_FUNCTION;
    m_hostDeque.getStream().wait_for( const_cast< camp::resources::Event * >( &m_popFromDeviceEvents[ id + m_deviceDeque.capacity() ] ) );
    m_deviceDeque.emplaceBackFromFront( m_hostDeque );

    // enqueue diskToHost on worker #2 if needed
    if( id >= m_hostDeque.capacity() )
    {
      // This buffer will go to host then to disk
      std::packaged_task< void() > task( std::bind( &lifoStorage< T >::diskToHost, this, id - m_hostDeque.capacity() ) );
      std::unique_lock< std::mutex > lock( m_task_queue_mutex[1] );
      m_task_queue[1].emplace_back( std::move( task ) );
      lock.unlock();
      m_task_queue_not_empty_cond[1].notify_all();
    }
  }

  /**
   * Copy data from disk to host memory
   *
   * @param id ID of the buffer to read on disk.
   */
  void diskToHost( int id )
  {
    LIFO_MARK_FUNCTION;
    {
      std::lock( m_hostDeque.m_emplaceMutex, m_hostDeque.m_backMutex );
      std::unique_lock< std::mutex > l1( m_hostDeque.m_emplaceMutex, std::adopt_lock );
      std::unique_lock< std::mutex > l2( m_hostDeque.m_backMutex, std::adopt_lock );
      m_hostDeque.m_notFullCond.wait( l1, [ this ]  { return !( m_hostDeque.full() ); } );
      readOnDisk( const_cast< T * >(m_hostDeque.next_back().dataIfContiguous()), id );
      m_hostDeque.inc_back();
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

    std::ofstream outfile;

    int const rank = MpiWrapper::initialized()?MpiWrapper::commRank( MPI_COMM_GEOSX ):0;
    std::string fileName = GEOSX_FMT( "{}_{:08}_{:04}.dat", m_name, id, rank );
    int lastDirSeparator = fileName.find_last_of( "/\\" );
    std::string dirName = fileName.substr( 0, lastDirSeparator );
    if( string::npos != lastDirSeparator && !dirExists( dirName ))
      makeDirsForPath( dirName );
    {
      LIFO_MARK_SCOPE( ofstreamWrite );
      const int fileDesc = open( fileName.c_str(), O_CREAT | O_WRONLY | O_DIRECT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH );
      GEOSX_ERROR_IF( fileDesc == -1,
                      "Could not open file "<< fileName << " for writting: " << strerror( errno ) );
      write( fileDesc, (char *)d, m_bufferSize );
      close( fileDesc );
    }
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
    int const rank = MpiWrapper::initialized()?MpiWrapper::commRank( MPI_COMM_GEOSX ):0;
    std::string fileName = GEOSX_FMT( "{}_{:08}_{:04}.dat", m_name, id, rank );
    const int fileDesc = open( fileName.c_str(), O_RDONLY | O_DIRECT );
    GEOSX_ERROR_IF( fileDesc == -1,
                    "Could not open file "<< fileName << " for reading: " << strerror( errno ) );
    read( fileDesc, (char *)d, m_bufferSize );
    close( fileDesc );
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
      {
        LIFO_MARK_SCOPE( runningTask );
        task();
      }
    }
  }
};
}
#endif // LIFOSTORAGE_HPP
