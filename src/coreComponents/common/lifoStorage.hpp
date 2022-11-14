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

#include "common/fixedSizeDeque.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{
/**
 * This class is used to store in a LIFO way buffers, first on device, then on host, then on disk.
 */
template< typename T >
class lifoStorage
{

private:

  /// Number of buffers to be inserted into the LIFO
  int m_maxNumberOfBuffers;
  /// Size of one buffer in bytes
  size_t m_bufferSize;
  /// Name used to store data on disk
  std::string m_name;
  /// Queue of data stored on device
  fixedSizeDeque< T, int > m_deviceDeque;

  /// Queue of data stored on host memory
  fixedSizeDeque< T, int > m_hostDeque;

  /**
   * Mutex to protect access to m_deviceDeque
   *
   * Only deviceDeque is accessed from user thread, other only from m_worker.
   */
  std::mutex m_deviceDequeMutex;
  /**
   * Future to be aware data are effectively copied on host memory from device.
   */
  std::deque< std::future< void > > m_deviceToHostFutures;
  /**
   * Future to be aware data are effectively copied on disk from host memory.
   */
  std::deque< std::future< void > > m_hostToDiskFutures;
  /**
   * Future to be aware data are effectively copied on host memory from disk.
   */
  std::deque< std::future< void > > m_diskToHostFutures;
  /**
   * Future to be aware data are effectively copied on device from host memory.
   */
  std::deque< std::future< void > > m_hostToDeviceFutures;
  /**
   * Future to be aware data push has been completed
   */
  std::future< void > m_pushFuture;
  /**
   * Future to be aware data pop has been completed
   */
  std::future< void > m_popFuture;

  /**
   * Counter of buffer stored in LIFO
   */
  int m_bufferCount;

  /**
   * Counter of buffer stored on disk
   */
  int m_bufferOnDiskCount;

  /// Condition used to tell m_worker queue has been filled or processed is stopped.
  std::condition_variable m_task_queue_not_empty_cond;
  /// Mutex to protect access to m_task_queue.
  std::mutex m_task_queue_mutex;
  /// Queue of task to be executed by m_worker.
  std::deque< std::packaged_task< void() > > m_task_queue;
  /// Thread to execute tasks.
  std::thread m_worker;
  /// Boolean to keep m_worker alive.
  bool m_continue;

public:
  /**
   * A LIFO storage will store numberOfBuffersToStoreDevice buffer on
   * device, numberOfBuffersToStoreHost on host and the rest on disk.
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
    m_worker( &lifoStorage< T >::wait_and_consume_tasks, this ),
    m_continue( true )
  {
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
    m_task_queue_not_empty_cond.notify_all();
    m_worker.join();
    GEOSX_ASSERT( m_deviceStorage.empty() && m_hostStorage.empty() && m_diskStorage.empty() && m_task_queue.empty() );
  }

  /**
   * Asynchroneously push a copy of the given LvArray into the LIFO
   *
   * @param array The LvArray to store in the LIFO, should match the size of the data used in constructor.
   */
  void pushAsync( arrayView1d< T > array )
  {
    pushWait();
    m_pushFuture = std::async( std::launch::async, [ array, this ] { push( array ); } );
  }

  /**
   * Waits for last push to be terminated
   */
  void pushWait()
  {
    if( m_pushFuture.valid() ) m_pushFuture.get();
  }

  /**
   * Push a copy of the given LvArray into the LIFO
   *
   * @param array The LvArray to store in the LIFO, should match the size of the data used in constructor.
   */
  void push( arrayView1d< T > array )
  {
    int id = m_bufferCount++;
    GEOSX_ERROR_IF( m_deviceDeque.capacity() == 0,
                    "Cannot save on a Lifo without device storage (please set lifoSize, lifoOnDevice and lifoOnHost in xml file)" );
    GEOSX_ERROR_IF( m_hostDeque.capacity() == 0,
                    "Cannot save on a Lifo without host storage (please set lifoSize, lifoOnDevice and lifoOnHost in xml file)" );
    m_deviceDequeMutex.lock();
    while( m_deviceDeque.full() )
    {
      m_deviceDequeMutex.unlock();
      m_deviceToHostFutures.front().wait();
      m_deviceToHostFutures.pop_front();
      m_deviceDequeMutex.lock();
    }
    m_deviceDeque.emplace_front( array.toSliceConst() );
    m_deviceDequeMutex.unlock();

    if( m_maxNumberOfBuffers - id > m_deviceDeque.capacity() )
    {
      // This buffer will go to host memory
      std::packaged_task< void() > task( std::bind( &lifoStorage< T >::deviceToHost, this ) );
      m_deviceToHostFutures.emplace_back( task.get_future() );
      m_task_queue.emplace_back( std::move( task ) );
      m_task_queue_not_empty_cond.notify_all();
    }

    if( m_maxNumberOfBuffers - id > m_deviceDeque.capacity() + m_hostDeque.capacity() )
    {
      // This buffer will go to host then to disk
      std::packaged_task< void() > task( std::bind( &lifoStorage< T >::hostToDisk, this ) );
      m_hostToDiskFutures.emplace_back( task.get_future() );
      m_task_queue.emplace_back( std::move( task ) );
      m_task_queue_not_empty_cond.notify_all();
    }
  }

  /**
   * Asynchroneously copy last data from the LIFO into the LvArray.
   *
   * @param array LvArray to store data from the LIFO into it.
   */
  void popAsync( arrayView1d< T > array )
  {
    popWait();
    m_popFuture = std::async( std::launch::async, [ array, this ] { pop( array ); } );
  }

  /**
   * Waits for last pop to be terminated
   */
  void popWait()
  {
    if( m_popFuture.valid() ) m_popFuture.get();
  }

  /**
   * Copy last data from the LIFO into the LvArray.
   *
   * @param array LvArray to store data from the LIFO into it.
   */
  void pop( arrayView1d< T > array )
  {
    int id = --m_bufferCount;

    m_deviceDequeMutex.lock();
    while( m_deviceDeque.empty() )
    {
      m_deviceDequeMutex.unlock();
      m_hostToDeviceFutures.front().wait();
      m_hostToDeviceFutures.pop_front();
      m_deviceDequeMutex.lock();
    }
    LvArray::memcpy( array.toSlice(), m_deviceDeque.front( ) );
    m_deviceDeque.pop_front();
    m_deviceDequeMutex.unlock();

    if( id >= m_deviceDeque.capacity() )
    {
      // Trigger pull one buffer from host
      std::packaged_task< void() > task( std::bind( &lifoStorage< T >::hostToDevice, this ) );
      m_hostToDeviceFutures.emplace_back( task.get_future() );
      m_task_queue.emplace_back( std::move( task ) );
      m_task_queue_not_empty_cond.notify_all();
    }

    if( id >= m_deviceDeque.capacity() + m_hostDeque.capacity() )
    {
      // Trigger pull one buffer from disk
      std::packaged_task< void() > task( std::bind( &lifoStorage< T >::diskToHost, this ) );
      m_diskToHostFutures.emplace_back( task.get_future() );
      m_task_queue.emplace_back( std::move( task ) );
      m_task_queue_not_empty_cond.notify_all();
    }
  }


private:
  /**
   * Copy data from device memory to host memory
   */
  void deviceToHost()
  {
    while( m_hostDeque.full() )
    {
      m_hostToDiskFutures.front().wait();
      m_hostToDiskFutures.pop_front();
    }
    m_hostDeque.emplace_front( m_deviceDeque.back() );
    m_deviceDequeMutex.lock();
    m_deviceDeque.pop_back();
    m_deviceDequeMutex.unlock();
  }

  /**
   * Copy data from host memory to disk
   */
  void hostToDisk()
  {
    writeOnDisk( m_hostDeque.back().dataIfContiguous() );
    m_hostDeque.pop_back();
  }

  /**
   * Copy data from host memory to device memory
   */
  void hostToDevice()
  {
    while( m_hostDeque.empty() )
    {
      m_diskToHostFutures.front().wait();
      m_diskToHostFutures.pop_front();
    }

    m_deviceDequeMutex.lock();
    m_deviceDeque.emplace_back( m_hostDeque.front() );
    m_deviceDequeMutex.unlock();
    m_hostDeque.pop_front();
  }

  /**
   * Copy data from disk to host memory
   */
  void diskToHost()
  {
    while( m_hostDeque.full() )
    {
      m_hostToDeviceFutures.front().wait();
      m_hostToDeviceFutures.pop_front();
    }
    m_hostDeque.inc_back();
    readOnDisk( const_cast< T * >(m_hostDeque.back().dataIfContiguous()) );
  }

  /**
   * Write data on disk
   *
   * @param d Data to store on disk.
   */
  void writeOnDisk( const T * d )
  {


    std::ofstream outfile;

    int id = m_bufferOnDiskCount++;
    int const rank = MpiWrapper::initialized()?MpiWrapper::commRank( MPI_COMM_GEOSX ):0;
    std::string fileName = GEOSX_FMT( "{}_{:08}_{:04}.dat", m_name, id, rank );
    makeDirsForPath( fileName.substr( 0, fileName.find_last_of( "/\\" ) ) );
    std::ofstream wf( fileName, std::ios::out | std::ios::binary );
    GEOSX_ERROR_IF( !wf || wf.fail() || !wf.is_open(),
                    "Could not open file "<< fileName << " for writting: " << strerror(errno) );
    if (wf)
      wf.write( (char *)d, m_bufferSize );
    GEOSX_ERROR_IF( wf.bad() || wf.fail(),
                    "An error occured while writting "<< fileName << ": " << strerror(errno) );
    wf.close();
  }

  /**
   * Read data from disk
   *
   * @param Handler to store datta read from disk.
   */
  void readOnDisk( T * d )
  {
    int id = --m_bufferOnDiskCount;
    int const rank = MpiWrapper::initialized()?MpiWrapper::commRank( MPI_COMM_GEOSX ):0;
    std::string fileName = GEOSX_FMT( "{}_{:08}_{:04}.dat", m_name, id, rank );
    std::ifstream wf( fileName, std::ios::in | std::ios::binary );
    GEOSX_ERROR_IF( !wf,
                    "Could not open file "<< fileName << " for reading: " << strerror(errno) );
    wf.read( (char *)d, m_bufferSize );
    wf.close();
    remove( fileName.c_str() );
  }

  /**
   * Execute the tasks pushed in m_task_queue in the given order.
   */
  void wait_and_consume_tasks()
  {
    while( m_continue )
    {
      std::unique_lock< std::mutex > lock( m_task_queue_mutex );
      m_task_queue_not_empty_cond.wait( lock, [ this ] { return !( m_task_queue.empty()  && m_continue ); } );
      if( m_continue == false ) break;
      m_task_queue.front()();
      m_task_queue.pop_front();
    }
  }
};
}
#endif // LIFOSTORAGE_HPP
