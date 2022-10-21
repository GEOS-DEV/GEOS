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

#if defined( GEOSX_USE_CUDA )
#  include <cuda.h>
#  include <cuda_runtime.h>
#endif

#include <chai/ArrayManager.hpp>
#include <chai/ChaiMacros.hpp>
#include "common/fixedSizeDeque.hpp"
#include <deque>
#include <future>
#include <string>
#include <functional>
#include "LvArray/src/Array.hpp"
#include "LvArray/src/ArrayView.hpp"
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
  /**
   * Handler to identify and store data stored inn the LIFO.
   */
  struct storageData
  {
    char * m_data;
  };

  /// Number of buffers to be inserted into the LIFO
  int m_maxNumberOfBuffers;
  /// Number of buffers that will be stored on device (at least 1)
  int m_numberOfBuffersToStoreOnDevice;
  /// Number of buffers that will be stored on host memory (at least 1)
  int m_numberOfBuffersToStoreOnHost;
  /// Size of one buffer in bytes
  size_t m_bufferSize;
  /// Name used to store data on disk
  std::string m_name;
  /// Queue of data stored on device
  std::deque< storageData > m_deviceStorage;
  fixedSizeDeque< T, int > m_deviceDeque;

  /// Queue of data stored on host memory
  std::deque< storageData > m_hostStorage;
  /// Queue of data stored on disk
  std::deque< storageData > m_diskStorage;
  /**
   * Mutex to protect access to m_deviceStorage.
   *
   * Only deviceStorage is accessed from user thread, other only from m_worker.
   */
  std::mutex m_deviceStorageMutex;
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
    m_numberOfBuffersToStoreOnDevice( numberOfBuffersToStoreOnDevice ), m_numberOfBuffersToStoreOnHost( numberOfBuffersToStoreOnHost ),
    m_bufferSize( elemCnt*sizeof( T ) ),
    m_name( name ),
    m_deviceDeque( maxNumberOfBuffers, elemCnt, LvArray::MemorySpace::cuda ),
    m_bufferCount( 0 ), m_bufferOnDiskCount( 0 ),
    m_worker( &lifoStorage< T >::wait_and_consume_tasks, this ),
    m_continue( true )
  {
    GEOSX_ASSERT( numberOfBuffersToStoreOnDevice > 0 && numberOfBuffersToStoreOnHost >= 0 && maxNumberOfBuffers >= numberOfBuffersToStoreOnHost * numberOfBuffersToStoreOnDevice );
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
   * Push a copy of the given LvArray into the LIFO
   *
   * @param array The LvArray to store in the LIFO, should match the size of the data used in constructor.
   */
  void push( arrayView1d< T > array )
  {
    std::cout <<__FILE__<<":"<<__LINE__<<" array[i] " <<  array[1] << std::endl;
    array.move( LvArray::MemorySpace::cuda, false );
    m_deviceDeque.emplace_back( array.toSliceConst() );
    push( array.data(), array.size() );
  }

  /**
   * push elemCnt T elements from buffer data into the LIFO.
   *
   * @param data    buffer to store data from in the LIFO.
   * @param elemCnt Number of element to store in the LIFO, should match the one used in constructor.
   */
  void push( T * data, size_t elemCnt )
  {
    GEOSX_ASSERT( m_bufferSize == sizeof(T)*elemCnt; );

    storageData d;
    int id = m_bufferCount++;
    chai::gpuMalloc((void * *)&d.m_data, m_bufferSize );
    chai::gpuMemcpy( d.m_data, data, m_bufferSize, gpuMemcpyDeviceToDevice );
    m_deviceStorageMutex.lock();
    while( m_deviceStorage.size() >= (size_t)m_numberOfBuffersToStoreOnDevice )
    {
      m_deviceStorageMutex.unlock();
      m_deviceToHostFutures.front().wait();
      m_deviceToHostFutures.pop_front();
      m_deviceStorageMutex.lock();
    }
    m_deviceStorage.emplace_front( d );
    m_deviceStorageMutex.unlock();

    if( m_maxNumberOfBuffers - id > m_numberOfBuffersToStoreOnDevice + m_numberOfBuffersToStoreOnHost )
    {
      // This buffer will go to host then to disk
      std::packaged_task< void() > task1 ( std::bind( &lifoStorage< T >::deviceToHost, this ) );
      m_deviceToHostFutures.emplace_back( task1.get_future() );
      m_task_queue.emplace_back( std::move( task1 ) );
      m_task_queue_not_empty_cond.notify_all();
      std::packaged_task< void() > task = std::packaged_task< void() >( std::bind( &lifoStorage< T >::hostToDisk, this ) );
      m_hostToDiskFutures.emplace_back( task.get_future() );
      m_task_queue.emplace_back( std::move( task ) );
      m_task_queue_not_empty_cond.notify_all();
    }
    else if( m_maxNumberOfBuffers - id > m_numberOfBuffersToStoreOnDevice )
    {
      // This buffer will go to host memory
      std::packaged_task< void() > task = std::packaged_task< void() >( std::bind( &lifoStorage< T >::deviceToHost, this ) );
      m_deviceToHostFutures.emplace_back( task.get_future() );
      m_task_queue.emplace_back( std::move( task ) );
      m_task_queue_not_empty_cond.notify_all();
    }
  }

  /**
   * Copy last data from the LIFO into the LvArray.
   *
   * @param array LvArray to store data from the LIFO into it.
   */
  void pop( arrayView1d< T > array )
  {
    array.move( LvArray::MemorySpace::cuda, true );
    pop( array.data(), array.size() );
    LvArray::memcpy( array.toSlice(), m_deviceDeque.front( ) );
    m_deviceDeque.pop_back();
  }

  /**
   * Copylast data from the LIFO into the given buffer
   *
   * @param data    The buffer to copy data into.
   * @param elemCnt Number of elements to store into data, should match size used in constructor.
   */
  void pop( T * data, size_t elemCnt )
  {
    GEOSX_ASSERT( m_bufferSize == sizeof(T)*elemCnt; );

    int id = --m_bufferCount;

    m_deviceStorageMutex.lock();
    while( m_deviceStorage.empty() )
    {
      m_deviceStorageMutex.unlock();
      m_hostToDeviceFutures.front().wait();
      m_hostToDeviceFutures.pop_front();
      m_deviceStorageMutex.lock();
    }
    storageData d = m_deviceStorage.front();
    m_deviceStorage.pop_front();
    m_deviceStorageMutex.unlock();

    chai::gpuMemcpy( data, d.m_data, m_bufferSize, gpuMemcpyDeviceToDevice );
    std::packaged_task< void() > task;
    if( id >= m_numberOfBuffersToStoreOnDevice + m_numberOfBuffersToStoreOnHost )
    {
      // Trigger pull one buffer from disk
      task = std::packaged_task< void() >( std::bind( &lifoStorage< T >::diskToHost, this ) );
      m_diskToHostFutures.emplace_back( task.get_future() );
      m_task_queue.emplace_back( std::move( task ) );
      m_task_queue_not_empty_cond.notify_all();
    }
    if( id >= m_numberOfBuffersToStoreOnDevice )
    {
      // Trigger pull one buffer from host
      task = std::packaged_task< void() >( std::bind( &lifoStorage< T >::hostToDevice, this ) );
      m_hostToDeviceFutures.emplace_back( task.get_future() );
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
    m_deviceStorageMutex.lock();
    storageData d = m_deviceStorage.back();
    m_deviceStorage.pop_back();

    m_deviceStorageMutex.unlock();

    char * ptr = new char[m_bufferSize];
    chai::gpuMemcpy( ptr, d.m_data, m_bufferSize, gpuMemcpyDeviceToHost );
    chai::gpuFree( d.m_data );
    d.m_data = ptr;

    while( m_hostStorage.size() >= (size_t)m_numberOfBuffersToStoreOnHost )
    {
      m_hostToDiskFutures.front().wait();
      m_hostToDiskFutures.pop_front();
    }
    m_hostStorage.emplace_front( d );
  }

  /**
   * Copy data from host memory to disk
   */
  void hostToDisk()
  {
    storageData d = std::move( m_hostStorage.back() );
    m_hostStorage.pop_back();
    writeOnDisk( d );
    delete[] d.m_data;
    d.m_data = nullptr;

    m_diskStorage.emplace_front( d );
  }

  /**
   * Copy data from host memory to device memory
   */
  void hostToDevice()
  {
    storageData d = std::move( m_hostStorage.front());
    m_hostStorage.pop_front();
    char * ptr;
    chai::gpuMalloc((void * *)&ptr, m_bufferSize );
    chai::gpuMemcpy( ptr, d.m_data, m_bufferSize, gpuMemcpyHostToDevice );
    delete[] d.m_data;
    d.m_data = ptr;
    m_deviceStorageMutex.lock();
    m_deviceStorage.emplace_back( d );
    m_deviceStorageMutex.unlock();
  }

  /**
   * Copy data from disk to host memory
   */
  void diskToHost()
  {
    storageData d = std::move( m_diskStorage.front() );
    m_diskStorage.pop_front();
    readOnDisk( d );
    m_hostStorage.emplace_back( d );
  }

  /**
   * Write data on disk
   *
   * @param d Data to store on disk.
   */
  void writeOnDisk( storageData d )
  {
    int id = m_bufferOnDiskCount++;
    int const rank = MpiWrapper::initialized()?MpiWrapper::commRank( MPI_COMM_GEOSX ):0;
    std::string fileName = GEOSX_FMT( "{}_{:08}_{:04}.dat", m_name, id, rank );
    std::ofstream wf( fileName, std::ios::out | std::ios::binary );
    GEOSX_THROW_IF( !wf,
                    "Could not open file "<< fileName << " for writting",
                    InputError );
    wf.write( d.m_data, m_bufferSize );
    wf.close();
    GEOSX_THROW_IF( !wf.good(),
                    "An error occured while writting "<< fileName,
                    InputError );
  }

  /**
   * Read data from disk
   *
   * @param Handler to store datta read from disk.
   */
  void readOnDisk( storageData & d )
  {
    int id = --m_bufferOnDiskCount;
    int const rank = MpiWrapper::initialized()?MpiWrapper::commRank( MPI_COMM_GEOSX ):0;
    std::string fileName = GEOSX_FMT( "{}_{:08}_{:04}.dat", m_name, id, rank );
    std::ifstream wf( fileName, std::ios::in | std::ios::binary );
    GEOSX_THROW_IF( !wf,
                    "Could not open file "<< fileName << " for reading",
                    InputError );
    d.m_data = new char[m_bufferSize];
    wf.read( d.m_data, m_bufferSize );
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
