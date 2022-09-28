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
  struct storageData
  {
    char * m_data;
    int m_id;
  };

  int m_maxNumberOfBuffers;
  int m_numberOfBuffersToStoreOnDevice;
  int m_numberOfBuffersToStoreOnHost;
  size_t m_bufferSize;
  std::string m_name;
  std::deque< storageData > m_deviceStorage;
  std::deque< storageData > m_hostStorage;
  std::deque< storageData > m_diskStorage;
  // only deviceStorage is accessed from user thread, other only from m_worker.
  std::mutex m_deviceStorageMutex;
  std::deque< std::future< void > > m_deviceToHostFutures;
  std::deque< std::future< void > > m_hostToDiskFutures;
  std::deque< std::future< void > > m_diskToHostFutures;
  std::deque< std::future< void > > m_hostToDeviceFutures;
  int m_bufferCount;

  std::condition_variable m_task_queue_not_empty_cond;
  std::mutex m_task_queue_mutex;
  std::deque< std::packaged_task< void() > > m_task_queue;
  std::thread m_worker;
  bool m_continue;

public:
  /**
   * A LIFO storage will store numberOfBuffersToStoreDevice buffer on
   * device, numberOfBuffersToStoreHost on host and the rest on disk.
   */
  lifoStorage( std::string name, size_t elemCnt, int numberOfBuffersToStoreOnDevice, int numberOfBuffersToStoreOnHost, int maxNumberOfBuffers ):
    m_maxNumberOfBuffers( maxNumberOfBuffers ),
    m_numberOfBuffersToStoreOnDevice( numberOfBuffersToStoreOnDevice ), m_numberOfBuffersToStoreOnHost( numberOfBuffersToStoreOnHost ),
    m_bufferSize( elemCnt*sizeof( T ) ),
    m_name( name ), m_bufferCount( 0 ), m_worker( &lifoStorage< T >::wait_and_consume_tasks, this ),
    m_continue( true )
  {
    assert( numberOfBuffersToStoreOnDevice > 0 && numberOfBuffersToStoreOnHost >= 0 && maxNumberOfBuffers >= numberOfBuffersToStoreOnHost * numberOfBuffersToStoreOnDevice );
  }

  lifoStorage( std::string name, arrayView1d< T > array, int numberOfBuffersToStoreOnDevice, int numberOfBuffersToStoreOnHost, int maxNumberOfBuffers ):
    lifoStorage( name, array.size(), numberOfBuffersToStoreOnDevice, numberOfBuffersToStoreOnHost, maxNumberOfBuffers ) {}

  ~lifoStorage()
  {
    m_continue = false;
    m_task_queue_not_empty_cond.notify_all();
    m_worker.join();
    assert( m_deviceStorage.empty() && m_hostStorage.empty() && m_diskStorage.empty() && m_task_queue.empty() );
  }

  void push( arrayView1d< T > array )
  {
    array.move( LvArray::MemorySpace::cuda, false );
    push( array.data(), array.size() );
  }

  void push( T * data, size_t elemCnt )
  {
    assert( m_bufferSize == sizeof(T)*elemCnt; );

    storageData d;
    d.m_id = m_bufferCount++;
    chai::gpuMalloc((void * *)&d.m_data, m_bufferSize );
    chai::gpuMemcpy( d.m_data, data, m_bufferSize, gpuMemcpyDeviceToDevice );

    m_deviceStorageMutex.lock();
    while( m_deviceStorage.size() >= m_numberOfBuffersToStoreOnDevice )
    {
      m_deviceStorageMutex.unlock();
      m_deviceToHostFutures.front().wait();
      m_deviceToHostFutures.pop_front();
      m_deviceStorageMutex.lock();
    }
    m_deviceStorage.emplace_front( d );
    m_deviceStorageMutex.unlock();

    if( m_maxNumberOfBuffers - d.m_id > m_numberOfBuffersToStoreOnDevice + m_numberOfBuffersToStoreOnHost )
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
    else if( m_maxNumberOfBuffers - d.m_id > m_numberOfBuffersToStoreOnDevice )
    {
      // This buffer will go to host memory
      std::packaged_task< void() > task = std::packaged_task< void() >( std::bind( &lifoStorage< T >::deviceToHost, this ) );
      m_deviceToHostFutures.emplace_back( task.get_future() );
      m_task_queue.emplace_back( std::move( task ) );
      m_task_queue_not_empty_cond.notify_all();
    }
  }

  void pop( arrayView1d< T > array )
  {
    array.move( LvArray::MemorySpace::cuda, true );
    pop( array.data(), array.size() );
  }

  void pop( T * data, size_t elemCnt )
  {
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
    if( d.m_id >= m_numberOfBuffersToStoreOnDevice + m_numberOfBuffersToStoreOnHost )
    {
      // Trigger pull one buffer from disk
      task = std::packaged_task< void() >( std::bind( &lifoStorage< T >::diskToHost, this ) );
      m_diskToHostFutures.emplace_back( task.get_future() );
      m_task_queue.emplace_back( std::move( task ) );
      m_task_queue_not_empty_cond.notify_all();
    }
    if( d.m_id >= m_numberOfBuffersToStoreOnDevice )
    {
      // Trigger pull one buffer from host
      task = std::packaged_task< void() >( std::bind( &lifoStorage< T >::hostToDevice, this ) );
      m_hostToDeviceFutures.emplace_back( task.get_future() );
      m_task_queue.emplace_back( std::move( task ) );
      m_task_queue_not_empty_cond.notify_all();
    }
  }

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

    while( m_hostStorage.size() >= m_numberOfBuffersToStoreOnHost )
    {
      m_hostToDiskFutures.front().wait();
      m_hostToDiskFutures.pop_front();
    }
    m_hostStorage.emplace_front( d );
  }

  void hostToDisk()
  {
    storageData d = std::move( m_hostStorage.back() );
    m_hostStorage.pop_back();
    writeOnDisk( d );
    delete[] d.m_data;
    d.m_data = nullptr;

    m_diskStorage.emplace_front( d );
  }

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

  void diskToHost()
  {
    storageData d = std::move( m_diskStorage.front() );
    m_diskStorage.pop_front();
    readOnDisk( d );
    m_hostStorage.emplace_back( d );
  }

  void writeOnDisk( storageData d )
  {
    int const rank = MpiWrapper::initialized()?MpiWrapper::commRank( MPI_COMM_GEOSX ):0;
    std::string fileName = GEOSX_FMT( "{}_{:08}_{:04}.dat", m_name, d.m_id, rank );
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

  void readOnDisk( storageData & d )
  {
    int const rank = MpiWrapper::initialized()?MpiWrapper::commRank( MPI_COMM_GEOSX ):0;
    std::string fileName = GEOSX_FMT( "{}_{:08}_{:04}.dat", m_name, d.m_id, rank );
    std::ifstream wf( fileName, std::ios::in | std::ios::binary );
    std::cerr << "Reading " << fileName << std::endl;
    GEOSX_THROW_IF( !wf,
                    "Could not open file "<< fileName << " for reading",
                    InputError );
    d.m_data = new char[m_bufferSize];
    wf.read( d.m_data, m_bufferSize );
    wf.close();
    remove( fileName.c_str() );
  }

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
