/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2020-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


#include "TestingBuffer.hpp"

#include <numeric>

namespace geosx
{

TestingBuffer::TestingBuffer( string name, std::vector< localIndex > const & dims ):
  m_name( name ),
  m_initialDims( dims ),
  m_count( 1 ),
  m_numStorages( 0 )
{

}

buffer_unit_type * TestingBuffer::getBufferHead()
{
  std::vector< localIndex > tmp = m_initialDims;
  tmp[0] = m_count;
  localIndex const n = std::accumulate( tmp.cbegin(), tmp.cend(), 1, std::multiplies< localIndex >() );
  m_dims.push_back( tmp );

  m_buffers.emplace_back( n, 0 );

  return reinterpret_cast<buffer_unit_type *>(m_buffers.back().data());
}

void TestingBuffer::init( bool existsOkay )
{
  GEOSX_LOG( "TestingBuffer::init" );
}

void TestingBuffer::write()
{
  if( m_name == "fluid_phaseDensity" )
  {
    GEOSX_ASSERT( m_dims.size() == m_buffers.size() );

    for( size_t collection = 0; collection < m_dims.size(); ++collection )
    {
      GEOSX_LOG( "TestingBuffer::write " + m_name + " for buffer " + std::to_string( m_numStorages + collection ) );
      std::vector< localIndex > const & dims = m_dims[collection];
      std::vector< double > const & buffer = m_buffers[collection];
      GEOSX_ASSERT( dims.size() == 2 );

      for( size_t tot = 0, i = 0; i < size_t( dims[0] ); ++i )
      {
        for( size_t j = 0; j < size_t( dims[1] ); ++j, ++tot )
        {
          std::cout << buffer[tot] << " ";
        }
        std::cout << std::endl;
      }
    }

    m_numStorages += m_dims.size();
    m_dims.clear();
    m_buffers.clear();
  }
}

void TestingBuffer::compressInFile()
{

}

void TestingBuffer::updateCollectingCount( localIndex count )
{
  m_count = count;
}

}