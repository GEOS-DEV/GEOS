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

#ifndef GEOSX_TESTINGBUFFER_HPP
#define GEOSX_TESTINGBUFFER_HPP

#include "HistoryIO.hpp"

namespace geosx
{

class TestingBuffer: public BufferedHistoryIO
{
public:

  TestingBuffer( string name, std::vector< localIndex > const & dims );

  buffer_unit_type * getBufferHead() override;

  void init( bool existsOkay ) override;

  void write() override;

  void compressInFile() override;

  void updateCollectingCount( localIndex count ) override;

private:
  string const m_name;

  std::vector< localIndex > const m_initialDims;

  std::vector< std::vector< localIndex > > m_dims;

  localIndex m_count;

  std::vector< std::vector< double > > m_buffers;

  localIndex m_numStorages;
};

}



#endif // include guard
