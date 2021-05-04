/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_MESH_MPICOMMUNICATIONS_MPI_ICOMMDATA_HPP_
#define GEOSX_MESH_MPICOMMUNICATIONS_MPI_ICOMMDATA_HPP_

#include "CommID.hpp"

#include "common/DataTypes.hpp"

namespace geosx
{

class MPI_iCommData
{
public:

  MPI_iCommData( int const commIdForBuffer,
                 int const commIdForSize );

  ~MPI_iCommData() = default;

  void resize( localIndex const numMessages );

  int size;
  int commID;
  int sizeCommID;
  std::map< string, string_array > fieldNames;

  array1d< MPI_Request > mpiSendBufferRequest;
  array1d< MPI_Request > mpiRecvBufferRequest;
  array1d< MPI_Status >  mpiSendBufferStatus;
  array1d< MPI_Status >  mpiRecvBufferStatus;

  array1d< MPI_Request > mpiSizeSendBufferRequest;
  array1d< MPI_Request > mpiSizeRecvBufferRequest;
  array1d< MPI_Status >  mpiSizeSendBufferStatus;
  array1d< MPI_Status >  mpiSizeRecvBufferStatus;
};
} /* namespace geosx */

#endif /* GEOSX_MESH_MPICOMMUNICATIONS_MPI_ICOMMDATA_HPP_ */
