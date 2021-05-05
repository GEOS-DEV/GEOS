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

#include "MPI_iCommData.hpp"

namespace geosx
{

MPI_iCommData::MPI_iCommData( int const commIdForBuffer,
                              int const commIdForSize ):
  size( 0 ),
  commID( commIdForBuffer ),      // CommunicationTools::getInstance().getCommID() ),
  sizeCommID( commIdForSize ),      //CommunicationTools::getInstance().getCommID() ),
  fieldNames(),
  mpiSendBufferRequest(),
  mpiRecvBufferRequest(),
  mpiSendBufferStatus(),
  mpiRecvBufferStatus()
{}

void MPI_iCommData::resize( localIndex numMessages )
{
  mpiSendBufferRequest.resize( numMessages );
  mpiRecvBufferRequest.resize( numMessages );
  mpiSendBufferStatus.resize( numMessages );
  mpiRecvBufferStatus.resize( numMessages );
  mpiSendBufferSizeRequest.resize( numMessages );
  mpiRecvBufferSizeRequest.resize( numMessages );
  mpiSendBufferSizeStatus.resize( numMessages );
  mpiRecvBufferSizeStatus.resize( numMessages );
  size = static_cast< int >(numMessages);
}

} /* namespace geosx */
