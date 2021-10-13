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

#include "MPI_iCommData.hpp"

namespace geosx
{

MPI_iCommData::MPI_iCommData( int const inputCommID ):
  m_size( 0 ),
  m_commID( inputCommID ),      // CommunicationTools::getInstance().getCommID() ),
  m_fieldNames(),
  m_mpiSendBufferRequest(),
  m_mpiRecvBufferRequest(),
  m_mpiSendBufferStatus(),
  m_mpiRecvBufferStatus(),
  m_mpiSendBufferSizeRequest(),
  m_mpiRecvBufferSizeRequest(),
  m_mpiSendBufferSizeStatus(),
  m_mpiRecvBufferSizeStatus()
{}

void MPI_iCommData::resize( localIndex numMessages )
{
  m_mpiSendBufferRequest.resize( numMessages );
  m_mpiRecvBufferRequest.resize( numMessages );
  m_mpiSendBufferStatus.resize( numMessages );
  m_mpiRecvBufferStatus.resize( numMessages );
  m_mpiSendBufferSizeRequest.resize( numMessages );
  m_mpiRecvBufferSizeRequest.resize( numMessages );
  m_mpiSendBufferSizeStatus.resize( numMessages );
  m_mpiRecvBufferSizeStatus.resize( numMessages );
  m_size = static_cast< int >(numMessages);
}

} /* namespace geosx */
