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

#ifndef GEOSX_MESH_MPICOMMUNICATIONS_MPI_ICOMMDATA_HPP_
#define GEOSX_MESH_MPICOMMUNICATIONS_MPI_ICOMMDATA_HPP_

#include "CommID.hpp"

#include "common/DataTypes.hpp"

namespace geosx
{

/**
 * Class to manage the MPI requests and status data for a collection
 * of neighbor communication pipelines.
 */
class MPI_iCommData
{
public:

  /**
   * Constructor
   * @param inputCommID The CommID integer that indicates what communication
   *   pipeline to use for a set of neighbor communications.
   */
  MPI_iCommData( int const inputCommID );

  /// Default destructor
  ~MPI_iCommData() = default;

  /**
   * Resize all the arrays of requests and status'.
   * @param numMessages The number of messages/neighbors to manage requests
   *   and statuses.
   */
  void resize( localIndex const numMessages );

  /**
   *
   * @return The size of allocation (i.e. number of request/status allocated)
   */
  int size() const { return m_size; }

  /**
   *
   * @return The integer commID for this set of communications
   */
  int commID() const { return m_commID; }

  /**
   * @return Reference to the field names registered with the communication data.
   */
  std::map< string, string_array > const & getFieldNames() const { return m_fieldNames;}

  /**
   * @return Reference to the field names registered with the communication data.
   */
  std::map< string, string_array > & fieldNames() { return m_fieldNames;}


  MPI_Request * mpiSendBufferRequest() { return m_mpiSendBufferRequest.data(); }
  MPI_Request * mpiRecvBufferRequest() { return m_mpiRecvBufferRequest.data(); }
  MPI_Status * mpiSendBufferStatus()   { return m_mpiSendBufferStatus.data(); }
  MPI_Status * mpiRecvBufferStatus()   { return m_mpiRecvBufferStatus.data(); }
  MPI_Request * mpiSendBufferSizeRequest() { return m_mpiSendBufferSizeRequest.data(); }
  MPI_Request * mpiRecvBufferSizeRequest() { return m_mpiRecvBufferSizeRequest.data(); }
  MPI_Status * mpiSendBufferSizeStatus()   { return m_mpiSendBufferSizeStatus.data(); }
  MPI_Status * mpiRecvBufferSizeStatus()   { return m_mpiRecvBufferSizeStatus.data(); }


  MPI_Request & mpiSendBufferRequest( localIndex const idx ) { return m_mpiSendBufferRequest[idx]; }
  MPI_Request & mpiRecvBufferRequest( localIndex const idx ) { return m_mpiRecvBufferRequest[idx]; }
  MPI_Status & mpiSendBufferStatus( localIndex const idx )   { return m_mpiSendBufferStatus[idx]; }
  MPI_Status & mpiRecvBufferStatus( localIndex const idx )   { return m_mpiRecvBufferStatus[idx]; }

  MPI_Request & mpiSendBufferSizeRequest( localIndex const idx ) { return m_mpiSendBufferSizeRequest[idx]; }
  MPI_Request & mpiRecvBufferSizeRequest( localIndex const idx ) { return m_mpiRecvBufferSizeRequest[idx]; }
  MPI_Status & mpiSendBufferSizeStatus( localIndex const idx )   { return m_mpiSendBufferSizeStatus[idx]; }
  MPI_Status & mpiRecvBufferSizeStatus( localIndex const idx )   { return m_mpiRecvBufferSizeStatus[idx]; }

private:
  /// The number of message pipelines to manage
  int m_size;

  /// The integer ID for the set of communication pipelines
  int m_commID;

  /// A collection of field names keyed on object keys to pack/unpack from
  /// communication pipeline.
  std::map< string, string_array > m_fieldNames;

  array1d< MPI_Request > m_mpiSendBufferRequest;
  array1d< MPI_Request > m_mpiRecvBufferRequest;
  array1d< MPI_Status >  m_mpiSendBufferStatus;
  array1d< MPI_Status >  m_mpiRecvBufferStatus;
  array1d< MPI_Request > m_mpiSendBufferSizeRequest;
  array1d< MPI_Request > m_mpiRecvBufferSizeRequest;
  array1d< MPI_Status >  m_mpiSendBufferSizeStatus;
  array1d< MPI_Status >  m_mpiRecvBufferSizeStatus;
};
} /* namespace geosx */

#endif /* GEOSX_MESH_MPICOMMUNICATIONS_MPI_ICOMMDATA_HPP_ */
