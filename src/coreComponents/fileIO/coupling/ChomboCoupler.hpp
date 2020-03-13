/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_FILEIO_COUPLING_CHOMBOCOUPLER_HPP_
#define GEOSX_FILEIO_COUPLING_CHOMBOCOUPLER_HPP_

#include "common/DataTypes.hpp"
#include "mesh/MeshLevel.hpp"

#include <string>

namespace geosx
{

class ChomboCoupler
{
public:

  ChomboCoupler( MPI_Comm const comm, const std::string & outputPath, const std::string & inputPath, MeshLevel & mesh );

  void write( double dt );

  void read( bool usePressures );

private:
  MPI_Comm const m_comm;
  std::string const m_outputPath;
  std::string const m_inputPath;
  std::int64_t m_face_offset;
  std::int64_t m_n_faces_written;
  std::int64_t m_node_offset;
  std::int64_t m_n_nodes_written;
  MeshLevel & m_mesh;
  int m_counter;
};

} /* namespace geosx */

#endif /* GEOSX_FILEIO_COUPLING_CHOMBOCOUPLER_HPP_ */
