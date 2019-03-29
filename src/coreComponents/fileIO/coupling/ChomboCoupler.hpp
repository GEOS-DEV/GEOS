/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef SRC_CORECOMPONENTS_FILEIO_COUPLING_CHOMBOCOUPLER_HPP_
#define SRC_CORECOMPONENTS_FILEIO_COUPLING_CHOMBOCOUPLER_HPP_

#include "common/DataTypes.hpp"
#include "mesh/MeshLevel.hpp"

#include <string>

namespace geosx 
{

class ChomboCoupler
{
public:

  ChomboCoupler(MPI_Comm const comm, const std::string& outputPath, const std::string& inputPath, MeshLevel& mesh);

  void write(double dt);

  void read(bool usePressures);

private:
  MPI_Comm const m_comm;
  std::string const m_outputPath;
  std::string const m_inputPath;
  std::int64_t m_face_offset;
  std::int64_t m_n_faces_written;
  std::int64_t m_node_offset;
  std::int64_t m_n_nodes_written;
  MeshLevel& m_mesh;
  int m_counter;
};

} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_FILEIO_COUPLING_CHOMBOCOUPLER_HPP_ */
