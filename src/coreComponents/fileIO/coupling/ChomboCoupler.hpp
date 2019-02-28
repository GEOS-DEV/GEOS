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

#include "hdf5_interface/coupler.hpp"
#include "common/DataTypes.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/FaceManager.hpp"
#include "mesh/MeshLevel.hpp"

#include <cstdint>
#include <tuple>
#include <string>

namespace geosx 
{

class ChomboCoupler
{
public:

  ChomboCoupler(MPI_Comm comm, const std::string& path, MeshLevel& mesh):
    m_comm(comm),
    m_file_path(path),
    m_face_offset(-1),
    m_n_faces_written(-1),
    m_node_offset(-1),
    m_n_nodes_written(-1),
    m_mesh(mesh),
    m_counter(0)
  {
    /* Create a dummy pressure field */
    FaceManager* faces = m_mesh.getFaceManager();
    faces->RegisterViewWrapper< array1d<double> >("Pressure");
  }

  const std::string& getPath() const
  { return m_file_path; }

  void setPath(const std::string& path)
  { m_file_path = path; }

  void write(double dt)
  {
    ++m_counter;
    FaceManager* faces = m_mesh.getFaceManager();
    NodeManager* nodes = m_mesh.getNodeManager();

    const OrderedVariableOneToManyRelation& face_connectivity = faces->nodeList();
    const localIndex n_faces = face_connectivity.size();

    /* Copy the face connectivity into a contiguous array. */
    std::int64_t* connectivity_array = new std::int64_t[4 * n_faces];
    for (localIndex i = 0; i < n_faces; ++i)
    {
      for (localIndex j = 0; j < 4; ++j)
      {
        connectivity_array[4 * i + j] = face_connectivity[i][j];
      }
    }

    integer const* ruptureState = faces->getReference<integer_array>("ruptureState").data();

    /* Build the face FieldMap. */
    FieldMap_in face_fields;
    real64* pressure_ptr = faces->getReference<real64_array>("Pressure").data();
    face_fields["Pressure"] = std::make_tuple(H5T_NATIVE_DOUBLE, 1, pressure_ptr);
    
    /* Update the dummy pressure field */
    for (localIndex i = 0; i < n_faces; ++i)
    {
      pressure_ptr[ i ] = (i + 1) * m_counter * 3.14159265359;
    }

    /* Build the node FieldMap. */
    r1_array const& reference_pos = nodes->getReference<r1_array>(nodes->viewKeys.referencePosition);
    localIndex const n_nodes = reference_pos.size();
    R1Tensor const* reference_pos_ptr = reference_pos.data();

    R1Tensor const* displacement_ptr = nodes->getReference<r1_array>(nodes->viewKeys.totalDisplacement).data();            

    FieldMap_in node_fields;
    node_fields["position"] = std::make_tuple(H5T_NATIVE_DOUBLE, 3, reference_pos_ptr);
    node_fields["displacement"] = std::make_tuple(H5T_NATIVE_DOUBLE, 3, displacement_ptr);

    writeBoundaryFile(m_comm, m_file_path.c_str(), m_face_offset,
                      m_n_faces_written, m_node_offset, 
                      m_n_nodes_written, dt, n_faces,
                      n_nodes, connectivity_array, ruptureState,
                      face_fields, node_fields);

    delete[] connectivity_array;
  }

private:
  MPI_Comm m_comm;
  std::string m_file_path;
  std::int64_t m_face_offset;
  std::int64_t m_n_faces_written;
  std::int64_t m_node_offset;
  std::int64_t m_n_nodes_written;
  MeshLevel& m_mesh;
  int m_counter;
};

} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_FILEIO_COUPLING_CHOMBOCOUPLER_HPP_ */
