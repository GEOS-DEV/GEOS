/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

#ifndef FILEIO_BLUEPRINT_BLUEPRINT_HPP
#define FILEIO_BLUEPRINT_BLUEPRINT_HPP

#include "common/DataTypes.hpp"

#include <unordered_map>
#include <string>
#include <mpi.h>

/* Forward declarations */
namespace axom
{
namespace sidre
{
class Group;
}
}

namespace geosx
{


/* Forward declarations */
class NodeManager;
class ElementRegionManager;

class Blueprint
{
public:
  Blueprint( const NodeManager& node_manager, const ElementRegionManager& elem_reg_manager,
             const std::string& output_path, MPI_Comm comm,
             const std::string& coord_name="coords", const std::string& topo_name="mesh");

  ~Blueprint()
  {}

  void write(int cycle) const;

private:
  void addNodes(axom::sidre::Group* coords, axom::sidre::Group* fields) const;

  void addCells(axom::sidre::Group* topo, axom::sidre::Group* fields) const;


  const static std::unordered_map<localIndex, const std::string> numNodesToElemName;

#ifdef USE_ATK
  const NodeManager& m_node_manager;
  const ElementRegionManager& m_elem_reg_manager;
  const MPI_Comm m_comm;
#endif
  const std::string m_output_path;
  const std::string m_coord_name;
  const std::string m_topo_name;
};

} /* end namespace geosx */

#endif /* FILEIO_BLUEPRINT_BLUEPRINT_HPP */
