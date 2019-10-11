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

#ifndef GEOSX_FILEIO_BLUEPRINT_BLUEPRINT_HPP
#define GEOSX_FILEIO_BLUEPRINT_BLUEPRINT_HPP

#include "common/DataTypes.hpp"
#include "mpiCommunications/MpiWrapper.hpp"

#include <unordered_map>
#include <string>

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

  void write(int cycle, integer const eventCounter ) const;

private:
  void addNodes(axom::sidre::Group* coords, axom::sidre::Group* fields) const;

  void addCells(axom::sidre::Group* topo, axom::sidre::Group* fields) const;


  const static std::unordered_map<localIndex, const std::string> numNodesToElemName;

#ifdef GEOSX_USE_ATK
  const NodeManager& m_node_manager;
  const ElementRegionManager& m_elem_reg_manager;
//  const MPI_Comm m_comm;
#endif
  const std::string m_output_path;
  const std::string m_coord_name;
  const std::string m_topo_name;
};

} /* end namespace geosx */

#endif /* GEOSX_FILEIO_BLUEPRINT_BLUEPRINT_HPP */
