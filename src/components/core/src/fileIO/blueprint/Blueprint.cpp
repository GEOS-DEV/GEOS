#include "common/GeosxConfig.hpp"

#include "Blueprint.hpp"
#include "common/DataTypes.hpp"
#include "common/Logger.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/ElementRegionManager.hpp"

#ifdef USE_ATK
#include "sidre/sidre.hpp"
#include "sidre/DataStore.hpp"
#include "sidre/SidreTypes.hpp"
#include "spio/IOManager.hpp"

#include "conduit_blueprint.hpp"
#include "conduit_relay.hpp"
#endif

#include <cstring>
#include <unordered_map>
#include <utility>
#include <mpi.h>
#include <iostream>

using namespace axom::sidre;

namespace geosx {

const std::unordered_map<int, const std::string> Blueprint::numNodesToElemName = 
{ 
  {1, "point"},
  {2, "line"},
  {3, "triangle"},
  {4, "tet"},
  {8, "hex"}
};



Blueprint::Blueprint( const NodeManager& node_manager, const ElementRegionManager& elem_reg_manager, 
                      const std::string& output_path, MPI_Comm comm, const std::string& coord_name,
                      const std::string& topo_name):
  m_node_manager(node_manager),
  m_elem_reg_manager(elem_reg_manager),
  m_output_path(output_path),
  m_comm(comm),
  m_coord_name(coord_name),
  m_topo_name(topo_name)
{}



void Blueprint::write(int cycle) const
{
#ifdef USE_ATK
  const string mesh_name = "bp_mesh";
  
  DataStore ds;
  Group* root = ds.getRoot()->createGroup(mesh_name);
  Group* coords = root->createGroup("coordsets/" + m_coord_name);
  Group* topo = root->createGroup("topologies/" + m_topo_name);
  Group* fields = root->createGroup("fields");

  addNodes(coords, fields);
  addCells(topo, fields);

  conduit::Node mesh_node;
  conduit::Node info;
  ds.getRoot()->createNativeLayout(mesh_node);
  if (!conduit::blueprint::verify("mesh", mesh_node[mesh_name], info))
  {
    std::cout << "does not conform to the blueprint:";
    info.print();
    std::cout << std::endl;

    GEOS_ERROR("Does not conform to the blueprint. See above errors");
  }

  conduit::Node root_node;
  Node & index = root_node["blueprint_index"];

  conduit::blueprint::mesh::generate_index(mesh_node[mesh_name], mesh_name, 1, index[mesh_name]);

  info.reset();
  if (!conduit::blueprint::mesh::index::verify(index[mesh_name], info))
  {
    std::cout << "index does not conform to the blueprint:";
    info.print();
    std::cout << std::endl;

    GEOS_ERROR("Does not conform to the blueprint. See above errors");
  }


  const std::string root_output_path = m_output_path + "_" + std::to_string(cycle) + ".root";
  const std::string output_path = m_output_path + "_" + std::to_string(cycle) + ".hdf5";

  root_node["protocol/name"] = "conduit_hdf5";
  root_node["protocol/version"] = "0.1";
        
  root_node["number_of_files"] = 1;
  root_node["number_of_trees"] = 1;
  root_node["file_pattern"] = output_path;
  root_node["tree_pattern"] = "/";

  conduit::relay::io::save(root_node, root_output_path, "hdf5");
  conduit::relay::io::save(mesh_node, output_path);
#endif /* USE_ATK */
}


void Blueprint::addNodes(Group* coords, Group* fields) const
{
#ifdef USE_ATK
  coords->createView("type")->setString("explicit");

  localIndex n_nodes = m_node_manager.getWrapperBase(m_node_manager.viewKeys.referencePosition)->size();

  View * x_view = coords->createView("values/x", axom::sidre::TypeID::DOUBLE_ID, n_nodes);
  View * y_view = coords->createView("values/y", axom::sidre::TypeID::DOUBLE_ID, n_nodes);
  View * z_view = coords->createView("values/z", axom::sidre::TypeID::DOUBLE_ID, n_nodes);

  double * x = x_view->allocate()->getData();
  double * y = y_view->allocate()->getData();
  double * z = z_view->allocate()->getData();
  view_rtype_const<r1_array> position = m_node_manager.referencePosition();

  for (localIndex i = 0; i < n_nodes; i++)
  {
    x[i] = position[i][0];
    y[i] = position[i][1];
    z[i] = position[i][2];
  }

  for (const std::pair<const std::string, const ViewWrapperBase*>& pair : m_node_manager.wrappers())
  {
    const ViewWrapperBase *view = pair.second;
    if (view->sizedFromParent() == 1 &&
        view->size() > 0 &&
        view->getName() != m_node_manager.viewKeys.referencePosition.Key()) 
    {
      Group* curField = fields->createGroup(view->getName());
      curField->createView("association")->setString("vertex");
      curField->createView("topology")->setString(m_topo_name);
      curField->createView("volume_dependent")->setString("false");
      View* data = curField->createGroup("values")->createView("data");
      view->registerDataPtr(data);
    }
  }
#endif /* USE_ATK */
}


void Blueprint::addCells(Group* topo, Group* fields) const
{
#ifdef USE_ATK
  std::unordered_map<std::string, int> numberOfNodesOfType = 
  {
    { "point", 0 },
    { "line", 0 },
    { "tet", 0 },
    { "hex", 0 }
  };

  m_elem_reg_manager.forCellBlocks([this, &numberOfNodesOfType](const CellBlock* cell_block) 
  {
    const std::string& elem_name = this->numNodesToElemName.at(cell_block->numNodesPerElement());
    numberOfNodesOfType[elem_name] += cell_block->m_toNodesRelation.size();
  });

  std::unordered_map<std::string, std::pair<localIndex*, localIndex>> elementsToNodes;

  topo->createView("coordset")->setString(m_coord_name);
  topo->createView("type")->setString("unstructured");
  Group* elements = topo->createGroup("elements");

  for (std::pair<const std::string, int>& value : numberOfNodesOfType)
  {
    const std::string& elem_type = value.first;
    const localIndex num_nodes = value.second;
    if (num_nodes > 0)
    {
      Group* elem_group = elements->createGroup(elem_type);
      elem_group->createView("shape")->setString(elem_type);
      View* coordinates = elem_group->createView("connectivity", detail::SidreTT<localIndex>::id, num_nodes);
      coordinates->allocate();
      elementsToNodes[elem_type] = std::pair<localIndex*, localIndex>(coordinates->getArray(), 0);
    }
  }

  m_elem_reg_manager.forCellBlocks([this, &elementsToNodes](const CellBlock* cell_block) 
  {
    const std::string& elem_name = this->numNodesToElemName.at(cell_block->numNodesPerElement());
    std::pair<localIndex*, localIndex>& pair = elementsToNodes[elem_name];
    localIndex* to = pair.first + pair.second;
    const localIndex* from = cell_block->m_toNodesRelation[0];
    const localIndex num_nodes = cell_block->m_toNodesRelation.size();
    const localIndex byte_size = num_nodes * sizeof(localIndex);
    std::memcpy(to, from, byte_size);
    pair.second += num_nodes;
  });
#endif /* USE_ATK */
}


} /* end namespace geosx */
