#include "Blueprint.hpp"
#include "managers/NodeManager.hpp"
#include "managers/ElementRegionManager.hpp"

#include "sidre/sidre.hpp"
#include "sidre/DataStore.hpp"
#include "sidre/SidreTypes.hpp"
#include "spio/IOManager.hpp"

#include <cstring>
#include <unordered_map>
#include <utility>
#include <mpi.h>

using namespace axom::sidre;

namespace geosx {

const std::string Blueprint::numNodesToElemName[] = 
{ 
  "point",
  "line",
  "tet",
  "hex"
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
  DataStore ds;
  Group* root = ds.getRoot();

  Group* coords = root->createGroup("coordsets/" + m_coord_name);
  Group* topo = root->createGroup("topologies/" + m_topo_name);
  Group* fields = root->createGroup("fields");
  addNodes(coords, fields);
  addCells(topo, fields);

  axom::spio::IOManager ioManager(m_comm);
  ioManager.write(root, 1, m_output_path + "_" + std::to_string(cycle), "sidre_hdf5");
}


void Blueprint::addNodes(Group* coords, Group* fields) const
{
  coords->createView("type")->setString("explicit");
  View* xyz = coords->createGroup("values")->createView("xyz");
  m_node_manager.getWrapperBase(m_node_manager.viewKeys.referencePosition)->registerDataPtr(xyz);

  for (const std::pair<const std::string, const ViewWrapperBase*>& pair : m_node_manager.wrappers())
  {
    const ViewWrapperBase *view = pair.second;
    if (view->sizedFromParent() == 1 && view->getName() != m_node_manager.viewKeys.referencePosition.Key()) 
    {
      Group* curField = fields->createGroup(view->getName());
      curField->createView("association")->setString("vertex");
      curField->createView("topology")->setString(m_topo_name);
      curField->createView("volume_dependent")->setString("false");
      View* data = curField->createGroup("values")->createView("data");
      view->registerDataPtr(data);
    }
  }
}


void Blueprint::addCells(Group* topo, Group* fields) const
{
  std::unordered_map<std::string, int> numberOfNodesOfType = 
  {
    { "point", 0 },
    { "line", 0 },
    { "tet", 0 },
    { "hex", 0 }
  };

  m_elem_reg_manager.forCellBlocks([this, &numberOfNodesOfType](const CellBlock* cell_block) 
  {
    const std::string& elem_name = this->numNodesToElemName[cell_block->numNodesPerElement()];
    numberOfNodesOfType[elem_name] += cell_block->m_toNodesRelation.size();
  });

  std::unordered_map<std::string, std::pair<int32*, int32>> elementsToNodes;

  for (std::pair<const std::string, int>& value : numberOfNodesOfType)
  {
    const std::string& elem_type = value.first;
    const int32 num_nodes = value.second;
    if (num_nodes > 0)
    {
      Group* elem_group = topo->createGroup(elem_type);
      elem_group->createView("coordset")->setString(m_coord_name);
      elem_group->createView("type")->setString("unstructured");
      elem_group->createView("elements/shape")->setString(elem_type);
      View* coordinates = elem_group->createView("elements/connectivity", INT32_ID, num_nodes);
      coordinates->allocate();
      elementsToNodes[elem_type] = std::pair<int32*, int32>(coordinates->getArray(), 0);
    }
  }

  m_elem_reg_manager.forCellBlocks([this, &elementsToNodes](const CellBlock* cell_block) 
  {
    const std::string& elem_name = this->numNodesToElemName[cell_block->numNodesPerElement()];
    std::pair<int32*, int32>& pair = elementsToNodes[elem_name];
    int32* to = pair.first + pair.second;
    const int32* from = cell_block->m_toNodesRelation[0];
    const int32 num_nodes = cell_block->m_toNodesRelation.size();
    const std_size_t byte_size = num_nodes * sizeof(int32);
    std::memcpy(to, from, byte_size);
    pair.second += num_nodes;
  });
}


} /* end namespace geosx */ 