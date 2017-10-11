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


  const static std::string numNodesToElemName[4];

  const NodeManager& m_node_manager;
  const ElementRegionManager& m_elem_reg_manager;
  const std::string m_output_path;
  const MPI_Comm m_comm;
  const std::string m_coord_name;
  const std::string m_topo_name;
};

} /* end namespace geosx */