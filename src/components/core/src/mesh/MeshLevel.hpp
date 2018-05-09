/*
 * MeshLevel.hpp
 *
 *  Created on: Sep 13, 2017
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MANAGERS_MESHLEVEL_HPP_
#define SRC_COMPONENTS_CORE_SRC_MANAGERS_MESHLEVEL_HPP_

#include "NodeManager.hpp"
#include "EdgeManager.hpp"
#include "FaceManager.hpp"
#include "ElementRegionManager.hpp"

namespace geosx
{
class ElementRegionManager;

class MeshLevel : public dataRepository::ManagedGroup
{
public:
  MeshLevel( string const & name,
             ManagedGroup * const parent );
  virtual ~MeshLevel() override;

  void InitializePostSubGroups( ManagedGroup * const ) override;

  void GenerateAdjacencyLists( localIndex_array & seedNodeList,
                               localIndex_array & nodeAdjacencyList,
                               localIndex_array & edgeAdjacencyList,
                               localIndex_array & faceAdjacencyList,
                               ElementRegionManager::ElementViewAccessor<localIndex_array> & elementAdjacencyList,
                               integer const depth );

  struct viewStructKeys
  {
    dataRepository::ViewKey meshLevel                = { "meshLevel" };
  } viewKeys;

  struct groupStructKeys
  {
    dataRepository::GroupKey vertexManager  = { "vertexManager" };
    dataRepository::GroupKey cellManager    = { "cellManager" };

    static constexpr auto nodeManagerString = "nodeManager";
    static constexpr auto edgeManagerString = "edgeManager";
    static constexpr auto faceManagerString = "faceManager";

    // This key is defined in problem manager:
    static constexpr auto elemManagerString = "ElementRegions";

    dataRepository::GroupKey nodeManager = {nodeManagerString};
    dataRepository::GroupKey edgeManager = {edgeManagerString};
    dataRepository::GroupKey faceManager = {faceManagerString};
    dataRepository::GroupKey elemManager = {elemManagerString};
  } groupKeys;

  NodeManager const * getNodeManager() const { return &m_nodeManager; }
  NodeManager * getNodeManager()             { return &m_nodeManager; }

  EdgeManager const * getEdgeManager() const { return this->GetGroup<EdgeManager>(groupKeys.edgeManager); }
  EdgeManager * getEdgeManager()             { return this->GetGroup<EdgeManager>(groupKeys.edgeManager); }

  FaceManager const * getFaceManager() const { return &m_faceManager; }
  FaceManager * getFaceManager()             { return &m_faceManager; }

  ElementRegionManager const * getElemManager() const { return &m_elementManager; }
  ElementRegionManager * getElemManager()             { return &m_elementManager; }

private:

  NodeManager m_nodeManager;
//  EdgeManager m_edgeManager;
  FaceManager m_faceManager;
  ElementRegionManager m_elementManager;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_MESHLEVEL_HPP_ */
