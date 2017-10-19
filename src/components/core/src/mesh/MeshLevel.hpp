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

class MeshLevel : public dataRepository::ManagedGroup
{
public:
  MeshLevel( string const & name,
             ManagedGroup * const parent );
  virtual ~MeshLevel();

  struct viewStructKeys
  {
    dataRepository::ViewKey meshLevel                = { "meshLevel" };
  }viewKeys;

  struct groupStructKeys
  {
    dataRepository::GroupKey vertexManager  = { "vertexManager" };
    dataRepository::GroupKey cellManager    = { "cellManager" };

    static constexpr auto nodeManagerString = "nodeManager";
    static constexpr auto edgeManagerString = "edgeManager";
    static constexpr auto faceManagerString = "faceManager";
    static constexpr auto elemManagerString = "elementManager";

    dataRepository::GroupKey nodeManager = {nodeManagerString};
    dataRepository::GroupKey edgeManager = {edgeManagerString};
    dataRepository::GroupKey faceManager = {faceManagerString};
    dataRepository::GroupKey elemManager = {elemManagerString};
  }groupKeys;

  NodeManager const * getNodeManager() const { return this->GetGroup<NodeManager>(groupKeys.nodeManager); }
  NodeManager * getNodeManager()             { return this->GetGroup<NodeManager>(groupKeys.nodeManager); }

  EdgeManager const * getEdgeManager() const { return this->GetGroup<EdgeManager>(groupKeys.edgeManager); }
  EdgeManager * getEdgeManager()             { return this->GetGroup<EdgeManager>(groupKeys.edgeManager); }

  FaceManager const * getFaceManager() const { return this->GetGroup<FaceManager>(groupKeys.faceManager); }
  FaceManager * getFaceManager()             { return this->GetGroup<FaceManager>(groupKeys.faceManager); }

  ElementRegionManager const * getElemManager() const { return this->GetGroup<ElementRegionManager>(groupKeys.elemManager); }
  ElementRegionManager * getElemManager()             { return this->GetGroup<ElementRegionManager>(groupKeys.elemManager); }

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_MESHLEVEL_HPP_ */
