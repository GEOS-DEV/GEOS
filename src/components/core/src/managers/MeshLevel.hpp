/*
 * MeshLevel.hpp
 *
 *  Created on: Sep 13, 2017
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MANAGERS_MESHLEVEL_HPP_
#define SRC_COMPONENTS_CORE_SRC_MANAGERS_MESHLEVEL_HPP_

#include "dataRepository/ManagedGroup.hpp"

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

    dataRepository::GroupKey nodeManager    = { "nodeManager" };
    dataRepository::GroupKey edgeManager    = { "edgeManager" };
    dataRepository::GroupKey faceManager    = { "faceManager" };
    dataRepository::GroupKey elementManager = { "elementManager" };
  }groupKeys;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_MESHLEVEL_HPP_ */
