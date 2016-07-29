/*
 * MeshUtilities.h
 *
 *  Created on: Dec 5, 2012
 *      Author: settgast1
 */

#ifndef MESHUTILITIES_H_
#define MESHUTILITIES_H_

#include "Common/Common.h"

class NodeManagerT;
class FaceManagerT;
class ElementManagerT;
namespace TICPP
{
class HierarchicalDataNode;
}

class MeshUtilities
{
public:
  MeshUtilities();
  virtual ~MeshUtilities();




  static void GenerateNodesets( TICPP::HierarchicalDataNode& hdn,
                                NodeManagerT& nodeManager );

  static void GenerateFasesetsAndAssociatedNodesets( TICPP::HierarchicalDataNode& hdn,
                                                     FaceManagerT& faceManager,
                                                     NodeManagerT& nodeManager);

  static void GenerateElementsets ( TICPP::HierarchicalDataNode& hdn,
                                    const NodeManagerT& nodeManager,
                                    ElementManagerT& elementManager);

};

#endif /* MESHUTILITIES_H_ */
