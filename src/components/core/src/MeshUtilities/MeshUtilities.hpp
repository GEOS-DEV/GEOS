/*
 * MeshUtilities.h
 *
 *  Created on: Dec 5, 2012
 *      Author: settgast1
 */

#ifndef MESHUTILITIES_H_
#define MESHUTILITIES_H_

#include "common/DataTypes.hpp"

namespace geosx
{

namespace dataRepository
{
class ManagedGroup;
}
class xmlWrapper;

class MeshUtilities
{
public:
  MeshUtilities();
  virtual ~MeshUtilities();




  static void GenerateNodesets( dataRepository::ManagedGroup const * geometry,
                                dataRepository::ManagedGroup * nodeManager );

//  static void GenerateFasesetsAndAssociatedNodesets( xmlWrapper const & hdn,
//                                                     ManagedGroup& faceManager,
//                                                     ManagedGroup& nodeManager);
//
//  static void GenerateElementsets ( xmlWrapper const & hdn,
//                                    const ManagedGroup& nodeManager,
//                                    ManagedGroup& elementManager);

};

}

#endif /* MESHUTILITIES_H_ */
