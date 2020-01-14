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

/**
 * @file MeshUtilities.h
 *
 */

#ifndef MESHUTILITIES_H_
#define MESHUTILITIES_H_

#include "common/DataTypes.hpp"

namespace geosx
{

namespace dataRepository
{
class Group;
}

class ObjectManagerBase;
class xmlWrapper;
class NodeManager;

class MeshUtilities
{
public:
  MeshUtilities();
  virtual ~MeshUtilities();



  static void GenerateNodesets( dataRepository::Group const * geometry,
                                NodeManager * const nodeManager );

//  static void GenerateFasesetsAndAssociatedNodesets( xmlWrapper const & hdn,
//                                                     Group&
// faceManager,
//                                                     Group&
// nodeManager);
//
//  static void GenerateElementsets ( xmlWrapper const & hdn,
//                                    const Group& nodeManager,
//                                    Group& elementManager);

};

}

#endif /* MESHUTILITIES_H_ */
