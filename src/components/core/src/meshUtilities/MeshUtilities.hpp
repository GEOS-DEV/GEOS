/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

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
//                                                     ManagedGroup&
// faceManager,
//                                                     ManagedGroup&
// nodeManager);
//
//  static void GenerateElementsets ( xmlWrapper const & hdn,
//                                    const ManagedGroup& nodeManager,
//                                    ManagedGroup& elementManager);

};

}

#endif /* MESHUTILITIES_H_ */
