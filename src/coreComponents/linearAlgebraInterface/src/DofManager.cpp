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

/**
 * @file DofManager.cpp
 */

#include "DofManager.hpp"

namespace geosx
{

// Constructor
DofManager::DofManager(MeshLevel const & meshLevel)
  :
  m_meshLevel(meshLevel)
{}


DofManager::close()
{
  NodeManager const * nodeManager = meshLevel.getNodeManager();
  localIndex n_ghost_rows  = nodeManager->GetNumberOfGhosts();
  localIndex n_local_rows  = nodeManager->size()-n_ghost_rows;
  std::cout << "No local rows" << n_local_rows << std::endl; 
}


}

