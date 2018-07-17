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
 * EpetraBlockSystem.cpp
 *
 *  Created on: Sep 19, 2017
 *      Author: settgast
 */

#include "EpetraBlockSystem.hpp"

namespace geosx
{
namespace systemSolverInterface
{

EpetraBlockSystem::EpetraBlockSystem():
  m_blockID(),
  m_blockIndex(),
  m_solverNames(),
  m_solverNameMap(),
  m_numBlocks(0),
  m_rowMap(),
  m_solution(),
  m_lastSolution(),
  m_rhs(),
  m_sparsity(),
  m_matrix()
{
  for( int i=0 ; i<MAX_NUM_BLOCKS ; ++i )
  {
    m_blockID[i] = BlockIDs::invalidBlock;
  }
}

EpetraBlockSystem::~EpetraBlockSystem()
{
  // TODO Auto-generated destructor stub
}

}
} /* namespace geosx */
