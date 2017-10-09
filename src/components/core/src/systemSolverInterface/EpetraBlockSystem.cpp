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
