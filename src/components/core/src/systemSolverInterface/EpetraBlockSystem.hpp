/*
 * EpetraBlockSystem.hpp
 *
 *  Created on: Sep 19, 2017
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_LEGACY_PHYSICSSOLVERS_EPETRABLOCKSYSTEM_HPP_
#define SRC_COMPONENTS_CORE_SRC_LEGACY_PHYSICSSOLVERS_EPETRABLOCKSYSTEM_HPP_

#include "Epetra_Map.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_LinearProblem.h"

#include "Teuchos_RCP.hpp"
#include "AztecOO.h"
#include "ml_MultiLevelPreconditioner.h"
#include "EpetraExt_MatrixMatrix.h"

#include "common/Logger.hpp"

namespace geosx
{
namespace systemSolverInterface
{


/**
 * @author settgast
 * @note class to hold the epetra system matrices and vectors.
 */
class EpetraBlockSystem
{
public:
  constexpr static int MAX_NUM_BLOCKS = 3;
  constexpr static int invalidIndex=-1;
  enum class BlockIDs
  {
    displacementBlock,
    fluidPressureBlock,
    temperatureBlock,
    invalidBlock
  };

  EpetraBlockSystem();

  ~EpetraBlockSystem();



  void SetBlockID( const BlockIDs id, std::string const& name )
  {

    if( m_blockIndex.find(id) == m_blockIndex.end() )
    {
      m_blockIndex[id] = m_numBlocks;
    }
    else
    {
      throw GEOS_ERROR("error in EpetraBlockSystem::SetBlockID(). BlockIDs ("+toString(id)+") has already been registered to index "+ toString(index) );
    }

    if( m_blockID[m_numBlocks] == BlockIDs::invalidBlock )
    {
      m_blockID[m_numBlocks] = id;
    }
    else
    {
      throw GEOS_ERROR("error in EpetraBlockSystem::SetBlockID(). BlockIDs ("+toString(id)+") has already been registered to index "+ toString(index) );
    }


    if( m_solverNames[m_numBlocks].empty() )
    {
      m_solverNames[m_numBlocks] = name;
    }
    else
    {
      throw GEOS_ERROR("error in EpetraBlockSystem::SetBlockID(). Block ("+toString(block)+") has already been used to register a solver named "+ m_solverNames[block] );
    }


    if( m_solverNameMap.find(name)==m_solverNameMap.end() )
    {
      m_solverNameMap[name] = m_numBlocks;
    }
    else
    {
      throw GEOS_ERROR("error in EpetraBlockSystem::SetBlockID(). Solver Name ("+name+") has already been used to register index "+ toString(index)+ ". Ya betta check yoself befo ya wreck yoself.\n" );
    }


    ++m_numBlocks;

  }

  BlockIDs GetBlockID( int index ) const
  {
    return m_blockID[index];
  }

  int GetIndex( BlockIDs id ) const
  {
    return this->m_blockIndex.at(id);
  }

  std::string GetSolverName( int index ) const
  {
    return m_solverNames[index];
  }


  Epetra_Map * GetRowMap( const int index )
  {
    if( m_blockID[index]==BlockIDs::invalidBlock )
    {
      throw GEOS_ERROR("SolverBase.h:EpetraBlockSystem::GetRowMap():m_blockID isn't set \n");
    }
    return m_rowMap[ index ].get();
  }

  Epetra_Map * GetRowMap( const BlockIDs dofID )
  {
    int index = m_blockIndex[dofID];
    return GetRowMap(index);
  }


  Epetra_Map * SetRowMap( const BlockIDs dofID, std::unique_ptr<Epetra_Map> rowMap )
  {
    int index = m_blockIndex[dofID];
    m_rowMap[index] = std::move( rowMap );
    return m_rowMap[index].get();
  }


  Epetra_FEVector * GetSolutionVector( int const index )
  {
    if( m_blockID[index]==BlockIDs::invalidBlock )
    {
      throw GEOS_ERROR("SolverBase.h:EpetraBlockSystem::GetSolutionVector():m_blockID isn't set \n");
    }
    return m_solution[ index ].get();
  }

  Epetra_FEVector * GetSolutionVector( const BlockIDs dofID )
  {
    int index = m_blockIndex[dofID];
    return GetSolutionVector(index);
  }

  Epetra_FEVector * SetSolutionVector( const BlockIDs dofID, std::unique_ptr<Epetra_FEVector> solution )
  {
    int index = m_blockIndex[dofID];
    m_solution[index] = std::move( solution );
    return m_solution[index].get();
  }



  Epetra_FEVector * GetResidualVector( int const index )
  {
    if( m_blockID[index]==BlockIDs::invalidBlock )
    {
      throw GEOS_ERROR("SolverBase.h:EpetraBlockSystem::GetResidualVector():m_blockID isn't set \n");
    }
    return m_rhs[ index ].get();
  }
  Epetra_FEVector * GetResidualVector( const BlockIDs dofID )
  {
    int index = m_blockIndex[dofID];
    return GetResidualVector(index);
  }

  Epetra_FEVector * SetResidualVector( const BlockIDs dofID, std::unique_ptr<Epetra_FEVector> rhs )
  {
    int index = m_blockIndex[dofID];
    m_rhs[index] = std::move( rhs );
    return m_rhs[index].get();
  }




  Epetra_FECrsGraph * GetSparsity( int const rowIndex,
                                   int const colIndex )
  {
    if( m_blockID[rowIndex]==BlockIDs::invalidBlock && m_blockID[colIndex]==BlockIDs::invalidBlock )
    {
      throw GEOS_ERROR("SolverBase.h:EpetraBlockSystem::GetSparsity():m_blockID isn't set \n");
    }
    return m_sparsity[ rowIndex ][ colIndex ].get() ;
  }
  Epetra_FECrsGraph * GetSparsity( const BlockIDs rowDofID,
                                   const BlockIDs colDofID )
  {
    int rowIndex = m_blockIndex[rowDofID];
    int colIndex = m_blockIndex[colDofID];
    return GetSparsity( rowIndex, colIndex );
  }

  Epetra_FECrsGraph * SetSparsity( const BlockIDs rowDofID,
                                   const BlockIDs colDofID,
                                   std::unique_ptr<Epetra_FECrsGraph> sparsity )
  {
    int rowIndex = m_blockIndex[rowDofID];
    int colIndex = m_blockIndex[colDofID];

    m_sparsity[rowIndex][colIndex] = std::move( sparsity );
    return m_sparsity[rowIndex][colIndex].get();
  }



  Epetra_FECrsMatrix * GetMatrix( int const rowIndex,
                                  int const colIndex )
  {
    if( m_blockID[rowIndex]==BlockIDs::invalidBlock && m_blockID[colIndex]==BlockIDs::invalidBlock )
    {
      throw GEOS_ERROR("SolverBase.h:EpetraBlockSystem::GetMatrix():m_blockID isn't set \n");
    }
    return m_matrix[ rowIndex ][ colIndex ].get() ;
  }

  Epetra_FECrsMatrix * GetMatrix( const BlockIDs rowDofID,
                                  const BlockIDs colDofID )
  {
    int rowIndex = m_blockIndex[rowDofID];
    int colIndex = m_blockIndex[colDofID];
    return GetMatrix( rowIndex, colIndex );
  }

  Epetra_FECrsMatrix * SetMatrix( const BlockIDs rowDofID,
                                   const BlockIDs colDofID,
                                   std::unique_ptr<Epetra_FECrsMatrix> matrix )
  {
    int rowIndex = m_blockIndex[rowDofID];
    int colIndex = m_blockIndex[colDofID];

    m_matrix[rowIndex][colIndex] = std::move( matrix );
    return m_matrix[rowIndex][colIndex].get();
  }




private:
  /// m_blockID
  BlockIDs m_blockID[MAX_NUM_BLOCKS];
  std::map<BlockIDs, int> m_blockIndex;
  std::string m_solverNames[MAX_NUM_BLOCKS];
  std::map<std::string, int> m_solverNameMap;
  int m_numBlocks = 0;

  std::unique_ptr<Epetra_Map>         m_rowMap[MAX_NUM_BLOCKS];
  std::unique_ptr<Epetra_FEVector>    m_solution[MAX_NUM_BLOCKS];
  std::unique_ptr<Epetra_FEVector>    m_lastSolution[MAX_NUM_BLOCKS];
  std::unique_ptr<Epetra_FEVector>    m_rhs[MAX_NUM_BLOCKS];
  std::unique_ptr<Epetra_FECrsGraph>  m_sparsity[MAX_NUM_BLOCKS][MAX_NUM_BLOCKS];
  std::unique_ptr<Epetra_FECrsMatrix> m_matrix[MAX_NUM_BLOCKS][MAX_NUM_BLOCKS];


  EpetraBlockSystem( EpetraBlockSystem& );
  EpetraBlockSystem& operator=(EpetraBlockSystem&);

};

}
} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_LEGACY_PHYSICSSOLVERS_EPETRABLOCKSYSTEM_HPP_ */
