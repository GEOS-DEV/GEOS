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
#include "codingUtilities/StringUtilities.hpp"

namespace geosx
{
namespace systemSolverInterface
{

enum class BlockIDs
{
  dummyScalarBlock,
  displacementBlock,
  fluidPressureBlock,
  temperatureBlock,
  compositionalBlock,
  invalidBlock
};

static double ClearRow ( Epetra_FECrsMatrix * matrix,
                         globalIndex const row,
                         const double factor )
{
  long long int rowTmp = static_cast<long long int>(row);
  int local_row = matrix->LRID(rowTmp);
  double LARGE = 0.0;

  if (local_row >= 0)
  {
    double *values = nullptr;
    int *col_indices = nullptr;
    int num_entries;

    matrix->ExtractMyRowView( local_row, num_entries, values, col_indices);


    if( values!=nullptr && col_indices!=nullptr && num_entries>0 )
    {
      int* diag_find = std::find(col_indices,col_indices+num_entries-1, local_row);
      long int diag_index = (diag_find - col_indices);

      for (int j=0 ; j<num_entries ; ++j)
      {
        if (diag_index != j )
        {
          values[j] = 0.;
        }
      }
      values[diag_index] *= factor;
      LARGE = values[diag_index];
    }
  }
  return (LARGE);
}

/**
 * @author settgast
 * @note class to hold the epetra system matrices and vectors.
 */
class EpetraBlockSystem
{
public:
  constexpr static int MAX_NUM_BLOCKS = 3;
  constexpr static int invalidIndex=-1;

  string BlockIDString( BlockIDs const id ) const
  {
    string rval;
    switch( id )
    {
      case BlockIDs::dummyScalarBlock:
        rval = "dummyScalarBlock";
        break;
      case BlockIDs::displacementBlock:
        rval = "displacementBlock";
        break;
    case BlockIDs::fluidPressureBlock:
      rval = "fluidPressureBlock";
      break;
    case BlockIDs::temperatureBlock:
      rval = "temperatureBlock";
      break;
    default:
      rval = "invalidBlock";
      break;
    }
    return rval;
  }


  EpetraBlockSystem();

  EpetraBlockSystem( EpetraBlockSystem const & ) = delete;
  EpetraBlockSystem( EpetraBlockSystem && ) = delete;
  EpetraBlockSystem& operator=( EpetraBlockSystem const & ) = delete;
  EpetraBlockSystem& operator=( EpetraBlockSystem && ) = delete;


  ~EpetraBlockSystem();



  void SetBlockID( const BlockIDs id, std::string const& name )
  {

    if( m_blockIndex.find(id) == m_blockIndex.end() )
    {
      m_blockIndex[id] = m_numBlocks;
    }
    else
    {
      GEOS_ERROR("error in EpetraBlockSystem::SetBlockID(). BlockIDs ("+BlockIDString(id)+") has already been registered to index "+
                 std::to_string(m_blockIndex[id]) );
    }

    if( m_blockID[m_numBlocks] == BlockIDs::invalidBlock )
    {
      m_blockID[m_numBlocks] = id;
    }
    else
    {
      GEOS_ERROR("error in EpetraBlockSystem::SetBlockID(). BlockIDs ("+BlockIDString(id)+") has already been registered to index "+
                 std::to_string(m_numBlocks) );
    }


    if( m_solverNames[m_numBlocks].empty() )
    {
      m_solverNames[m_numBlocks] = name;
    }
    else
    {
      GEOS_ERROR("error in EpetraBlockSystem::SetBlockID(). Block ("+std::to_string(
                   m_numBlocks)+") has already been used to register a solver named "+ m_solverNames[m_numBlocks] );
    }


    if( m_solverNameMap.find(name)==m_solverNameMap.end() )
    {
      m_solverNameMap[name] = m_numBlocks;
    }
    else
    {
      GEOS_ERROR("error in EpetraBlockSystem::SetBlockID(). Solver Name ("+name+") has already been used to register index "+ std::to_string(m_numBlocks) );
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


  Epetra_Map const * GetRowMap( const int index ) const
  {
    if( m_blockID[index]==BlockIDs::invalidBlock )
    {
      GEOS_ERROR("SolverBase.h:EpetraBlockSystem::GetRowMap():m_blockID isn't set \n");
    }
    return m_rowMap[ index ].get();
  }

  Epetra_Map * GetRowMap( const int index )
  {
    return const_cast<Epetra_Map *>( const_cast<EpetraBlockSystem const *>(this)->GetRowMap(index) );
  }


  Epetra_Map const * GetRowMap( const BlockIDs dofID ) const
  {
    int const index = m_blockIndex.find(dofID)->second;
    return GetRowMap(index);
  }

  Epetra_Map * GetRowMap( const BlockIDs dofID )
  {
    return const_cast<Epetra_Map *>( const_cast<EpetraBlockSystem const *>(this)->GetRowMap(dofID) );
  }



  Epetra_Map * SetRowMap( const BlockIDs dofID, std::unique_ptr<Epetra_Map> rowMap )
  {
    int index = m_blockIndex[dofID];
    m_rowMap[index] = std::move( rowMap );
    return m_rowMap[index].get();
  }


  Epetra_FEVector const * GetSolutionVector( int const index ) const
  {
    if( m_blockID[index]==BlockIDs::invalidBlock )
    {
      GEOS_ERROR("SolverBase.h:EpetraBlockSystem::GetSolutionVector():m_blockID isn't set \n");
    }
    return m_solution[ index ].get();
  }

  Epetra_FEVector * GetSolutionVector( int const index )
  {
    return const_cast<Epetra_FEVector *>( const_cast<EpetraBlockSystem const *>(this)->GetSolutionVector(index) );
  }


  Epetra_FEVector const * GetSolutionVector( const BlockIDs dofID ) const
  {
    int const index = m_blockIndex.find(dofID)->second;
    return GetSolutionVector(index);
  }

  Epetra_FEVector * GetSolutionVector( BlockIDs const dofID )
  {
    return const_cast<Epetra_FEVector *>( const_cast<EpetraBlockSystem const *>(this)->GetSolutionVector(dofID) );
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
      GEOS_ERROR("SolverBase.h:EpetraBlockSystem::GetResidualVector():m_blockID isn't set \n");
    }
    return m_rhs[ index ].get();
  }
  Epetra_FEVector const * GetResidualVector( int const index ) const
  {
    if( m_blockID[index]==BlockIDs::invalidBlock )
    {
      GEOS_ERROR("SolverBase.h:EpetraBlockSystem::GetResidualVector():m_blockID isn't set \n");
    }
    return m_rhs[ index ].get();
  }
  Epetra_FEVector * GetResidualVector( const BlockIDs dofID )
  {
    int index = m_blockIndex[dofID];
    return GetResidualVector(index);
  }
  Epetra_FEVector const * GetResidualVector( const BlockIDs dofID ) const
  {
    int index = m_blockIndex.at(dofID);
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
      GEOS_ERROR("SolverBase.h:EpetraBlockSystem::GetSparsity():m_blockID isn't set \n");
    }
    return m_sparsity[ rowIndex ][ colIndex ].get();
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
      GEOS_ERROR("SolverBase.h:EpetraBlockSystem::GetMatrix():m_blockID isn't set \n");
    }
    return m_matrix[ rowIndex ][ colIndex ].get();
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

  double ClearSystemRow ( int const blockRow,
                          globalIndex const rowIndex,
                          const double factor )
  {
    double LARGE = 0;
    for( int col=0 ; col<m_numBlocks ; ++col )
    {
      if( blockRow == col )
      {
        LARGE=ClearRow( m_matrix[blockRow][col].get(), rowIndex, factor );
      }
      else
      {
        ClearRow( m_matrix[blockRow][col].get(), rowIndex, 0.0 );
      }
    }
    return (LARGE);
  }

  double ClearSystemRow ( const BlockIDs rowDofID,
                          globalIndex const rowIndex,
                          const double factor )
  {
    int rowBlockIndex = m_blockIndex[rowDofID];
    return ClearSystemRow( rowBlockIndex, rowIndex, factor );
  }


  int numBlocks() const { return m_numBlocks; }


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

#endif /* SRC_COMPONENTS_CORE_SRC_LEGACY_PHYSICSSOLVERS_EPETRABLOCKSYSTEM_HPP_
        */
