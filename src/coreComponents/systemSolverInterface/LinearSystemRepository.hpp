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

#include "common/DataTypes.hpp"
#include "codingUtilities/StringUtilities.hpp"
#include "dataRepository/ManagedGroup.hpp"

namespace geosx
{
namespace systemSolverInterface
{

//enum class BlockIDs
//{
//  dummyScalarBlock,
//  displacementBlock,
////  fluidPressureBlock,
//  temperatureBlock,
//  invalidBlock
//};

struct BlockIDs
{
  constexpr static auto dummyScalarBlock = "dummy";
  constexpr static auto fluidPressureBlock = "fluidPressure";
  constexpr static auto displacementBlock = "displacement";
  constexpr static auto invalidBlock = "invalid";
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
class LinearSystemRepository : public dataRepository::ManagedGroup
{
public:
  constexpr static int MAX_NUM_BLOCKS = 3;
  constexpr static int invalidIndex=-1;

//  string BlockIDString( BlockIDs const id ) const
//  {
//    string rval;
//    switch( id )
//    {
//      case BlockIDs::dummyScalarBlock:
//        rval = "dummyScalarBlock";
//        break;
//      case BlockIDs::displacementBlock:
//        rval = "displacementBlock";
//        break;
//    case BlockIDs::fluidPressureBlock:
//      rval = "fluidPressureBlock";
//      break;
//    case BlockIDs::temperatureBlock:
//      rval = "temperatureBlock";
//      break;
//    default:
//      rval = "invalidBlock";
//      break;
//    }
//    return rval;
//  }


  LinearSystemRepository();

  LinearSystemRepository( LinearSystemRepository const & ) = delete;
  LinearSystemRepository( LinearSystemRepository && ) = delete;
  LinearSystemRepository& operator=( LinearSystemRepository const & ) = delete;
  LinearSystemRepository& operator=( LinearSystemRepository && ) = delete;


  ~LinearSystemRepository();



//  void SetBlockID( const BlockIDs id, std::string const& name )
//  {
//
//    if( m_blockIndex.find(id) == m_blockIndex.end() )
//    {
//      m_blockIndex[id] = m_numBlocks;
//    }
//    else
//    {
//      GEOS_ERROR("error in EpetraBlockSystem::SetBlockID(). BlockIDs ("+BlockIDString(id)+") has already been registered to index "+
//                 std::to_string(m_blockIndex[id]) );
//    }
//
//    if( m_blockID[m_numBlocks] == BlockIDs::invalidBlock )
//    {
//      m_blockID[m_numBlocks] = id;
//    }
//    else
//    {
//      GEOS_ERROR("error in EpetraBlockSystem::SetBlockID(). BlockIDs ("+BlockIDString(id)+") has already been registered to index "+
//                 std::to_string(m_numBlocks) );
//    }
//
//
//    if( m_solverNames[m_numBlocks].empty() )
//    {
//      m_solverNames[m_numBlocks] = name;
//    }
//    else
//    {
//      GEOS_ERROR("error in EpetraBlockSystem::SetBlockID(). Block ("+std::to_string(
//                   m_numBlocks)+") has already been used to register a solver named "+ m_solverNames[m_numBlocks] );
//    }
//
//
//    if( m_solverNameMap.find(name)==m_solverNameMap.end() )
//    {
//      m_solverNameMap[name] = m_numBlocks;
//    }
//    else
//    {
//      GEOS_ERROR("error in EpetraBlockSystem::SetBlockID(). Solver Name ("+name+") has already been used to register index "+ std::to_string(m_numBlocks) );
//    }
//
//
//    ++m_numBlocks;
//
//  }


  template< typename T >
  T const * GetRowMap( string const & name ) const
  {
    string fullName = name +"_rowMap";
    return &(this->getReference<T>( fullName ));
  }

  template< typename T >
  T * GetRowMap( string const & name )
  {
    string fullName = name +"_rowMap";
    return &(this->getReference<T>( fullName ));
  }

  template< typename T >
  T * SetRowMap( string const & name, std::unique_ptr<T> rowMap )
  {
    string fullName = name +"_rowMap";
    return RegisterViewWrapper( fullName, std::move(rowMap) )->
           setRestartFlag(dataRepository::RestartFlags::NO_WRITE)->getPointer();
  }



  template< typename T >
  T const * GetSolutionVector( string const & name ) const
  {
    string fullName = name +"_solution";
    return &(this->getReference<T>( fullName ));
  }

  template< typename T >
  T * GetSolutionVector( string const & name )
  {
    string fullName = name +"_solution";
    return &(this->getReference<T>( fullName ));
  }

  template< typename T >
  T * SetSolutionVector( string const & name, std::unique_ptr<T> solution )
  {
    string fullName = name +"_solution";
    return RegisterViewWrapper( fullName, std::move(solution) )->
           setRestartFlag(dataRepository::RestartFlags::NO_WRITE)->getPointer();
  }

  template< typename T >
  T * GetResidualVector( string const & name )
  {
    string fullName = name +"_residual";
    return &(this->getReference<T>( fullName ));
  }
  template< typename T >
  T const * GetResidualVector( string const & name ) const
  {
    string fullName = name +"_residual";
    return &(this->getReference<T>( fullName ));
  }

  template< typename T >
  T * SetResidualVector( string const & name, std::unique_ptr<T> rhs )
  {
    string fullName = name +"_residual";
    return RegisterViewWrapper( fullName, std::move(rhs) )->
           setRestartFlag(dataRepository::RestartFlags::NO_WRITE)->getPointer();
  }


  template< typename T >
  T * GetSparsity( string const & rowDofID,
                   string const & colDofID )
  {
    string fullName = rowDofID + "-" + colDofID +"_sparsity";
    return &(this->getReference<T>( fullName ));
  }

  template< typename T >
  T * SetSparsity( string const & rowDofID,
                                   string const & colDofID,
                                   std::unique_ptr<Epetra_FECrsGraph> sparsity )
  {
    string fullName = rowDofID + "-" + colDofID +"_sparsity";
    return RegisterViewWrapper( fullName, std::move(sparsity) )->
           setRestartFlag(dataRepository::RestartFlags::NO_WRITE)->getPointer();
  }




  template< typename T >
  T * GetMatrix( string const & rowDofID,
                 string const & colDofID )
  {
    string fullName = rowDofID + "-" + colDofID +"_matrix";
    return &(this->getReference<T>( fullName ));
  }

  template< typename T >
  T * SetMatrix( string const & rowDofID,
                 string const & colDofID,
                 std::unique_ptr<T> matrix )
  {
    string fullName = rowDofID + "-" + colDofID +"_matrix";
    return RegisterViewWrapper( fullName, std::move(matrix) )->
           setRestartFlag(dataRepository::RestartFlags::NO_WRITE)->getPointer();

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
        if( m_matrix[blockRow][col].get() != nullptr )
        {
          ClearRow( m_matrix[blockRow][col].get(), rowIndex, 0.0 );
        }
      }
    }
    return (LARGE);
  }

  double ClearSystemRow ( const BlockIDs rowDofID,
                          globalIndex const rowIndex,
                          const double factor )
  {
//    int rowBlockIndex = m_blockIndex[rowDofID];
    return 0;//ClearSystemRow( rowBlockIndex, rowIndex, factor );
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


  LinearSystemRepository( LinearSystemRepository& );
  LinearSystemRepository& operator=(LinearSystemRepository&);

};

}
} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_LEGACY_PHYSICSSOLVERS_EPETRABLOCKSYSTEM_HPP_
        */
