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

namespace geosx
{
namespace systemSolverInterface
{


/**
 * @author settgast
 * @note class to hold the epetra system matrices and vectors.
 */
template< int T_DIM >
class EpetraBlockSystem
{
public:
  enum class BlockIDs
  {
    displacementBlock,
    fluidPressureBlock,
    temperatureBlock,
    numBlocks,
    invalidBlock
  };

    EpetraBlockSystem():
      m_blockID(),
      m_solverNames(),
      m_solverNameMap()
    {
      for( int i=0 ; i<T_DIM ; ++i )
      {
        m_blockID[i] = BlockIDs::invalidBlock;
      }
    }

    ~EpetraBlockSystem()
    {}

    void SetBlockID( int const block, const BlockIDs id, std::string const& name )
    {
      if( block>=T_DIM )
      {
        throw GEOS_ERROR("SolverBase.h:EpetraBlockSystem::SetBlockID(): block index is larger than number of blocks\n");
      }
      else if( block<0 )
      {
        throw GEOS_ERROR("SolverBase.h:EpetraBlockSystem::SetBlockID(): block index is less than 0\n");
      }

      if( m_blockID[block]==BlockIDs::invalidBlock )
      {
        m_blockID[block] = id;
      }
      else
      {
        throw GEOS_ERROR("error in EpetraBlockSystem::SetBlockID(). BlockIDs ("+toString(id)+") has already been registered to block "+ toString(block) );
      }



      if( m_solverNames[block].empty() )
      {
        m_solverNames[block] = name;
      }
      else
      {
        throw GEOS_ERROR("error in EpetraBlockSystem::SetBlockID(). Block ("+toString(block)+") has already been used to register a solver named "+ m_solverNames[block] );
      }



      std::map<std::string, int>::iterator iterSolverNameMap = m_solverNameMap.find(name);
      if( iterSolverNameMap==m_solverNameMap.end() )
      {
        m_solverNameMap[name] = block;
      }
      else
      {
        throw GEOS_ERROR("error in EpetraBlockSystem::SetBlockID(). Solver Name ("+name+") has already been used to register block "+ toString(block)+ ". Ya betta check yoself befo ya wreck yoself.\n" );
      }
    }

    int GetBlockID( BlockIDs id ) const
    {
      return m_blockID[id];
    }

    std::string GetSolverName( BlockIDs id ) const
    {
      return m_solverNames[id];
    }


    std::unique_ptr<Epetra_Map> & GetRowMap( const BlockIDs dofID )
    {
      if( m_blockID[dofID]==BlockIDs::invalidBlock )
      {
        throw GEOS_ERROR("SolverBase.h:EpetraBlockSystem::GetRowMap():m_blockID isn't set \n");
      }
      return m_rowMap[ m_blockID[dofID] ];
    }

    std::unique_ptr<Epetra_FEVector> & GetSolutionVector( const BlockIDs dofID )
    {
      if( m_blockID[dofID]==BlockIDs::invalidBlock )
      {
        throw GEOS_ERROR("SolverBase.h:EpetraBlockSystem::GetSolutionVector():m_blockID isn't set \n");
      }
      return m_solution[ m_blockID[dofID] ];
    }

    std::unique_ptr<Epetra_FEVector> & GetResidualVector( const BlockIDs dofID )
    {
      if( m_blockID[dofID]==BlockIDs::invalidBlock )
      {
        throw GEOS_ERROR("SolverBase.h:EpetraBlockSystem::GetResidualVector():m_blockID isn't set \n");
      }
      return m_rhs[ m_blockID[dofID] ];
    }

    std::unique_ptr<Epetra_FECrsGraph> & GetSparsity( const BlockIDs rowDofID,
                                                      const BlockIDs colDofID )
    {
      if( m_blockID[rowDofID]==BlockIDs::invalidBlock && m_blockID[colDofID]==BlockIDs::invalidBlock )
      {
        throw GEOS_ERROR("SolverBase.h:EpetraBlockSystem::GetSparsity():m_blockID isn't set \n");
      }
      return m_sparsity[ m_blockID[rowDofID] ][ m_blockID[colDofID] ] ;
    }

    std::unique_ptr<Epetra_FECrsMatrix> & GetMatrix( const BlockIDs rowDofID,
                                                     const BlockIDs colDofID )
    {
      if( m_blockID[rowDofID]==BlockIDs::invalidBlock && m_blockID[colDofID]==BlockIDs::invalidBlock )
      {
        throw GEOS_ERROR("SolverBase.h:EpetraBlockSystem::GetMatrix():m_blockID isn't set \n");
      }
      return m_matrix[ m_blockID[rowDofID] ][ m_blockID[colDofID] ] ;
    }

    bool HasMatrixBlock( int const rowDofID,
                         int const colDofID )
    {
      return ( m_blockID[rowDofID]!=BlockIDs::invalidBlock && m_blockID[colDofID]!=BlockIDs::invalidBlock ) ? true : false;
    }

    std::unique_ptr<Epetra_Map>         m_rowMap[T_DIM];
    std::unique_ptr<Epetra_FEVector>    m_solution[T_DIM];
    std::unique_ptr<Epetra_FEVector>    m_lastSolution[T_DIM];
    std::unique_ptr<Epetra_FEVector>    m_rhs[T_DIM];
    std::unique_ptr<Epetra_FECrsGraph>  m_sparsity[T_DIM][T_DIM];
    std::unique_ptr<Epetra_FECrsMatrix> m_matrix[T_DIM][T_DIM];


  private:
    /// m_blockID
    BlockIDs m_blockID[T_DIM];
    std::string m_solverNames[T_DIM];
    std::map<std::string, int> m_solverNameMap;

    EpetraBlockSystem( EpetraBlockSystem& );
    EpetraBlockSystem& operator=(EpetraBlockSystem&);

  };

}

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_LEGACY_PHYSICSSOLVERS_EPETRABLOCKSYSTEM_HPP_ */
