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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file SolverBase.h
 * @author settgast1
 * @date Feb 10, 2011
 */

#ifndef SOLVERBASE_H_
#define SOLVERBASE_H_


#include "Common/typedefs.h"
#include "IO/ticpp/TinyXMLParser.h"
#include "MPI_Communications/SpatialPartition.h"
#include "ObjectManagers/PhysicalDomainT.h"
//#include "ObjectManagers/ProblemManagerT.h"

#if GPAC_MPI
class Epetra_MpiComm;
#else
class Epetra_SerialComm;
#endif

#include "Epetra_Map.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_LinearProblem.h"

#include "EpetraExt_SolverMap_CrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

#include "AztecOO.h"

#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"
#include "ml_RowMatrix.h"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

#include "Thyra_OperatorVectorClientSupport.hpp"
#include "Thyra_AztecOOLinearOpWithSolveFactory.hpp"
#include "Thyra_AztecOOLinearOpWithSolve.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraLinearOpBase.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_DefaultInverseLinearOp.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_MLPreconditionerFactory.hpp"


#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace ML_Epetra
{ class MultiLevelPreconditioner; }

class AztecOO;

class PhysicalDomainT;
class ProblemManagerT;
class FractunatorBase;

namespace EpetraBlock
{
enum ID
{
  solidBlock,
  fluidBlock,
  numBlockDof
};
}

/**
 * @author settgast
 * @note class to hold the epetra system matrices and vectors.
 */
class Epetra_System
{

public:

  Epetra_System():
    m_nblocks(1),
    m_blockID(),
    m_solverNames(),
    m_solverNameMap()
  {
    for( int i=0 ; i<EpetraBlock::numBlockDof ; ++i )
    {
      m_blockID[i] = -1;
    }


    m_rowMap.resize(m_nblocks);
    m_solution.resize(m_nblocks);
    m_rhs.resize(m_nblocks);
    m_sparsity.resize2(m_nblocks,m_nblocks);
    m_matrix.resize2(m_nblocks,m_nblocks);
//    m_scratch.resize2(m_nblocks,m_nblocks);

#if USECPP11!=1
    for( int i=0 ; i<m_nblocks ; ++i )
    {
      m_rowMap[i] = NULL;
      m_solution[i] = NULL;
      m_rhs[i] = NULL;
      for( int j=0 ; j<m_nblocks ; ++j )
      {
        m_sparsity[i][j] = NULL;
        m_matrix[i][j] = NULL;
      }
    }
#endif

  }


  ~Epetra_System()
  {
#if USECPP11!=1
    for( int i=0 ; i<m_nblocks ; ++i )
    {
      delete m_rowMap[i];
      delete m_solution[i];
      delete m_rhs[i];
      for( int j=0 ; j<m_nblocks ; ++j )
      {
        delete m_sparsity[i][j];
        delete m_matrix[i][j];
      }
    }
#endif

  }

  void SetNumBlocks( const int numBlocks )
  {
    m_nblocks = numBlocks;
    m_rowMap.resize(m_nblocks);
    m_solution.resize(m_nblocks);
    m_rhs.resize(m_nblocks);
    m_sparsity.resize2(m_nblocks,m_nblocks);
    m_matrix.resize2(m_nblocks,m_nblocks);
//    m_scratch.resize2(m_nblocks,m_nblocks);
  }

  void SetBlockID( const int block, const EpetraBlock::ID id, const std::string& name )
  {
    if( block>=m_nblocks )
    {
      throw GPException("SolverBase.h:Epetra_System::SetBlockID(): block index is larger than number of blocks\n");
    }
    else if( block<0 )
    {
      throw GPException("SolverBase.h:Epetra_System::SetBlockID(): block index is less than 0\n");
    }

    if( m_blockID[id]==-1 )
    {
      m_blockID[id] = block;
    }
    else
    {
      throw GPException("error in Epetra_System::SetBlockID(). EpetraBlock::ID ("+toString(id)+") has already been used to register block "+ toString(
                          block)+ ". Ya betta check yoself befo ya wreck yoself.\n" );
    }

    if( m_solverNames[id].empty() )
    {
      m_solverNames[id] = name;
    }
    else
    {
      throw GPException("error in Epetra_System::SetBlockID(). EpetraBlock::ID ("+toString(id)+") has already been used to register solvername "+
                        m_solverNames[id]+ ". Ya betta check yoself befo ya wreck yoself.\n" );
    }

    std::map<std::string, int>::iterator iterSolverNameMap = m_solverNameMap.find(name);
    if( iterSolverNameMap==m_solverNameMap.end() )
    {
      m_solverNameMap[name] = block;
    }
    else
    {
      throw GPException("error in Epetra_System::SetBlockID(). Solver Name ("+name+") has already been used to register block "+ toString(
                          block)+ ". Ya betta check yoself befo ya wreck yoself.\n" );
    }
  }

  int GetBlockID( EpetraBlock::ID id ) const
  {
    return m_blockID[id];
  }

  std::string GetSolverName( EpetraBlock::ID id ) const
  {
    return m_solverNames[id];
  }


#if USECPP11==1
  std::shared_ptr<Epetra_Map>&
#else
  Epetra_Map*&
#endif
  GetRowMap( const EpetraBlock::ID dofID )
  {
    if( m_blockID[dofID]==-1 )
    {
      throw GPException("SolverBase.h:Epetra_System::GetRowMap():m_blockID isn't set \n");
    }
    return m_rowMap[ m_blockID[dofID] ];
  }

#if USECPP11==1
  std::shared_ptr<Epetra_FEVector>&
#else
  Epetra_FEVector*&
#endif
  GetSolutionVector( const EpetraBlock::ID dofID )
  {
    if( m_blockID[dofID]==-1 )
    {
      throw GPException("SolverBase.h:Epetra_System::GetSolutionVector():m_blockID isn't set \n");
    }
    return m_solution[ m_blockID[dofID] ];
  }

#if USECPP11==1
  std::shared_ptr<Epetra_FEVector>&
#else
  Epetra_FEVector*&
#endif
  GetResidualVector( const EpetraBlock::ID dofID )
  {
    if( m_blockID[dofID]==-1 )
    {
      throw GPException("SolverBase.h:Epetra_System::GetResidualVector():m_blockID isn't set \n");
    }
    return m_rhs[ m_blockID[dofID] ];
  }

#if USECPP11==1
  std::shared_ptr<Epetra_FECrsGraph>&
#else
  Epetra_FECrsGraph*&
#endif
  GetSparsity( const EpetraBlock::ID rowDofID,
               const EpetraBlock::ID colDofID )
  {
    if( m_blockID[rowDofID]==-1 && m_blockID[colDofID]==-1 )
    {
      throw GPException("SolverBase.h:Epetra_System::GetSparsity():m_blockID isn't set \n");
    }
    return m_sparsity[ m_blockID[rowDofID] ][ m_blockID[colDofID] ];
  }

#if USECPP11==1
  std::shared_ptr<Epetra_FECrsMatrix>&
#else
  Epetra_FECrsMatrix*&
#endif
  GetMatrix( const EpetraBlock::ID rowDofID,
             const EpetraBlock::ID colDofID )
  {
    if( m_blockID[rowDofID]==-1 && m_blockID[colDofID]==-1 )
    {
      throw GPException("SolverBase.h:Epetra_System::GetMatrix():m_blockID isn't set \n");
    }
    return m_matrix[ m_blockID[rowDofID] ][ m_blockID[colDofID] ];
  }

  bool HasMatrixBlock( const EpetraBlock::ID rowDofID,
                       const EpetraBlock::ID colDofID )
  {
    bool rval = false;
    if( m_blockID[rowDofID]!=-1 && m_blockID[colDofID]!=-1 )
    {
      rval = true;
    }
    return rval;
  }

/*
 #if USECPP11==1
   std::shared_ptr<Epetra_FECrsMatrix>
 #else
   Epetra_FECrsMatrix*&
 #endif
   GetScratch( const EpetraBlock::ID rowDofID, const EpetraBlock::ID colDofID )
   {
    if( m_blockID[rowDofID]==-1 && m_blockID[colDofID]==-1 )
    {
      throw GPException("SolverBase.h:Epetra_System::GetScratch():m_blockID
         isn't set \n");
    }
    return m_scratch[ m_blockID[rowDofID] ][ m_blockID[colDofID] ] ;
   }
 */
#if USECPP11==1
  array<std::shared_ptr<Epetra_Map> >         m_rowMap;
  array<std::shared_ptr<Epetra_FEVector> >    m_solution;
  array<std::shared_ptr<Epetra_FEVector> >    m_lastSolution;
  array<std::shared_ptr<Epetra_FEVector> >    m_rhs;
  Array2dT<std::shared_ptr<Epetra_FECrsGraph> >  m_sparsity;
  Array2dT<std::shared_ptr<Epetra_FECrsMatrix> > m_matrix;
//  Array2dT<std::shared_ptr<Epetra_FECrsMatrix> > m_scratch;
#else
  array<Epetra_Map*>         m_rowMap;
  array<Epetra_FEVector*>    m_solution;
  array<Epetra_FEVector*>    m_rhs;
  Array2dT<Epetra_FECrsGraph*>  m_sparsity;
  Array2dT<Epetra_FECrsMatrix*> m_matrix;
#endif



private:
  int m_nblocks;
  int m_blockID[EpetraBlock::numBlockDof];
  std::string m_solverNames[EpetraBlock::numBlockDof];
  std::map<std::string, int> m_solverNameMap;

  Epetra_System( Epetra_System& );
  Epetra_System& operator=(Epetra_System&);

};

class SolverBase
{
public:
  SolverBase(const std::string& name,
             ProblemManagerT* const pm );
  virtual ~SolverBase();

  virtual void TimeStepSetup( const realT& time,
                              const realT& dt,
                              PhysicalDomainT& domain,
                              SpatialPartition& partition,
                              const bool setupSystem )
  {
    (void)time;
    (void)dt;
    (void)domain;
    (void)partition;
    (void)setupSystem;
  }

  virtual double TimeStepExecute( const realT& time,
                                  const realT& dt,
                                  PhysicalDomainT& domain,
                                  SpatialPartition& partition )
  {
    (void)time;
    (void)dt;
    (void)domain;
    (void)partition;

    return dt;
  }

  virtual void TimeStepCleanup( PhysicalDomainT& domain, const realT& dt )
  {
    (void)domain;
  }

  virtual double TimeStep( const realT& time,
                           const realT& dt,
                           const int cycleNumber,
                           PhysicalDomainT& domain,
                           const array<string>& namesOfSolverRegions,
                           SpatialPartition& partition,
                           FractunatorBase* const fractunator ) = 0;

  virtual realT UpdateTimeStepMid(realT dt) { return dt; }

  virtual void PostProcess( PhysicalDomainT& domain,
                            SpatialPartition& partition,
                            const array<string>& namesOfSolverRegions);

  virtual void SetMaxStableTimeStep( const realT& time,
                                     PhysicalDomainT& domain,
                                     const array<string>& namesOfSolverRegions,
                                     SpatialPartition& partition);

  virtual void Initialize( PhysicalDomainT& domain, SpatialPartition& partition ) = 0;

  virtual void InitializeCommunications( PartitionBase& partition ) = 0;


  virtual void SetNumRowsAndTrilinosIndices( PhysicalDomainT& domain,
                                             SpatialPartition& partition,
                                             int& numLocalRows,
                                             int& numGlobalRows,
                                             array<integer>& localIndices,
                                             int offset )
  {
    (void)domain;
    (void)partition;
    (void)numLocalRows;
    (void)numGlobalRows;
    (void)localIndices;
    throw GPException("SolverBase::SetNumRowsAndTrilinosIndices() not overridden");
  }


  virtual void SetupSystem ( PhysicalDomainT& domain,
                             SpatialPartition& partition )
  {
    (void)domain;
    (void)partition;
    throw GPException("SolverBase::SetupSystem() not overridden");
  }

  virtual void SetSparsityPattern( PhysicalDomainT& domain )
  {
    (void)domain;
    throw GPException("SolverBase::SetSparsityPattern() not overridden");
  }

  /// returns name of specific solver instances
  std::string Name(){return m_name;};

  virtual void RegisterFields( PhysicalDomainT& domain ) = 0;

  virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn );

  virtual void WriteSilo( SiloFile& siloFile ) const;

  virtual void ReadSilo( const SiloFile& siloFile );

  StableTimeStep m_stabledt;
  realT m_courant, m_relaxationCoefficient;
  bool m_solverErrorFlag;


  std::map<PhysicalDomainT::ObjectDataStructureKeys, array<string> > m_syncedFields;

  std::string m_name;


#if USECPP11==1
  void SetRowMapPtr( std::shared_ptr<Epetra_Map>& rowMap )            { m_rowMap=rowMap; }
  void SetSparsityPtr( std::shared_ptr<Epetra_FECrsGraph>& sparsity ) { m_sparsity = sparsity; }
  void SetMatrixPtr( std::shared_ptr<Epetra_FECrsMatrix>& matrix )    { m_matrix = matrix; }
  void SetSolutionPtr( std::shared_ptr<Epetra_FEVector>& solution )   { m_solution = solution; }
  void SetRhsPtr( std::shared_ptr<Epetra_FEVector>& rhs )             { m_rhs = rhs; }
  void SetEpetraSystemPtr( std::shared_ptr<Epetra_System>& epetraSystem ) { m_epetraSystem = epetraSystem; }
#else
  void SetRowMapPtr( Epetra_Map* const rowMap )            { m_rowMap=rowMap; }
  void SetSparsityPtr( Epetra_FECrsGraph* const sparsity ) { m_sparsity = sparsity; }
  void SetMatrixPtr( Epetra_FECrsMatrix* const matrix )    { m_matrix = matrix; }
  void SetSolutionPtr( Epetra_FEVector* const solution )   { m_solution = solution; }
  void SetRhsPtr( Epetra_FEVector* const rhs )             { m_rhs = rhs; }
#endif


  void Solve ( PhysicalDomainT&  domain,
               SpatialPartition& partition,
               const realT time,
               const realT dt );


  void SolveBlock( PhysicalDomainT&  domain,
                   SpatialPartition& partition,
                   Epetra_System& epetraSystem,
                   const realT time,
                   const realT dt );

  // TODO these should be pure virtuals;
  virtual void SetInitialGuess( const PhysicalDomainT& domain,
                                realT* const local_solution )
  {
    (void)domain;
    (void)local_solution;
  }

  virtual void SetupMLPreconditioner( const PhysicalDomainT& domain, ML_Epetra::MultiLevelPreconditioner* MLPrec ) {}

  virtual void SetLinearSolverParameters( AztecOO& solver );

  realT ClearRow ( Epetra_FECrsMatrix* const matrix,
                   const unsigned int row,
                   const realT factor );


  virtual realT CheckSolution( const realT* const local_solution,
                               const PhysicalDomainT& domain,
                               const localIndex dofOffset )
  {
    (void)local_solution;
    (void)domain;
    (void)dofOffset;
    return 1.0;
  }

  virtual realT CheckSolutionBlock( const PhysicalDomainT& domain)
  {
    (void)domain;
    return 1.0;
  }


  virtual void PropagateSolution( const realT* const local_solution,
                                  const realT scalingFactor,
                                  PhysicalDomainT& domain,
                                  const localIndex dofOffset ) {}

  virtual void PropagateSolutionBlock( const realT scalingFactor,
                                       PhysicalDomainT& domain) {}

  virtual void SynchronizeFields( SpatialPartition& partition ) {}

  virtual void PostSyncConsistency( PhysicalDomainT& domain, SpatialPartition& partition ) {}

  const std::string& TrilinosIndexString() const
  {
    return m_trilinosIndexStr;
  }

  Epetra_System* m_system;

  struct
  {
    realT krylov_tol;       // Solver convergence criteria
    realT newton_tol;
    int   m_maxIters;       // Maximum number of solver iterations
    int   m_kspace;         // Number of krylov vectors before GMRES restart
    realT m_ilut_fill;      // Fill factor for ILUT preconditioner
    realT m_ilut_drop;      // Drop tolerance for ILUT preconditioner
    bool  m_useMLPrecond;   // Use ML preconditioner
    bool  m_useInnerSolver;  // Use row scaling
    int   m_scalingOption;  // Use row scaling
    bool  m_useBicgstab;    // Use bicgstab instead of gmres
    int   m_verbose;        // print extra info
    bool  m_useDirectSolver; // Use Direct solver
    bool  m_useNewtonSolve;  // Use Newton-Raphson iterations
    int   m_maxIterNewton;   // Maximum number of Newton-Raphson iterations
    int   m_numKrylovIter;
    realT m_KrylovResidualInit;
    realT m_KrylovResidualFinal;
  }
  m_numerics;

protected:

  std::string m_trilinosIndexStr;
  realT m_dtInit = 0.0;

  std::map<std::string,SolverBase*>& m_solvers;

public:
  array<string> m_commonFields;

protected:

#if GPAC_MPI
  const Epetra_MpiComm* epetra_comm;
#else
  const Epetra_SerialComm* epetra_comm;
#endif


#if USECPP11==1

  std::shared_ptr<Epetra_System> m_epetraSystem;

  std::shared_ptr<Epetra_Map>         m_rowMap;
  std::shared_ptr<Epetra_FECrsGraph>  m_sparsity;
  std::shared_ptr<Epetra_FECrsMatrix> m_matrix;
  std::shared_ptr<Epetra_FEVector>    m_solution;
  std::shared_ptr<Epetra_FEVector>    m_rhs;

#else

  Epetra_Map*         m_rowMap;
  Epetra_FECrsGraph*  m_sparsity;
  Epetra_FECrsMatrix* m_matrix;
  Epetra_FEVector*    m_solution;
  Epetra_FEVector*    m_rhs;

#endif

  double m_maxDofVal = 0.0;

  bool m_systemSetupFlag = true;
  bool m_systemAssemble = true;


private:

  std::string m_lsp;


  int m_timeIntegrationFlag;


  virtual void SetTimeIntegrationFlag( )  {}

  virtual void TimeStepDerived( const realT& time,
                                const realT& dt,
                                PhysicalDomainT& domain,
                                const array<string>& namesOfSolverRegions,
                                SpatialPartition& partition ) {};


  virtual void WriteSiloDerived( SiloFile& siloFile ) const {}
  virtual void ReadSiloDerived( const SiloFile& siloFile ) {}

};


#endif /* SOLVERBASE_H_ */
