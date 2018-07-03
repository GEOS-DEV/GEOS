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
 * @file ParallelPlateFlowSolverFEM.h
 * @author walsh24
 * @date June 1, 2011
 */

#ifndef PARALLELPLATEFLOWSOLVERFEM_H_
#define PARALLELPLATEFLOWSOLVERFEM_H_

#include "PhysicsSolvers/ParallelPlateFlowSolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ProblemManagerT.h"

#include "Common/Common.h"
#include "Utilities/TrilinosUtilities.h"
#include "Utilities/RCVSparse.h"

#include "PhysicsSolvers/PhysicsSolverStrings.h"

#if GPAC_MPI
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_LinearProblem.h"

#include "EpetraExt_RowMatrixOut.h"
#include "Teuchos_RCP.hpp"
#include "AztecOO.h"
#include "ml_MultiLevelPreconditioner.h"


class ParallelPlateFlowSolverFEM : public ParallelPlateFlowSolverBase
{
public:
  ParallelPlateFlowSolverFEM(  const std::string& name,
                               ProblemManagerT* const pm );
  virtual ~ParallelPlateFlowSolverFEM();

  void ReadXML( TICPP::HierarchicalDataNode* const hdn );
  void RegisterFields( PhysicalDomainT& domain );
  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);

  void InitializeCommunications( PartitionBase& partition );

  void TimeStep( const realT& time, const realT& dt,PhysicalDomainT& domain,
                 const array<string>& namesOfSolverRegions, SpatialPartition& partition,
                 FractunatorBase* const fractunator );

  void DefineFlowSets( PhysicalDomainT& domain );

  void CalculateMassRate(PhysicalDomainT& domain, SpatialPartition& partition,realT time,realT dt );
  void UpdateEOS( PhysicalDomainT& domain,const realT dt, const bool updateMass );
  void UpdateFlux( const realT time, PhysicalDomainT& domain);


  /// name of solver class
  /// NB: Parent class's Name() function returns name of specific solver objects
  static const char* SolverName(){return "ParallelPlateFlowSolverFEM";};



  // Flags
  bool m_doApertureUpdate;

private:


  void SetupSystem (PhysicalDomainT& domain, SpatialPartition& partition, const realT& time);
  void Assemble    (PhysicalDomainT& domain, SpatialPartition& partition, const realT& time, const realT& dt);
  void Solve       (PhysicalDomainT& domain, SpatialPartition& partition);

  lSet m_faceSet;


  realT m_phi; // Mixed Euler parameter (0 = forward difference, 0.5 = central
               // difference, 1.0 = backward difference)


  // MPI
  const int this_mpi_process;
  const int n_mpi_processes;

  #if GPAC_MPI
  const Epetra_MpiComm & m_epetra_comm;
  #else
  const Epetra_SerialComm & m_epetra_comm;
  #endif

  Teuchos::RCP<Epetra_Map>         row_map;
  Teuchos::RCP<Epetra_FECrsGraph>  sparsity;
  Teuchos::RCP<Epetra_FECrsMatrix> matrix;
  Teuchos::RCP<Epetra_FEVector>    solution;
  Teuchos::RCP<Epetra_FEVector>    rhs;


  std::string m_TrilinosIndexStr;



};


#endif /* ParallelPlateFlowSolverFEM_IM_H_ */
