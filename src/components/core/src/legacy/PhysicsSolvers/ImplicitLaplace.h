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

#ifndef IMPLICIT_LAPLACE_H
#define IMPLICIT_LAPLACE_H

/**
 * @file ImplicitLaplace.h
 * @author white230
 */

#include "Common/Common.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ProblemManagerT.h"
#include "ObjectManagers/FunctionManager.h"
#include "ObjectManagers/FaceManagerT.h"
#include "ElementLibrary/FiniteElement.h"
#include "ElementLibrary/LagrangeBasis.h"
#include "ElementLibrary/GaussQuadrature.h"
#include "PhysicsSolvers/SolverBase.h"
#include "PhysicsSolvers/SolverFactory.h"
#include "PhysicsSolverStrings.h"

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

/**
 * Physics solver for Laplace's equation
 */

template <int dim>
class ImplicitLaplaceSolver : public SolverBase
{
public:
  ImplicitLaplaceSolver( const std::string& name,
                         ProblemManagerT* const pm );

  ~ImplicitLaplaceSolver();

  static const char* SolverName();

  virtual void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);

  virtual void RegisterFields(PhysicalDomainT &domain);

  virtual void InitializeCommunications(PartitionBase& partition);

  virtual double TimeStep(const realT& time,
                          const realT& dt,
                          const int cycleNumber,
                          PhysicalDomainT * domain,
                          const array<string>& namesOfSolverRegions,
                          SpatialPartition& partition,
                          FractunatorBase* const fractunator);

private:

  void ReadXML(TICPP::HierarchicalDataNode* const hdn);

  void SetupSystem (PhysicalDomainT& domain, SpatialPartition& partition);
  void Assemble    (PhysicalDomainT& domain, SpatialPartition& partition);
  void Solve       (PhysicalDomainT& domain, SpatialPartition& partition);

    #if GPAC_MPI
  const Epetra_MpiComm & epetra_comm;
    #else
  const Epetra_SerialComm & epetra_comm;
    #endif

  const int this_mpi_process;
  const int n_mpi_processes;
  const bool     verbose;

  struct
  {
    double diffusion;
    double source;
  }
  equation_data;

  struct
  {
    double krylov_tol;
  }
  numerics;

  typedef std::map<ElementManagerT::RegKeyType, ElementRegionT > RegionMap;

  Epetra_Map*         row_map;
  Epetra_FECrsGraph*  sparsity;
  Epetra_FECrsMatrix* matrix;
  Epetra_FEVector*    solution;
  Epetra_FEVector*    rhs;

  std::map<PhysicalDomainT::ObjectDataStructureKeys, array<string> > syncedFields;
};


#endif
