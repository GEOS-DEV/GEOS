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
 * @file ImplicitMechanicsSolver.h
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */

#ifndef ImplicitMechanicsSolver_H_
#define ImplicitMechanicsSolver_H_

#include "SolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "IO/ticpp/TinyXMLParser.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"
#ifdef SRC_EXTERNAL
#include "Contact/CommonPlaneContact.h"
#endif

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



template <int dim>
class ImplicitMechanicsSolver : public SolverBase
{
public:
  ImplicitMechanicsSolver( const std::string& name,
                           ProblemManagerT* const pm );
  ~ImplicitMechanicsSolver();

  double
  TimeStep(const realT& time, const realT& dt,
           const int cycleNumber,
           PhysicalDomainT& domain,
           const array<string>& namesOfSolverRegions, SpatialPartition& partition,
           FractunatorBase* const fractunator);


  void Initialize(PhysicalDomainT& domain, SpatialPartition& partition  );

  void InitializeCommunications( PartitionBase& partition );

  virtual void
  RegisterFields(PhysicalDomainT& domain);

  virtual void
  ReadXML( TICPP::HierarchicalDataNode* const hdn );

  static const char*
  SolverName()
  {
    return "ImplicitMechanicsSolver3D";
  };

  // boundary conditions
  virtual void TractionBC(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                          BoundaryConditionBase* bc, const lSet& set, realT time);

  virtual void PressureBC(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                          BoundaryConditionBase* bc, const lSet& set, realT time);

  virtual void DisplacementBC(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                              BoundaryConditionBase* bc, const lSet& set, realT time);



  virtual void NonpenetratingBC_NeighborUpdate(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                                               BoundaryConditionBase* bc, realT time);
  virtual void NonpenetratingBC_DetectContact(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                                              BoundaryConditionBase* bc, realT time);

  virtual void NonpenetratingBC_Sparsity(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                                         BoundaryConditionBase* bc, realT time);

  virtual void NonpenetratingBC_Apply(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                                      BoundaryConditionBase* bc, realT time);

  virtual void NonpenetratingBC_UpdateAperture(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                                               BoundaryConditionBase* bc, realT time);


private:

  void SetupSystem (PhysicalDomainT& domain, SpatialPartition& partition);
  void Assemble    (PhysicalDomainT& domain, SpatialPartition& partition, realT time);
  void Solve       (PhysicalDomainT& domain, SpatialPartition& partition, realT time);

  #if GPAC_MPI
  const Epetra_MpiComm & epetra_comm;
  #else
  const Epetra_SerialComm & epetra_comm;
  #endif

  const unsigned this_mpi_process;
  const unsigned n_mpi_processes;
  const bool     verbose;

  bool m_useMLPrecond;

  struct
  {
    double lambda; // Lame's constant
    double G;      // Shear modulus
  }
  equation_data;

  struct
  {
    double krylov_tol;
  }
  numerics;

  typedef std::map<ElementManagerT::RegKeyType, ElementRegionT > RegionMap;

  Teuchos::RCP<Epetra_Map>         row_map;
  Teuchos::RCP<Epetra_FECrsGraph>  sparsity;
  Teuchos::RCP<Epetra_FECrsMatrix> matrix;
  Teuchos::RCP<Epetra_FEVector>    solution;
  Teuchos::RCP<Epetra_FEVector>    rhs;
  array<integer> dummyDof;

  std::string m_trilinosIndexStr;
  static int m_instances;


  realT m_cfl;
  realT m_nonContactModulus;
  bool m_doApertureUpdate;

  bool m_recordIncrementalDisplacement;

  std::map<PhysicalDomainT::ObjectDataStructureKeys, array<string> > syncedFields;
};

template <int dim>
int ImplicitMechanicsSolver<dim>::m_instances = 0;


#endif /* ImplicitMechanicsSolver_H_ */
