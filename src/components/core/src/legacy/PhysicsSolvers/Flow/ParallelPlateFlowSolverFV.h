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
 * @file ParallelPlateFlowSolverFV.h
 * @author walsh24
 * @date June 1, 2011
 */

#ifndef PARALLELPLATEFLOWSOLVER_IM_H_
#define PARALLELPLATEFLOWSOLVER_IM_H_

#include "../ParallelPlateFlowSolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ProblemManagerT.h"

#include "Common/Common.h"
#include "Common/typedefs.h"
#include "Utilities/TrilinosUtilities.h"
#include "Utilities/RCVSparse.h"

#include "../PhysicsSolverStrings.h"


#ifdef SRC_EXTERNAL
#include "BoundaryConditions/PerforatedCasedWellboreBoundaryCondition.h"
#include "BoundaryConditions/CavityPressureBoundaryCondition.h"
#endif

class ParallelPlateFlowSolverFV : public ParallelPlateFlowSolverBase
{
public:
  ParallelPlateFlowSolverFV(  const std::string& name,
                              ProblemManagerT* const pm );
  virtual ~ParallelPlateFlowSolverFV();

  void ReadXML( TICPP::HierarchicalDataNode* const hdn );
  void RegisterFields( PhysicalDomainT& domain );
  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);

  void InitializeCommunications( PartitionBase& partition );

  void RegisterTemporaryFields( PhysicalDomainT& domain );
  void DeregisterTemporaryFields( PhysicalDomainT& domain );
  void FillTemporaryFields( PhysicalDomainT& domain );
  void OverwriteFieldsWithTemporaryFields( PhysicalDomainT& domain );

  virtual void TimeStepSetup( const realT& time,
                              const realT& dt,
                              PhysicalDomainT& domain,
                              SpatialPartition& partition,
                              const bool setupSystem  );

  virtual double TimeStepExecute( const realT& time,
                                  const realT& dt,
                                  PhysicalDomainT& domain,
                                  SpatialPartition& partition );

  virtual void TimeStepCleanup( PhysicalDomainT& domain, const realT& dt );

  double TimeStep( const realT& time, const realT& dt,
                   const int cycleNumber,
                   PhysicalDomainT& domain,
                   const array<string>& namesOfSolverRegions, SpatialPartition& partition,
                   FractunatorBase* const fractunator );

  void DefineFlowSets( PhysicalDomainT& domain );

  void CalculateMassRate(PhysicalDomainT& domain, SpatialPartition& partition,realT time,realT dt );
  void UpdateEOS( const realT time, const realT dt, PhysicalDomainT& domain );
  void UpdateFlux( const realT time,
                   const realT dt,
                   PhysicalDomainT& domain, SpatialPartition& partition);

  // Pressure controlled boundary condition
  void PressureBoundaryCondition(PhysicalDomainT& domain,
                                 ObjectDataStructureBaseT& object,
                                 BoundaryConditionBase* const bc,
                                 const lSet& set,
                                 const realT time,
                                 const realT dt,
                                 const int dofOffset );

  void PressureBoundaryCondition_VelocityUpdate(PhysicalDomainT& domain,
                                                ObjectDataStructureBaseT& object,
                                                BoundaryConditionBase* const bc,
                                                const lSet& set,
                                                const realT time,
                                                const realT dt);

  void FaceFluxBoundaryCondition( PhysicalDomainT& domain,ObjectDataStructureBaseT& object,
                                  BoundaryConditionBase* const bc, const lSet& aset, const realT time, const realT dt );

  virtual void MassBC(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                      BoundaryConditionBase* bc, const lSet& set, realT time);

  virtual void MassRateBC( PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                           BoundaryConditionBase* bc, const lSet& set, realT time, realT dt);

  virtual void VolumeRateBC( PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                             BoundaryConditionBase* bc, const lSet& set, realT time, realT dt);

  void SetConcentrationField(std::string& concentrationFieldName,std::string reactionRateFieldName,std::string rrDerivFieldName);

  /// name of solver class
  /// NB: Parent class's Name() function returns name of specific solver objects
  static const char* SolverName(){return "ParallelPlateFlowSolverFV";};

  void GenerateParallelPlateGeometricQuantities( PhysicalDomainT& domain,realT time,realT dt );

  // Boundary conditions
//  realT m_dt; // current dt

  // Flags
  bool m_doApertureUpdate;

  void SetNumRowsAndTrilinosIndices( PhysicalDomainT& domain,
                                     SpatialPartition& partition,
                                     int& numLocalRows,
                                     int& numGlobalRows,
                                     array<integer>& localIndices,
                                     int offset );


  void Assemble    (PhysicalDomainT& domain, SpatialPartition& partition, const realT& time, const realT& dt)
  {
    throw GPException("ParallelPlateFlowSolverBase::Assemble() not overridden");
  }


  realT Assemble    ( PhysicalDomainT& domain,
                      Epetra_System& epetraSystem,
                      const realT& time,
                      const realT& dt );

  realT AssembleCD( PhysicalDomainT& domain,
                    Epetra_System& epetraSystem,
                    const realT& time,
                    const realT& dt );

  using SolverBase::SetInitialGuess;
  virtual void SetInitialGuess( const PhysicalDomainT& domain,
                                realT* const local_solution );

  virtual realT CheckSolution( const realT* const local_solution,
                               const PhysicalDomainT& domain,
                               const localIndex dofOffset );

  virtual void PropagateSolution( const realT* const local_solution,
                                  const realT scalingFactor,
                                  PhysicalDomainT& domain,
                                  const localIndex dofOffset  );

  virtual void PostSyncConsistency( PhysicalDomainT& domain, SpatialPartition& partition );


/*
   void Solve ( PhysicalDomainT&  domain,
               SpatialPartition& partition,
               const realT time,
               const realT dt );
 */


  virtual void CalculateAndApplyMassFlux( const realT dt, PhysicalDomainT& domain );


  virtual void CalculateCarterLeakOff( const realT time, const realT dt, PhysicalDomainT& domain );

  virtual void ApplyFluxBoundaryCondition( const realT time, const realT dt, const int cycleNumber, const int rank, PhysicalDomainT& domain );

  virtual void FlowControlledBoundaryCondition( PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                                                BoundaryConditionBase* bc,
                                                const lSet& set,
                                                realT time );


  virtual void CalculateApertureDerivatives( const FaceManagerT& faceManager,
                                             const NodeManager& nodeManager );
private:


  void SetupSystem (PhysicalDomainT& domain, SpatialPartition& partition );

  void InitializeDensity( PhysicalDomainT& domain);
  void UpdateAperture(PhysicalDomainT&  domain);

  realT TwoFacePermeability_FV(const array<R1Tensor>& edgeCenters,
                               const array<real64>& edgeLengths,
                               const array<R1Tensor>& faceCenters,
                               const array<real64>& apertures,
                               const localIndex eg,
                               const localIndex kf,
                               const localIndex kfb,
                               const array<array<real64> >* const dwdu,
                               array<real64>* const dkdu_r,
                               array<real64>* const dkdu_s);

  lSet m_faceSet;
  localIndex m_numFaces;
  std::map<localIndex,localIndex> m_faceDofMap;

  realT m_negativePressureSlope;


  realT m_phi; // Mixed Euler parameter (0 = forward difference, 0.5 = central
               // difference, 1.0 = backward difference)



  // MPI
  const int this_mpi_process;
  const int n_mpi_processes;

  std::map<PhysicalDomainT::ObjectDataStructureKeys, array<string> > syncedFields;

  static int m_instances;
  realT m_gamma;

  array<locallyIndexedReal> indexedFluxes;


};


#endif /* PARALLELPLATEFLOWSOLVER_IM_H_ */
