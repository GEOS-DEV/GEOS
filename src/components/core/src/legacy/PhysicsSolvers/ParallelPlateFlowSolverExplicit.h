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
 * @file ParallelPlateFlowSolverExplicit.h
 * @author settgast1
 * @date Feb 10, 2011
 */

#ifndef PARALLELPLATEFLOWSOLVEREXPLICIT_H_
#define PARALLELPLATEFLOWSOLVEREXPLICIT_H_

#include "PhysicsSolvers/ParallelPlateFlowSolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "IO/ticpp/TinyXMLParser.h"

class ParallelPlateFlowSolverExplicit : public ParallelPlateFlowSolverBase
{
public:
  ParallelPlateFlowSolverExplicit( const std::string& name,
                                   ProblemManagerT* const pm );
  virtual ~ParallelPlateFlowSolverExplicit();

  double TimeStep( const realT& time,
                   const realT& dt,
                   const int cycleNumber,
                   PhysicalDomainT& domain,
                   const array<string>& namesOfSolverRegions,
                   SpatialPartition& partition,
                   FractunatorBase* const fractunator );

  void PostProcess (PhysicalDomainT& domain,
                    SpatialPartition& partition,
                    const array<string>& namesOfSolverRegions);
  virtual void RegisterFields( PhysicalDomainT& domain );

  void GenerateParallelPlateGeometricQuantities( PhysicalDomainT& domain,
                                                 realT time,realT dt );

  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);

  void InitializeCommunications( PartitionBase& partition );

  virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn );


  /// name of solver class
  /// NB: Name() function of parent class returns name of specific solver
  // objects
  static const char* SolverName(){return "ParallelPlateFlowSolverExplicit";};

  void CalculateAndApplyMassFlux( const realT dt, PhysicalDomainT& domain );

  void CalculateCarterLeakOff( const realT time, const realT dt, PhysicalDomainT& domain );

  void CalculateMatrixFlowLeakOff( const realT time, const realT dt, PhysicalDomainT& domain );

  void ApplyFluxBoundaryCondition( const realT time, const realT dt, const int cycleNumber, const int rank, PhysicalDomainT& domain );

  void FlowControlledBoundaryCondition( PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                                        BoundaryConditionBase* bc,
                                        const lSet& set,
                                        realT time );

  void UpdateEOS( const realT time, const realT dt, PhysicalDomainT& domain );

  void CalculateNodalPressure ( PhysicalDomainT& domain, SpatialPartition& partition);

private:
  realT m_bBarton;
  realT m_aBarton;
  realT m_wZeroStress;
  realT m_dT;
  realT m_farFieldPorePressure;
  int m_pressureDependentLeakoff;
  realT m_apertureMovingAverageCoeff;


};

#endif /* PARALLELPLATEFLOWSOLVEREXPLICIT_H_ */
