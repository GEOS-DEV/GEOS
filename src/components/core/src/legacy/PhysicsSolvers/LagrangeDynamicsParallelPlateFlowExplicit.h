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
 * @file LagrangeDynamicsParallelPlateFlowExplicit.h
 * @author settgast1
 * @date Feb 28, 2011
 */

#ifndef LAGRANGEDYNAMICSPARALLELPLATEFLOWEXPLICIT_H_
#define LAGRANGEDYNAMICSPARALLELPLATEFLOWEXPLICIT_H_

#include "SolverBase.h"
#include "LagrangeExplicitDynamicsSolver.h"
#include "ParallelPlateFlowSolverExplicit.h"

class LagrangeDynamicsParallelPlateFlowExplicit : public SolverBase
{
public:
  LagrangeDynamicsParallelPlateFlowExplicit( const std::string& name,
                                             ProblemManagerT* const pm );
  virtual ~LagrangeDynamicsParallelPlateFlowExplicit();

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

  virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn );

  virtual void RegisterFields( PhysicalDomainT& domain );

  void Initialize(PhysicalDomainT& domain, SpatialPartition& partition );

  void InitializeCommunications( PartitionBase& partition );

  void ApplyForcesFromContact(PhysicalDomainT& domain,
                              StableTimeStep& timeStep,
                              const realT dt );

  /// name of solver class
  /// NB: Name() function of parent class returns name of specific solver
  // objects
  static const char* SolverName(){return "LagrangeDynamicsParallelPlateFlowExplicit";};


private:
  LagrangeExplicitDynamicsSolver m_ldSolve;
  ParallelPlateFlowSolverExplicit m_ppSolve;

  realT m_kJn, m_kJs;
  realT m_COFJ;
  realT m_fLockedInSIF;
  realT m_faceStrengthRandomFactor;

  void CalculateContactStress(PhysicalDomainT& domain,
                              const realT& time,
                              const realT& dt,
                              localIndex& kf,
                              const localIndex faceIndex[],
                              realT& stressPen,
                              R1Tensor& stressShear);

};



#endif /* LAGRANGEDYNAMICSPARALLELPLATEFLOWEXPLICIT_H_ */
