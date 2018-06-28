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
 * @file LagrangeExplicitDynamicsSolver.h
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */

#ifndef LAGRANGEEXPLICITDYNAMICSSOLVER_H_
#define LAGRANGEEXPLICITDYNAMICSSOLVER_H_

#include "SolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "IO/ticpp/TinyXMLParser.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"
#ifdef SRC_EXTERNAL
#include "Contact/CommonPlaneContact.h"
#endif

class LagrangeExplicitDynamicsSolver : public SolverBase
{
public:
  LagrangeExplicitDynamicsSolver( const std::string& name,
                                  ProblemManagerT* const pm );
  ~LagrangeExplicitDynamicsSolver();

  virtual double
  TimeStep(const realT& time, const realT& dt,
           const int cycleNumber,
           PhysicalDomainT& domain,
           const array<string>& namesOfSolverRegions, SpatialPartition& partition,
           FractunatorBase* const fractunator);

  void PostProcess (PhysicalDomainT& domain,
                    SpatialPartition& partition,
                    const array<string>& namesOfSolverRegions);

  virtual void SetMaxStableTimeStep( const realT& time,
                                     PhysicalDomainT& domain,
                                     const array<string>& namesOfSolverRegions,
                                     SpatialPartition& partition);

  void Initialize(PhysicalDomainT& domain, SpatialPartition& partition  );

  virtual void InitializeCommunications( PartitionBase& partition );

  virtual void
  RegisterFields(PhysicalDomainT& domain);

  //  void LagrangianNodeUpdatePart1of2( ObjectDataStructureBaseT& nodeManager,
  // const realT& time, const realT& dt );
  //  void LagrangianNodeUpdatePart2of2( ObjectDataStructureBaseT& nodeManager,
  // const realT& time, const realT& dt );


  virtual void
  ReadXML( TICPP::HierarchicalDataNode* const hdn);

  static const char*
  SolverName()
  {
    return "LagrangeExplicitDynamicsSolver";
  }
  ;


  void ApplyGapDamping( NodeManager& nodeManager,
                        const FaceManagerT& faceManager,
                        const realT dt );


  realT m_dampingM;
  realT m_gapdamping;
  array<lSet> m_KinematicConstraintNodes;
//  array<std::pair<localIndex,localIndex> > m_KinematicConstraintFaces;

  bool m_tiedNodesFlag;
  realT m_tiedNodeNormalRuptureStress;
  realT m_tiedNodeShearRuptureStress;
  realT m_tiedNodeTolerance;

  void ApplyForcesFromContact(PhysicalDomainT& domain,
                              StableTimeStep& timeStep,
                              const realT dt );


protected:


private:
  virtual void WriteSiloDerived( SiloFile& siloFile ) const;
  virtual void ReadSiloDerived( const SiloFile& siloFile );



};

namespace LagrangeExplicitDynamicsFunctions
{

void LinearPointUpdatePart1( ObjectDataStructureBaseT& objectManager,
                             const realT& time,
                             const realT& dt,
                             const bool clearForces = true);

void LinearPointUpdatePart2( ObjectDataStructureBaseT& objectManager,
                             const realT& time,
                             const realT& dt,
                             const realT& damping = 0.0 );

void RotationalPointUpdatePart1(DiscreteElementManagerBaseT& discreteElementManager,
                                const realT& time, const realT& dt);

void RotationalPointUpdatePart1b(DiscreteElementManagerT& discreteElementManager);


void RotationalPointUpdatePart2(ObjectDataStructureBaseT& objectManager,
                                const realT& time, const realT& dt);
}

#endif /* LAGRANGEEXPLICITDYNAMICSSOLVER_H_ */
