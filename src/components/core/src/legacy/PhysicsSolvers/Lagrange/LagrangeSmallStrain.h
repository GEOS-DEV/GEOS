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

#ifndef LAGRANGESMALLSTRAIN_H_
#define LAGRANGESMALLSTRAIN_H_

#include "LagrangeSolverBase.h"



class LagrangeSmallStrain : public LagrangeSolverBase
{
public:
  LagrangeSmallStrain(  const std::string& name,
                        ProblemManagerT* const pm );
  virtual ~LagrangeSmallStrain();


  virtual void Initialize(PhysicalDomainT& domain, SpatialPartition& partition  );

  void InitializeCommunications( PartitionBase& partition );

  virtual void
  RegisterFields(PhysicalDomainT& domain);

  virtual void
  ReadXML( TICPP::HierarchicalDataNode* const hdn );

  static const char*
  SolverName()
  {
    return "LagrangeSmallStrain";
  };

  // boundary conditions
  virtual void TractionBC(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                          BoundaryConditionBase* bc, const lSet& set, realT time);

  virtual void PressureBC(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                          BoundaryConditionBase* bc, const lSet& set, realT time);

  virtual void DisplacementBC(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                              BoundaryConditionBase* bc, const lSet& set, realT time);

  virtual void FixNodesBC( NodeManager const& nodeManager, const lSet& set );


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


  virtual void SetInitialGuess( const PhysicalDomainT& domain,
                                realT* const local_solution  );

  virtual void PropagateSolution( const realT* const local_solution,
                                  const realT scalingFactor,
                                  PhysicalDomainT& domain,
                                  const localIndex dofOffset  );
  bool m_staticKMatrix;


private:


  virtual void ProcessElementRegion( NodeManager& nodeManager,
                                     ElementRegionT& elemRegion,
                                     const realT dt )
  {
    throw GPException( "LagrangeSmallStrain::ProcessElementRegions() not implemented");
  }



  realT m_nonContactModulus;
  bool m_doApertureUpdate;
  bool m_recordIncrementalDisplacement;
  bool m_hasReferenceStress;

};



#endif /* LAGRANGESMALLSTRAIN_H_ */
