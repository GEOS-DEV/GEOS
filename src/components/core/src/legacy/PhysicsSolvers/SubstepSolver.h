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
 * @file SubstepSolver.h
 * @author walsh24
 * @date January 7, 2012
 */

#ifndef SUBSTEPSOLVER_H_
#define SUBSTEPSOLVER_H_

#include "SolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ProblemManagerT.h"

#include "Common/Common.h"

class SubstepSolver : public SolverBase
{
public:
  SubstepSolver(const std::string& name,
                ProblemManagerT* const pm );
  ~SubstepSolver();

  void ReadXML( TICPP::HierarchicalDataNode* const hdn);

  double TimeStep( const realT& time,
                   const realT& dt,
                   const int cycleNumber,
                   PhysicalDomainT& domain,
                   const array<string>& namesOfSolverRegions,
                   SpatialPartition& partition,
                   FractunatorBase* const fractunator );

  void RegisterFields( PhysicalDomainT& domain  ){};
  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition  ){};

  void InitializeCommunications( PartitionBase& partition  ) {};


  static const char* SolverName(){return "SubstepSolver";};
  static unsigned NumberOfInstances(void);


private:

  ProblemManagerT* m_problemManagerPtr;
  std::map<std::string,SolverBase*>* m_solverMapPtr;
  array<string> m_SolverNames;
  std::vector<array<string> >m_namesOfSolverRegions;

  realT m_dt;  // maximum substep timestep
  unsigned m_depth;
  unsigned m_subcycleNumber; // number of subcycles from start of simulation
  static unsigned m_instances;
};


#endif
