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
 * @file SteadyStateADRSolver.cpp
 * @author walsh24
 * @date January 7, 2012
 */


#include "SolverFactory.h"
#include "SubstepSolver.h"
#include "Common/Common.h"
#include "Common/intrinsic_typedefs.h"
#include "Utilities/StringUtilities.h"

#include <cstdlib>


unsigned SubstepSolver::m_instances = 0;

/// Steady-State Advection diffusion reaction solver
SubstepSolver::SubstepSolver( const std::string& name,
                              ProblemManagerT* const pm ):
  SolverBase(std::string("SubstepSolver_" + toString(m_instances)), pm ),
  m_depth(0),
  m_subcycleNumber(0)
{
  m_solverMapPtr = &(pm->m_solvers);
  m_problemManagerPtr = pm;
  ++m_instances;
}

SubstepSolver::~SubstepSolver()
{
  // TODO Auto-generated destructor stub
}

/*
   <SolverApplications>
    <Application name="1" begintime="0.0" endtime="2 us" dt="1us">
       <Apply solver="ledSolver" toregions="Fracture_Floor Fracture_Roof" />
       <Substep dt = "0.1 us">
         <Apply solver="pvSolver" toregions="Fracture_Floor" />
       <Substep/>
    </Application>
   </SolverApplications>
 */
void SubstepSolver::ReadXML(TICPP::HierarchicalDataNode* const hdn)
{
  SolverBase::ReadXML( hdn );

  m_dt = hdn->GetAttributeOrDefault<realT>("dt",std::numeric_limits<double>::max());
  m_depth = hdn->GetAttributeOrDefault<realT>("depth",0);
  std::map<std::string,SolverBase*>& solverMap = *m_solverMapPtr;

  for (HierarchicalDataNode* applicationNode = hdn->Next(true) ; applicationNode ; applicationNode = hdn->Next())
  {

    const std::string& apType = applicationNode->Heading();

    std::string solverName;
    array<string> regionNames;

    if( streq(apType,"Substep") )
    {
      // Recursive substep solver
      solverName = "SubstepSolver_" + toString(m_instances);
      applicationNode->AddAttributePair("name", solverName );
      applicationNode->AddAttributePair("depth", toString(m_depth+1) );
      solverMap[solverName] = SolverFactory::NewSolver(SubstepSolver::SolverName(), applicationNode, m_problemManagerPtr);

    }
    else
    {
      solverName = applicationNode->GetAttributeString("solver");
      regionNames = applicationNode->GetStringVector("toregions");
    }
    m_SolverNames.push_back(solverName);
    m_namesOfSolverRegions.push_back(regionNames);

  }
}

unsigned SubstepSolver::NumberOfInstances(void){return m_instances;}

/**
 *
 *
 **/
double SubstepSolver::TimeStep( const realT& time,
                                const realT& dt,
                                const int cycleNumber,
                                PhysicalDomainT& domain,
                                const array<string>& namesOfSolverRegions,
                                SpatialPartition& partition,
                                FractunatorBase* const fractunator )
{

  realT dt_return = dt;

  std::map<std::string,SolverBase*>& solverMap = *m_solverMapPtr;
  unsigned numSolvers = m_SolverNames.size();

  m_stabledt.m_maxdt = 0.9*std::numeric_limits<double>::max(); // Substep solver
                                                               // does not
                                                               // change macro
                                                               // timestep (only
                                                               // internal)

  realT t = time;
  realT endTime = time+dt;

  realT almostEndTime = endTime*(1.0 - std::numeric_limits<realT>::epsilon());
  while( t < almostEndTime)
  {
    // Determine stable substep timestep.
    realT subDt = endTime-t;  // ensures substep and macrostep are in sync


    if(m_dt < subDt)
      subDt = m_dt;
    for(unsigned i = 0 ; i < numSolvers ; ++i)
    {
      SolverBase* solverPtr = solverMap[m_SolverNames[i]];
      if( solverPtr->m_stabledt.m_maxdt < subDt )
        subDt = solverPtr->m_stabledt.m_maxdt;
    }

    realT dummy = subDt;
    MPI_Allreduce (&dummy,&subDt,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);

    if (partition.m_rank == 0)
    {
      for(unsigned i = 0 ; i < m_depth ; ++i)
        std::cout << "    ";
      std::cout<<"    Substep: Time=" << t << ", dt=" << subDt << std::endl;
    }

    // Perform substep
    for(unsigned i = 0 ; i < numSolvers ; ++i)
    {
      SolverBase* solverPtr = solverMap[m_SolverNames[i]];
      solverPtr->TimeStep( t, subDt, m_subcycleNumber, domain, m_namesOfSolverRegions[i], partition, fractunator );
    }
    t += subDt;
    m_subcycleNumber += 1;
  }
  return dt_return;
}

REGISTER_SOLVER( SubstepSolver )
