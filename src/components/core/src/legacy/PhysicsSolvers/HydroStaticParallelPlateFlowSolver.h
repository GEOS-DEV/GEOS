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
 * @file HydroStaticParallelPlateFlowSolver.h
 * @author fu4
 * @date Mar 26, 2014
 * @brief Assuming pressure gradient can be ignored.  Mainly for toughness
 * dominated simulations. The total volume and total mass are used to work out
 * the pressure.
 */

#ifndef HYDROSTATICPARALLELPLATEFLOWSOLVER_H_
#define HYDROSTATICPARALLELPLATEFLOWSOLVER_H_

#include "PhysicsSolvers/ParallelPlateFlowSolverExplicit.h"
#include "ObjectManagers/PhysicalDomainT.h"

class HydroStaticParallelPlateFlowSolver : public ParallelPlateFlowSolverExplicit
{
public:
  HydroStaticParallelPlateFlowSolver( const std::string& name,
                                      ProblemManagerT* const pm );
  virtual ~HydroStaticParallelPlateFlowSolver();

  void PostProcess (PhysicalDomainT& domain,
                    SpatialPartition& partition,
                    const array<string>& namesOfSolverRegions);

  virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn );

  /// name of solver class
  /// NB: Name() function of parent class returns name of specific solver
  // objects
  static const char* SolverName(){return "HydroStaticParallelPlateFlowSolver";};

  void CalculateAndApplyMassFlux( const realT dt, PhysicalDomainT& domain );


  void UpdateEOS( const realT time, const realT dt, PhysicalDomainT& domain );
  void ApplyBoreholePressure( const realT time, const realT dt, realT pressure, PhysicalDomainT& domain );

private:
  realT m_cavityVolume;
  realT m_cavityMass;
  array<string> m_boreholeSetNames;

};

#endif /* HYDROSTATICPARALLELPLATEFLOWSOLVER_H_ */
