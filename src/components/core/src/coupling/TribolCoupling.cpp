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
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file TribolCoupling.cpp
 */

#include "TribolCoupling.hpp"

namespace geosx
{

void TribolCoupling::Initialize(dataRepository::ManagedGroup * eventManager, dataRepository::ManagedGroup * domain)
{
  real64& dt = *(eventManager->getData<real64>("dt"));
  integer& cycle = *(eventManager->getData<integer>("cycle"));

  // Currently we need the previous dt, but we assume fixed dt for now.
  // We probably want to pass in prevDt in the update anyway.
  s_tribolProblem = new vista::View("geosx", 0) ;
  s_tribolDomain = s_tribolProblem->viewCreate("DomainWorld", 0)->viewCreate("domain", 0, 0) ;
  s_tribolDomain->viewCreate("extFaces", 0) ;
  s_tribolDomain->viewCreate("nodes", 0) ;
  s_tribolDomain->viewCreate("bricks", 0) ;
  SlideWorldAdapter::CreateWorld(MPI_COMM_GEOSX, MPI_COMM_WORLD,
                                 4000, // comm tag 
                                 3, // dimension
                                 0, // axisym
                                 1, // numSS
                                 4, // nodesPerFace
                                 8, // nodesPerElem
                                 &cycle,
                                 &dt,
                                 &dt, // prevDt
                                 1, // noParamInput
                                 s_tribolProblem) ;
  SlideWorldAdapter::SetSourceData(s_tribolProblem) ;
}

void TribolCoupling::SyncTermination(int* terminate)
{
   SlideWorldAdapter::SyncTermination(terminate) ;
}

void TribolCoupling::SyncTimestep(real64* newDt)
{
    // Currently does not change GEOS timestep.
    // Will also be incorrect if timestep is adjusted by maxTime.
    SlideWorldAdapter::GetRequiredTimestep(newDt) ;

    // Ignore GEOS requests
    int plotThisCycle = 0 ;
    int dumpThisCycle = 0 ;
    int advectCountThisCycle = 0 ;
    int slideDeleteThisCycle = 0 ;

    SlideWorldAdapter::SyncControl(&slideDeleteThisCycle, &advectCountThisCycle, &plotThisCycle, &dumpThisCycle) ;

    if (slideDeleteThisCycle) {
        SlideWorldAdapter::Cleanup();
        SlideWorldAdapter::SetSourceData(s_tribolProblem, true);
    }
}

void TribolCoupling::Cleanup()
{
   delete s_tribolProblem ;
   s_tribolDomain = nullptr ;
}

} /* namespace geosx */
