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
 * @file SolverApplicationSet.h
 * @author settgast1
 * @date Mar 2, 2011
 */

#ifndef SOLVERAPPLICATIONSET_H_
#define SOLVERAPPLICATIONSET_H_

#include "Common/Common.h"

typedef std::pair< std::string, array<string> > SolverApplication;

class SolverApplicationSet
{
public:
  SolverApplicationSet();
  ~SolverApplicationSet();

  std::vector< SolverApplication > m_solverAppliedToRegion;
  realT m_beginTime;
  realT m_endTime;
  realT m_deltaTime;

};

#endif /* SOLVERAPPLICATIONSET_H_ */
