/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

/*
 * Communication.h
 *
 *  Created on: Aug 24, 2012
 *      Author: johnson346
 */

#ifndef COMMUNICATION_H_
#define COMMUNICATION_H_

#include <mpi.h>

namespace CommRegistry
{
enum commID
{
  genericComm01,
  genericComm02,
  lagrangeSolver01,
  lagrangeSolver02,
  parallelPlateFlowSolver,
  parallelPlateFlowSolverClosedJoint,
  ReactionFrontSolver,
  //reactionFrontApertureUpdate, // removed until integrated into trunk and
  // tests updated.
  steadyStateParallelPlateFlowSolver,
  steadyStateParallelPlateFlowSolverB,
  immiscibleFluidSolver,
  twoDARDSolver,
  lagrangeParallelPlateFlowSolver,
  implicitLaplaceComm,
  bemMechanicsSolver,
  hydrofractureSolver,
  fractureFlowSolver,
  matrixFlowSolver,
  matrixFlowSolverFV,
  singlePhasePMFlowThermalSolverFV,
  processShareSendSize,
  processShareSendNodeData,
  processShareSendFaceData,
  processShareSendFacePairData,
  processShareReturn,

  maxComm
};

inline int CommTag( const int senderRank, const int receiverRank, const commID comm )
{
  int m_size;
  MPI_Comm_size( MPI_COMM_GEOSX, &m_size );
  return senderRank * m_size + receiverRank + m_size*m_size*comm;
}
}



class Communication
{
public:
  Communication();
  virtual ~Communication();
};

#endif /* COMMUNICATION_H_ */
