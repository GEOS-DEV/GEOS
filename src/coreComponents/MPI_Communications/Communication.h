/*
* ------------------------------------------------------------------------------------------------------------
* SPDX-License-Identifier: LGPL-2.1-only
*
* Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
* Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
* Copyright (c) 2018-2019 Total, S.A
* Copyright (c) 2019-     GEOSX Contributors
* All right reserved
*
* See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
* ------------------------------------------------------------------------------------------------------------
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
