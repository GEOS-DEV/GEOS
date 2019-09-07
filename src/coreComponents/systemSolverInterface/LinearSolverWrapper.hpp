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
* LinearSolverWrapper.hpp
*
*  Created on: Sep 12, 2017
*      Author: settgast
*/

#ifndef SRC_COMPONENTS_CORE_SRC_SYSTEMSOLVERINTERFACE_LINEARSOLVERWRAPPER_HPP_
#define SRC_COMPONENTS_CORE_SRC_SYSTEMSOLVERINTERFACE_LINEARSOLVERWRAPPER_HPP_

#include "physicsSolvers/SystemSolverParameters.hpp"
#include "systemSolverInterface/EpetraBlockSystem.hpp"

namespace geosx
{
namespace systemSolverInterface
{

class LinearSolverWrapper
{
public:
LinearSolverWrapper();
virtual ~LinearSolverWrapper();

void SolveSingleBlockSystem( EpetraBlockSystem * const system,
SystemSolverParameters const * const params,
BlockIDs const blockID );

#ifdef GEOSX_USE_MPI
Epetra_MpiComm m_epetraComm;
#else
Epetra_SerialComm m_epetraComm;
#endif

};



}
} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_SYSTEMSOLVERINTERFACE_LINEARSOLVERWRAPPER_HPP_
*/
