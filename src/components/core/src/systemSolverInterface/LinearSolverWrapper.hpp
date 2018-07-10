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

/*
 * LinearSolverWrapper.hpp
 *
 *  Created on: Sep 12, 2017
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_SYSTEMSOLVERINTERFACE_LINEARSOLVERWRAPPER_HPP_
#define SRC_COMPONENTS_CORE_SRC_SYSTEMSOLVERINTERFACE_LINEARSOLVERWRAPPER_HPP_

#include "systemSolverInterface/SystemSolverParameters.hpp"
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
                               EpetraBlockSystem::BlockIDs const blockID );

#if USE_MPI
  Epetra_MpiComm m_epetraComm;
#else
  Epetra_SerialComm m_epetraComm;
#endif

};



}
} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_SYSTEMSOLVERINTERFACE_LINEARSOLVERWRAPPER_HPP_
        */
