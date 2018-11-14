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
 * LinearSolverWrapper.cpp
 *
 *  Created on: Sep 12, 2017
 *      Author: settgast
 */

#include "LinearSolverWrapper.hpp"


#include "Epetra_Map.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_LinearProblem.h"

#include "EpetraExt_SolverMap_CrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

#include "AztecOO.h"

#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"
#include "ml_RowMatrix.h"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

#include "Amesos_BaseSolver.h"
#include "Amesos.h"
#include "Thyra_OperatorVectorClientSupport.hpp"
#include "Thyra_AztecOOLinearOpWithSolveFactory.hpp"
#include "Thyra_AztecOOLinearOpWithSolve.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraLinearOpBase.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_DefaultInverseLinearOp.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_MLPreconditionerFactory.hpp"


#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "linearAlgebraInterface/src/TrilinosInterface.hpp"

namespace geosx
{
namespace systemSolverInterface
{

LinearSolverWrapper::LinearSolverWrapper():
#ifdef GEOSX_USE_MPI
  m_epetraComm(MPI_COMM_GEOSX)
#else
  m_epetraComm()
#endif
{
  // TODO Auto-generated constructor stub

}

LinearSolverWrapper::~LinearSolverWrapper()
{
  // TODO Auto-generated destructor stub
}





}
} /* namespace geosx */
