/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhasePoromechanicsSolverDynamic.cpp
 *
 *  Created on: Jan 29, 2022
 *      Author: ron
 */

#include "SinglePhasePoromechanicsSolverDynamic.hpp"

#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"

namespace geosx {

using namespace dataRepository;

SinglePhasePoromechanicsSolverDynamic::SinglePhasePoromechanicsSolverDynamic( const string & name,
        Group * const parent ):
SinglePhasePoromechanicsSolver( name, parent )
{
	// TODO Auto-generated constructor stub
	registerWrapper( viewKeyStruct::couplingTypeOptionStringString(), &m_couplingTypeOption ).
	    setInputFlag( InputFlags::REQUIRED ).
	    setDescription( "Coupling method. Valid options:\n* " + EnumStrings< CouplingTypeOption >::concat( "\n* " ) );
}

SinglePhasePoromechanicsSolverDynamic::~SinglePhasePoromechanicsSolverDynamic()
{
	// TODO Auto-generated destructor stub
}

void SinglePhasePoromechanicsSolverDynamic::postProcessInput()
{
}

void SinglePhasePoromechanicsSolverDynamic::initializePostInitialConditionsPreSubGroups()
{}


real64 SinglePhasePoromechanicsSolverDynamic::solverStep( real64 const & time_n,
                                        real64 const & dt,
                                        int const cycleNumber,
                                        DomainPartition & domain )
{
  real64 dtReturn = dt;

  if( m_couplingTypeOption == CouplingTypeOption::FIM )
  {
	  dtReturn = SinglePhasePoromechanicsSolver::solverStep(time_n, dt, cycleNumber, domain);

  }
  else
  {
	  dtReturn = explicitStep( time_n, dt, cycleNumber, domain );
  }

  return dtReturn;
}


real64 SinglePhasePoromechanicsSolverDynamic::explicitStep( real64 const & time_n,
                                          real64 const & dt,
                                          const int cycleNumber,
                                          DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  if( m_couplingTypeOption == CouplingTypeOption::FEM_ExplicitTransient )
  {
	  m_solidSolver->explicitStep( time_n, dt, cycleNumber, domain );
  }
  else if( m_couplingTypeOption == CouplingTypeOption::FEM_ImplicitTransient )
  {
	  m_solidSolver->explicitStep( time_n, dt, cycleNumber, domain );
	  //Apply deformation to flowsolver
	  //this->UpdateDeformationForCoupling( domain );
	  m_flowSolver->solverStep( time_n, dt, cycleNumber, domain );
  }
  return dt;
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhasePoromechanicsSolverDynamic, string const &, Group * const )

} /* namespace geosx */
