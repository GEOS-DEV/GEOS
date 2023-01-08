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
 * @file ThermalSinglePhasePoromechanicsFixedStressSolver.cpp
 *
 */

#include "ThermalSinglePhasePoromechanicsFixedStressSolver.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mesh/DomainPartition.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

ThermalSinglePhasePoromechanicsFixedStressSolver::ThermalSinglePhasePoromechanicsFixedStressSolver( const string & name,
                                                                                                    Group * const parent ):
  Base( name, parent )
{}

void ThermalSinglePhasePoromechanicsFixedStressSolver::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  // Set m_updateThermalPoroElasticityFlag in SolidMechanicsLagrangianFEM to 1
  solidMechanicsSolver()->turnOnFixedStressThermoPoroElasticityFlag();
}

void ThermalSinglePhasePoromechanicsFixedStressSolver::postProcessInput()
{
  Base::postProcessInput();
  GEOSX_WARNING_IF( getNonlinearSolverParameters().m_couplingType == NonlinearSolverParameters::CouplingType::FullyImplicit,
                    "FullyImplicit coupling not implemented for this solver. A sequential coupling approach will be used." );
  getNonlinearSolverParameters().m_couplingType = NonlinearSolverParameters::CouplingType::Sequential;
}

ThermalSinglePhasePoromechanicsFixedStressSolver::~ThermalSinglePhasePoromechanicsFixedStressSolver()
{
  // TODO Auto-generated destructor stub
}

void ThermalSinglePhasePoromechanicsFixedStressSolver::mapSolutionBetweenSolvers( DomainPartition & domain, integer const solverType )
{
  GEOSX_MARK_FUNCTION;
  if( solverType == static_cast< integer >( SolverType::SolidMechanics ) )
  {
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            auto & subRegion )
      {
        flowSolver()->updatePorosityAndPermeability( subRegion );
      } );
    } );
  }
}

REGISTER_CATALOG_ENTRY( SolverBase, ThermalSinglePhasePoromechanicsFixedStressSolver, string const &, Group * const )

} /* namespace geosx */
