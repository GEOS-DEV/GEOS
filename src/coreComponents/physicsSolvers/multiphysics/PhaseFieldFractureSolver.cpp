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
 * @file PhaseFieldFractureSolver.cpp
 *
 */

#include "PhaseFieldFractureSolver.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "finiteElement/Kinematics.h"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

PhaseFieldFractureSolver::PhaseFieldFractureSolver( const string & name,
                                                    Group * const parent ):
  Base( name, parent )
{
  // Only sequential coupling implemented for this solver 
  getNonlinearSolverParameters().m_couplingType = NonlinearSolverParameters::CouplingType::Sequential;
}

PhaseFieldFractureSolver::~PhaseFieldFractureSolver()
{
  // TODO Auto-generated destructor stub
}

void PhaseFieldFractureSolver::mapSolutionBetweenSolvers( DomainPartition & domain, integer const idx )
{

  GEOSX_MARK_FUNCTION;
  if( idx ==  static_cast< integer >( SolverType::Damage ) )
  {
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                  MeshLevel & mesh,
                                                                  arrayView1d< string const > const & regionNames )
    {
      NodeManager & nodeManager = mesh.getNodeManager();

      string const & damageFieldName = damageSolver()->getFieldName();

      string const & discretizationName = damageSolver()->getDiscretizationName();

      //should get reference to damage field here.
      arrayView1d< real64 const > const nodalDamage = nodeManager.getReference< array1d< real64 > >( damageFieldName );

      ElementRegionManager & elemManager = mesh.getElemManager();

      // begin region loop
      elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [discretizationName, nodalDamage]
                                                                  ( localIndex const,
                                                                  CellElementSubRegion & elementSubRegion )
      {
        string const & solidModelName = elementSubRegion.getReference< string >( SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString());
        constitutive::SolidBase &
        solidModel = elementSubRegion.getConstitutiveModel< constitutive::SolidBase >( solidModelName );

        ConstitutivePassThru< DamageBase >::execute( solidModel, [&elementSubRegion, discretizationName, nodalDamage]( auto & damageModel )
        {
          using CONSTITUTIVE_TYPE = TYPEOFREF( damageModel );
          typename CONSTITUTIVE_TYPE::KernelWrapper constitutiveUpdate = damageModel.createKernelUpdates();

          arrayView2d< real64 > const damageFieldOnMaterial = constitutiveUpdate.m_damage;
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemNodes = elementSubRegion.nodeList();

          finiteElement::FiniteElementBase const &
          fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( discretizationName );

          integer const numElems = elementSubRegion.size();

          finiteElement::FiniteElementDispatchHandler< ALL_FE_TYPES >::dispatch3D( fe, [=] ( auto & finiteElement )
          {
            using FE_TYPE = TYPEOFREF( finiteElement );           

            forAll< parallelDevicePolicy<> >( numElems, [=] ( localIndex const k )
            {
              constexpr localIndex numNodesPerElement = FE_TYPE::numNodes;
              constexpr localIndex n_q_points = FE_TYPE::numQuadraturePoints;

              for( localIndex q = 0; q < n_q_points; ++q )
              {
                real64 N[ numNodesPerElement ];
                FE_TYPE::calcN( q, N );

                damageFieldOnMaterial( k, q ) = 0;
                for( localIndex a = 0; a < numNodesPerElement; ++a )
                {
                  damageFieldOnMaterial( k, q ) += N[a] * nodalDamage[elemNodes( k, a )];
                  //solution is probably not going to work because the solution of the coupled solver
                  //has both damage and displacements. Using the damageResult field from the Damage solver
                  //is probably better
                }
              }
            } );
          } );
        } );
      } );
    } );
  }

}

REGISTER_CATALOG_ENTRY( SolverBase, PhaseFieldFractureSolver, string const &, Group * const )

} /* namespace geosx */
