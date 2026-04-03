/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PhaseFieldPoromechanicsSolver.cpp
 *
 */

#include "PhaseFieldPoromechanicsSolver.hpp"

#include "fieldSpecification/TractionBoundaryCondition.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

template< typename FLOW_SOLVER >
PhaseFieldPoromechanicsSolver< FLOW_SOLVER >::PhaseFieldPoromechanicsSolver( const string & name,
                                                                             Group * const parent ):
  Base( name, parent )
{}

template< typename FLOW_SOLVER >
void PhaseFieldPoromechanicsSolver< FLOW_SOLVER >::resetStateToBeginningOfStep( DomainPartition & domain )
{
  Base::resetStateToBeginningOfStep( domain );

  applyDamageOnTractionBC( domain );
}

template< typename FLOW_SOLVER >
void PhaseFieldPoromechanicsSolver< FLOW_SOLVER >::postInputInitialization()
{
  Base::postInputInitialization();
  GEOS_WARNING_IF( this->getNonlinearSolverParameters().couplingType() == NonlinearSolverParameters::CouplingType::FullyImplicit,
                   "FullyImplicit coupling not implemented for this solver. A sequential coupling approach will be used." );
  this->getNonlinearSolverParameters().m_couplingType = NonlinearSolverParameters::CouplingType::Sequential;
}

template< typename FLOW_SOLVER >
void PhaseFieldPoromechanicsSolver< FLOW_SOLVER >::mapSolutionBetweenSolvers( real64 const & dt, DomainPartition & domain, integer const solverType )
{
  GEOS_UNUSED_VAR( dt );

  if( solverType ==  static_cast< integer >( SolverType::Damage ) )
  {
    GEOS_MARK_FUNCTION;
    this->template forDiscretizationOnMeshTargets<>( domain.getMeshBodies(), [&] ( string const &,
                                                                                   MeshLevel & mesh,
                                                                                   string_array const & regionNames )
    {
      NodeManager & nodeManager = mesh.getNodeManager();

      arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const xNodes = nodeManager.referencePosition();

      string const & damageFieldName = this->damageSolver()->getFieldName();

      string const & discretizationName = this->damageSolver()->getDiscretizationName();

      //should get reference to damage field here.
      arrayView1d< real64 const > const nodalDamage = nodeManager.getReference< array1d< real64 > >( damageFieldName );

      ElementRegionManager & elemManager = mesh.getElemManager();

      // begin region loop
      elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [discretizationName, xNodes, nodalDamage]
                                                                  ( localIndex const,
                                                                  CellElementSubRegion & elementSubRegion )
      {
        string const & solidModelName = elementSubRegion.getReference< string >( SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString());
        constitutive::SolidBase &
        solidModel = elementSubRegion.getConstitutiveModel< constitutive::SolidBase >( solidModelName );

        ConstitutivePassThru< DamageBase >::execute( solidModel, [&elementSubRegion, discretizationName, xNodes, nodalDamage]( auto & damageModel )
        {
          using CONSTITUTIVE_TYPE = TYPEOFREF( damageModel );
          typename CONSTITUTIVE_TYPE::KernelWrapper constitutiveUpdate = damageModel.createKernelUpdates();

          arrayView2d< real64 > const damageFieldOnMaterial = constitutiveUpdate.m_newDamage;
          arrayView3d< real64 > const damageGradOnMaterial = constitutiveUpdate.m_damageGrad;
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemToNodes = elementSubRegion.nodeList();

          finiteElement::FiniteElementBase const &
          fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( discretizationName );

          finiteElement::FiniteElementDispatchHandler< ALL_FE_TYPES >::dispatch3D( fe, [=, &elementSubRegion] ( auto & finiteElement )
          {
            using FE_TYPE = TYPEOFREF( finiteElement );

            DamageAndDamageGradientInterpolationKernel< FE_TYPE > interpolationKernel( elementSubRegion );

            interpolationKernel.interpolateDamageAndGradient( elemToNodes, xNodes, nodalDamage, damageFieldOnMaterial, damageGradOnMaterial );
          } );
        } );
      } );
    } );
  }
  else if( solverType == static_cast< integer >( SolverType::Poromechanics ) )
  {
    this->poromechancisSolver()->flowSolver()->updatePressureGradient( domain );
  }
}

template< typename FLOW_SOLVER >
void PhaseFieldPoromechanicsSolver< FLOW_SOLVER >::applyDamageOnTractionBC( DomainPartition & domain )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  GEOS_MARK_FUNCTION;
  this->template forDiscretizationOnMeshTargets<>( domain.getMeshBodies(), [&] ( string const &,
                                                                                 MeshLevel & mesh,
                                                                                 string_array const & )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    string const & damageFieldName = this->damageSolver()->getFieldName();

    // Get an array of nodal damage values
    arrayView1d< real64 const > const nodalDamage = nodeManager.getReference< array1d< real64 > >( damageFieldName );

    fsManager.forSubGroups< TractionBoundaryCondition >( [&] ( TractionBoundaryCondition & fs )
    {
      stdVector< string > const targetPath = stringutilities::tokenize( fs.getObjectPath(), "/" );

      dataRepository::Group * targetGroup = &mesh;

      dataRepository::Group * const elemRegionSubGroup = targetGroup->getGroupPointer( ElementRegionManager::groupKeyStruct::elementRegionsGroup() );

      if( elemRegionSubGroup != nullptr )
      {
        targetGroup = elemRegionSubGroup;
      }

      dataRepository::Group * const elemSubRegionSubGroup = targetGroup->getGroupPointer( ElementRegionBase::viewKeyStruct::elementSubRegions() );
      if( elemSubRegionSubGroup != nullptr )
      {
        targetGroup = elemSubRegionSubGroup;
      }

      targetGroup = &targetGroup->getGroup( targetPath[0] );

      Group & target = *targetGroup;

      dataRepository::Group const & setGroup = target.getGroup( ObjectManagerBase::groupKeyStruct::setsString() );
      string_array setNames = fs.getSetNames();
      for( auto & setName : setNames )
      {
        if( setGroup.hasWrapper( setName ) )
        {
          SortedArrayView< localIndex const > const & targetSet = setGroup.getReference< SortedArray< localIndex > >( setName );

          fs.reinitScaleSet( faceManager,
                             targetSet,
                             nodalDamage );
        }
      }
    } );
  } );

}

template class PhaseFieldPoromechanicsSolver<>;
template class PhaseFieldPoromechanicsSolver< SinglePhaseReactiveTransport >;

namespace
{
typedef PhaseFieldPoromechanicsSolver< SinglePhaseReactiveTransport > PhaseFieldReactivePoromechanics;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, PhaseFieldReactivePoromechanics, string const &, Group * const )
typedef PhaseFieldPoromechanicsSolver<> PhaseFieldPoromechanics;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, PhaseFieldPoromechanics, string const &, Group * const )
}

} /* namespace geos */
