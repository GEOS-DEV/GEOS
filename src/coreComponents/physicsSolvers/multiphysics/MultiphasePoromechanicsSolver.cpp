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
 * @file MultiphasePoroelasticSolver.cpp
 */

#include "MultiphasePoromechanicsSolver.hpp"

#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/solid/PorousSolid.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/multiphysics/MultiphasePoromechanicsKernel.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

MultiphasePoromechanicsSolver::MultiphasePoromechanicsSolver( const string & name,
                                                              Group * const parent )
  : Base( name, parent )
{
  registerWrapper( viewKeyStruct::stabilizationTypeString(), &m_stabilizationType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Stabilization type. Options are:\n" +
                    toString( StabilizationType::None ) + " - Add no stabilization to mass equation,\n" +
                    toString( StabilizationType::Global ) + " - Add stabilization to all faces,\n" +
                    toString( StabilizationType::Local ) + " - Add stabilization only to interiors of macro elements." );

  registerWrapper( viewKeyStruct::stabilizationRegionNamesString(), &m_stabilizationRegionNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Regions where stabilization is applied." );

  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::multiphasePoromechanics;
  m_linearSolverParameters.get().mgr.separateComponents = true;
  m_linearSolverParameters.get().mgr.displacementFieldName = keys::TotalDisplacement;
  m_linearSolverParameters.get().dofsPerNode = 3;
}

void MultiphasePoromechanicsSolver::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      subRegion.registerWrapper< string >( viewKeyStruct::porousMaterialNamesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setSizedFromParent( 0 );

      subRegion.registerExtrinsicData< extrinsicMeshData::flow::elementMacroID >( getName());
    } );
  } );
}

void MultiphasePoromechanicsSolver::initializePostInitialConditionsPreSubGroups()
{

  SolverBase::initializePostInitialConditionsPreSubGroups();

  if( m_stabilizationType == StabilizationType::Global )
  {
    SortedArray< localIndex > regionFilter;

    DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                  MeshLevel & mesh,
                                                                  arrayView1d< string const > const & )
    {
      ElementRegionManager & elemManager = mesh.getElemManager();
      for( string const & regionName : m_stabilizationRegionNames )
      {
        regionFilter.insert( elemManager.getRegions().getIndex( regionName ) );
      }
    } );

    SortedArrayView< const localIndex > const regionFilterView = regionFilter.toView();

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                  MeshLevel & mesh,
                                                                  arrayView1d< string const > const & )
    {
      ElementRegionManager & elemManager = mesh.getElemManager();
      NodeManager const & nodeManager = mesh.getNodeManager();

      ElementRegionManager::ElementViewAccessor< arrayView1d< integer > > elemMacroID =
        elemManager.constructViewAccessor< array1d< integer >, arrayView1d< integer > >( extrinsicMeshData::flow::elementMacroID::key() );

      ArrayOfArraysView< localIndex const > elemRegionList = nodeManager.elementRegionList();
      ArrayOfArraysView< localIndex const > elemSubRegionList = nodeManager.elementSubRegionList();
      ArrayOfArraysView< localIndex const > elemList = nodeManager.elementList();

      forAll< serialPolicy >( nodeManager.size(), [&] ( localIndex const a )
      {

        for( localIndex k = 0; k < elemRegionList[a].size(); ++k )
        {
          // collect the element number
          localIndex const er = elemRegionList[a][k];
          localIndex const esr = elemSubRegionList[a][k];
          localIndex const ei = elemList[a][k];

          if( regionFilterView.contains( er ))
          {
            elemMacroID[er][esr][ei] = 1;
          }
        }
      } );

    } );


  }

  if( m_stabilizationType == StabilizationType::Local )
  {

    SortedArray< localIndex > regionFilter;

    DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                  MeshLevel & mesh,
                                                                  arrayView1d< string const > const & )
    {
      ElementRegionManager & elemManager = mesh.getElemManager();
      for( string const & regionName : m_stabilizationRegionNames )
      {
        regionFilter.insert( elemManager.getRegions().getIndex( regionName ) );
      }
    } );

    SortedArrayView< const localIndex > const regionFilterView = regionFilter.toView();

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                  MeshLevel & mesh,
                                                                  arrayView1d< string const > const & regionNames )
    {

      ElementRegionManager & elemManager = mesh.getElemManager();

      GEOSX_UNUSED_VAR( regionNames );

      NodeManager const & nodeManager = mesh.getNodeManager();

      ElementRegionManager::ElementViewAccessor< arrayView1d< integer > > elemMacroID =
        elemManager.constructViewAccessor< array1d< integer >, arrayView1d< integer > >( extrinsicMeshData::flow::elementMacroID::key() );

      array1d< integer > nodeVisited( nodeManager.size() );
      nodeVisited.setValues< serialPolicy >( 0 );
      arrayView1d< integer > const nodeVisitedView = nodeVisited.toView();

      arrayView1d< integer const > const bdryNodes = nodeManager.getDomainBoundaryIndicator();

      ArrayOfArraysView< localIndex const > elemRegionList = nodeManager.elementRegionList();
      ArrayOfArraysView< localIndex const > elemSubRegionList = nodeManager.elementSubRegionList();
      ArrayOfArraysView< localIndex const > elemList = nodeManager.elementList();

      array1d< integer > mpiBdryNodes( nodeManager.size() );
      mpiBdryNodes.setValues< serialPolicy >( 0 );
      arrayView1d< integer > const mpiBdryNodesView = mpiBdryNodes.toView();
      ElementRegionManager::ElementViewAccessor< arrayView1d< integer > > elemGhostRank =
        elemManager.constructViewAccessor< array1d< integer >, arrayView1d< integer > >( extrinsicMeshData::ghostRank::key() );

      forAll< serialPolicy >( nodeManager.size(), [&] ( localIndex const a )
      {
        for( localIndex k = 0; k < elemRegionList[a].size(); ++k )
        {
          localIndex const er = elemRegionList[a][k];
          localIndex const esr = elemSubRegionList[a][k];
          localIndex const ei = elemList[a][k];

          if( elemGhostRank[er][esr][ei] >= 0 )
          {
            mpiBdryNodesView[a] = 1;
          }
        }
      } );

      integer currentID = 0;

      forAll< serialPolicy >( nodeManager.size(), [&] ( localIndex const a )
      {

        if( bdryNodes[a] == 1 )
        {
          // nodeVisitedView[a] = 1;
        }

        if( mpiBdryNodes[a] == 1 )
        {
          nodeVisitedView[a] = 1;
        }

        if( nodeVisitedView[a] != 1 )
        {

          for( localIndex k = 0; k < elemRegionList[a].size(); ++k )
          {
            // collect the element number
            localIndex const er = elemRegionList[a][k];
            localIndex const esr = elemSubRegionList[a][k];
            localIndex const ei = elemList[a][k];

            if( regionFilterView.contains( er ))
            {
              elemMacroID[er][esr][ei] = currentID;
            }

            // get the elemToNodes maps
            ElementRegionBase const & region = elemManager.getRegion( er );
            CellElementSubRegion const & subRegion = region.getSubRegion< CellElementSubRegion, localIndex >( esr );
            arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = subRegion.nodeList();

            // get the nodes connected to this element
            for( localIndex l = 0; l < elemsToNodes[ei].size(); ++l )
            {
              localIndex const iNode = elemsToNodes[ei][l];
              nodeVisitedView[iNode] = 1;

            }
          }

          ++currentID;
        }

      } );


      // part 2 - assign any unassigned cell
      // loop through all elements and check if each has been assigned a macroelement (ie, when you index into ID it should be >0)
      // If not, check elements that share a face and add it to one of their macroelements
      // If no neighbors in a macroelement, reconsider this element at the end of the loop


      FaceManager const & faceManager = mesh.getFaceManager();

      arrayView2d< localIndex const > elemRegionListFace = faceManager.elementRegionList();
      arrayView2d< localIndex const > elemSubRegionListFace = faceManager.elementSubRegionList();
      arrayView2d< localIndex const > elemListFace = faceManager.elementList();

      integer badFaces = 1;

      while( badFaces > 0 )
      {
        badFaces = 0;
        forAll< serialPolicy >( faceManager.size(), [&] ( localIndex const a )
        {
          localIndex er1 = elemRegionListFace[a][0];
          localIndex esr1 = elemSubRegionListFace[a][0];
          localIndex ei1 = elemListFace[a][0];

          localIndex er2 = elemRegionListFace[a][1];
          localIndex esr2 = elemSubRegionListFace[a][1];
          localIndex ei2 = elemListFace[a][1];

          if( er1 < 0 || esr1 < 0 || ei1 < 0 || er2 < 0 || esr2 < 0 || ei2 < 0 )
          {}
          else
          {

            if( elemMacroID[er1][esr1][ei1] == -1 && elemMacroID[er2][esr2][ei2] > -1 && elemGhostRank[er1][esr1][ei1] < 0 && regionFilterView.contains( er1 ) &&  regionFilterView.contains( er2 ) )
            {
              elemMacroID[er1][esr1][ei1] = elemMacroID[er2][esr2][ei2];
            }

            if( elemMacroID[er1][esr1][ei1] > -1 && elemMacroID[er2][esr2][ei2] == -1 && elemGhostRank[er2][esr2][ei2] < 0 && regionFilterView.contains( er1 ) &&  regionFilterView.contains( er2 ) )
            {
              elemMacroID[er2][esr2][ei2] = elemMacroID[er1][esr1][ei1];
            }

            if( elemMacroID[er1][esr1][ei1] == -1 && elemMacroID[er2][esr2][ei2] == -1 && elemGhostRank[er1][esr1][ei1] < 0 && elemGhostRank[er2][esr2][ei2] < 0 && regionFilterView.contains( er1 ) &&
                regionFilterView.contains( er2 ))
            {
              ++badFaces;
            }

          }

        } );

        std::cout << badFaces << std::endl;
      }

    } );

  }

}

void MultiphasePoromechanicsSolver::setupCoupling( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                                   DofManager & dofManager ) const
{
  dofManager.addCoupling( keys::TotalDisplacement,
                          CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem );
}

void MultiphasePoromechanicsSolver::assembleSystem( real64 const GEOSX_UNUSED_PARAM( time ),
                                                    real64 const dt,
                                                    DomainPartition & domain,
                                                    DofManager const & dofManager,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();

    string const displacementDofKey = dofManager.getKey( dataRepository::keys::TotalDisplacement );
    arrayView1d< globalIndex const > const & displacementDofNumber = nodeManager.getReference< globalIndex_array >( displacementDofKey );

    string const flowDofKey = dofManager.getKey( CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString() );

    localIndex const numComponents = flowSolver()->numFluidComponents();
    localIndex const numPhases = flowSolver()->numFluidPhases();

    real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

    poromechanicsKernels::MultiphaseKernelFactory kernelFactory( displacementDofNumber,
                                                                 flowDofKey,
                                                                 dofManager.rankOffset(),
                                                                 gravityVectorData,
                                                                 numComponents,
                                                                 numPhases,
                                                                 FlowSolverBase::viewKeyStruct::fluidNamesString(),
                                                                 localMatrix,
                                                                 localRhs );

    // Cell-based contributions
    solidMechanicsSolver()->getMaxForce() =
      finiteElement::
        regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                      constitutive::PorousSolidBase,
                                      CellElementSubRegion >( mesh,
                                                              regionNames,
                                                              solidMechanicsSolver()->getDiscretizationName(),
                                                              viewKeyStruct::porousMaterialNamesString(),
                                                              kernelFactory );
  } );

  // Face-based contributions
  if( m_stabilizationType == StabilizationType::Global ||  m_stabilizationType == StabilizationType::Local )
  {
    flowSolver()->assembleStabilizedFluxTerms( dt,
                                               domain,
                                               dofManager,
                                               localMatrix,
                                               localRhs );
  }
  else
  {
    flowSolver()->assembleFluxTerms( dt,
                                     domain,
                                     dofManager,
                                     localMatrix,
                                     localRhs );
  }
}

real64 MultiphasePoromechanicsSolver::solverStep( real64 const & time_n,
                                                  real64 const & dt,
                                                  int const cycleNumber,
                                                  DomainPartition & domain )
{
  real64 dt_return = dt;

  // setup monolithic coupled system
  SolverBase::setupSystem( domain, m_dofManager, m_localMatrix, m_rhs, m_solution );

  implicitStepSetup( time_n, dt, domain );

  dt_return = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  implicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

void MultiphasePoromechanicsSolver::updateState( DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   CellElementSubRegion & subRegion )
    {
      flowSolver()->updateFluidState( subRegion );
    } );
  } );
}

void MultiphasePoromechanicsSolver::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elementRegionManager = mesh.getElemManager();
    elementRegionManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                                       [&]( localIndex const,
                                                                            ElementSubRegionBase & subRegion )
    {
      string & porousName = subRegion.getReference< string >( viewKeyStruct::porousMaterialNamesString() );
      porousName = getConstitutiveName< CoupledSolidBase >( subRegion );
      GEOSX_ERROR_IF( porousName.empty(), GEOSX_FMT( "Solid model not found on subregion {}", subRegion.getName() ) );
    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, MultiphasePoromechanicsSolver, string const &, Group * const )

} /* namespace geosx */
