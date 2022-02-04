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

#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include <iostream>

namespace geosx {

using namespace dataRepository;
using namespace constitutive;

template< typename POROUSWRAPPER_TYPE >
void execute4( POROUSWRAPPER_TYPE porousWrapper,
		CellElementSubRegion & subRegion,
               arrayView1d< real64 const > const & volume,
               arrayView1d< real64 const > const & deltaVolume )
{
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < porousWrapper.numGauss(); ++q )
    {
      porousWrapper.updatePorosityFromVolume( k, q,
    		  volume[k],
			  deltaVolume[k] );
    }
  } );
}

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
	SinglePhasePoromechanicsSolver::postProcessInput();
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


void SinglePhasePoromechanicsSolverDynamic::updateDeformationForCoupling( DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  NodeManager const & nodeManager = mesh.getNodeManager();

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const u = nodeManager.totalDisplacement();

  if( m_couplingTypeOption == CouplingTypeOption::FEM_ExplicitTransient
		  || m_couplingTypeOption == CouplingTypeOption::FEM_ImplicitTransient)
  {
    // update matrix volume and porosity
    forTargetSubRegionsComplete< CellElementSubRegion >( mesh, [&]( localIndex const,
                                                                    localIndex const,
                                                                    localIndex const,
                                                                    ElementRegionBase &,
                                                                    CellElementSubRegion & elementSubRegion )
    {
     arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

     ElementType elementType = elementSubRegion.getElementType();

     arrayView1d< real64 > const & volume = elementSubRegion.getReference< array1d< real64 > >( CellElementSubRegion::viewKeyStruct::elementVolumeString() );;

     arrayView1d< real64 > const & deltaVolume = elementSubRegion.getExtrinsicData< extrinsicMeshData::flow::deltaVolume >();

     localIndex const numNodesPerElement = elemsToNodes.size( 1 );

     // TODO: remove use of R1Tensor and use device policy
     forAll< parallelDevicePolicy< 32 > >( elementSubRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
     {
       // update element volume
       real64 Xlocal[10][3];
       for( localIndex a = 0; a < numNodesPerElement; ++a )
       {
    	   for ( localIndex aa = 0; aa < 3; ++aa )
    	   {
    		   Xlocal[a][aa] = nodePosition[elemsToNodes[ei][a]][aa];
    		   Xlocal[a][aa] += u[elemsToNodes[ei][a]][aa];
    	   }
       }
       switch( elementType )
       {
       	   case ElementType::Hexahedron:
       	   {
       		   deltaVolume[ei] = computationalGeometry::hexVolume( Xlocal ) - volume[ei];
               break;
       	   }
       	   case ElementType::Tetrahedron:
       	   {
               deltaVolume[ei] = computationalGeometry::tetVolume( Xlocal ) - volume[ei];
               break;
           }
           case ElementType::Prism:
           {
        	   deltaVolume[ei] = computationalGeometry::wedgeVolume( Xlocal ) - volume[ei];
               break;
           }
           case ElementType::Pyramid:
           {
        	   deltaVolume[ei] = computationalGeometry::pyramidVolume( Xlocal ) - volume[ei];
               break;
           }
           default:
           {
               GEOSX_ERROR( "Volume calculation not supported for element type " << elementType << " and for elementSubRegion " );
           }
       }

       volume[ei] += deltaVolume[ei];
     } );
   } );


    forTargetSubRegions< CellElementSubRegion >( mesh, [&]( localIndex const targetIndex,
    		CellElementSubRegion & subRegion )
        {

    	  arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
    	  arrayView1d< real64 const > const & deltaVolume = subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaVolume >();

    	  CoupledSolidBase & porousSolid = subRegion.template getConstitutiveModel< CoupledSolidBase >( m_flowSolver->solidModelNames()[targetIndex] );

    	  constitutive::ConstitutivePassThru< CoupledSolidBase >::execute( porousSolid, [=, &subRegion] ( auto & castedPorousSolid )
    	  {
    	    typename TYPEOFREF( castedPorousSolid ) ::KernelWrapper porousWrapper = castedPorousSolid.createKernelUpdates();

    	    execute4( porousWrapper, subRegion, volume, deltaVolume );
    	  } );

        });
  }  // end of if

}


real64 SinglePhasePoromechanicsSolverDynamic::explicitStep( real64 const & time_n,
                                          real64 const & dt,
                                          const int cycleNumber,
                                          DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  if( m_couplingTypeOption == CouplingTypeOption::FEM_ExplicitTransient )
  {
	  //m_flowSolver->calculateAndApplyMassFlux( time_n, dt, domain );
	  m_solidSolver->explicitStepDisplacementUpdate( time_n, dt, cycleNumber, domain );
	  this->updateDeformationForCoupling( domain );
	  //m_flowSolver->updateEOSExplicit( time_n, dt, domain );
	  m_solidSolver->explicitStepVelocityUpdate( time_n, dt, cycleNumber, domain );
  }
  else if( m_couplingTypeOption == CouplingTypeOption::FEM_ImplicitTransient )
  {
	  m_solidSolver->explicitStep( time_n, dt, cycleNumber, domain );
	  //Apply deformation to flowsolver
	  this->updateDeformationForCoupling( domain );
	  m_flowSolver->solverStep( time_n, dt, cycleNumber, domain );
  }
  return dt;
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhasePoromechanicsSolverDynamic, string const &, Group * const )

} /* namespace geosx */
