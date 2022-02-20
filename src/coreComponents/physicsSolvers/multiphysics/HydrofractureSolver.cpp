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

/**
 * @file HydrofractureSolver.cpp
 *
 */


#include "HydrofractureSolver.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactRelationBase.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "finiteElement/Kinematics.h"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/FaceElementRegion.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

//TJ: access tip quantities in the SurfaceGenerator class
#include "physicsSolvers/surfaceGeneration/SurfaceGenerator.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "managers/FieldSpecification/SourceFluxBoundaryCondition.hpp"

//TJ: used to write (time, tipLoc) pair into a file
#include <iostream>
#include <fstream>

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

HydrofractureSolver::HydrofractureSolver( const std::string & name,
                                          Group * const parent ):
  SolverBase( name, parent ),
  m_solidSolverName(),
  m_flowSolverName(),
  m_couplingTypeOptionString( "FIM" ),
  m_couplingTypeOption(),
  m_solidSolver( nullptr ),
  m_flowSolver( nullptr ),
  m_maxNumResolves( 10 )
{
  registerWrapper( viewKeyStruct::solidSolverNameString, &m_solidSolverName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the solid mechanics solver to use in the poroelastic solver" );

  registerWrapper( viewKeyStruct::fluidSolverNameString, &m_flowSolverName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the fluid mechanics solver to use in the poroelastic solver" );

  registerWrapper( viewKeyStruct::couplingTypeOptionStringString, &m_couplingTypeOptionString )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Coupling option: (FIM, SIM_FixedStress)" );

  registerWrapper( viewKeyStruct::contactRelationNameString, &m_contactRelationName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of contact relation to enforce constraints on fracture boundary." );

  registerWrapper( viewKeyStruct::maxNumResolvesString, &m_maxNumResolves )->
    setApplyDefaultValue( 10 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Value to indicate how many resolves may be executed to perform surface generation after the execution of flow and mechanics solver. " );

  m_numResolves[0] = 0;
}

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
void HydrofractureSolver::RegisterDataOnMesh( dataRepository::Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * meshLevel = Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    ElementRegionManager * const elemManager = meshLevel->getElemManager();
    elemManager->forElementRegions< FaceElementRegion >( [&] ( FaceElementRegion * const region )
    {
      region->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion * const subRegion )
      {
        subRegion->registerWrapper< array1d< real64 > >( viewKeyStruct::separationCoeff0String )->
          setRestartFlags( RestartFlags::NO_WRITE );
        subRegion->registerWrapper< array1d< real64 > >( viewKeyStruct::apertureAtFailureString )->
          setApplyDefaultValue( -1.0 )->
          setPlotLevel( PlotLevel::LEVEL_0 );

        subRegion->registerWrapper< array1d< real64 > >( FaceElementSubRegion::viewKeyStruct::dSeparationCoeffdAperString )->
          setRestartFlags( RestartFlags::NO_WRITE );
      } );
    } );
  }
}
#endif

void HydrofractureSolver::ImplicitStepSetup( real64 const & time_n,
                                             real64 const & dt,
                                             DomainPartition * const domain,
                                             DofManager & GEOSX_UNUSED_PARAM( dofManager ),
                                             ParallelMatrix & GEOSX_UNUSED_PARAM( matrix ),
                                             ParallelVector & GEOSX_UNUSED_PARAM( rhs ),
                                             ParallelVector & GEOSX_UNUSED_PARAM( solution ) )
{
  /* TJ: calculate the face element effective aperture
   *     and deltaVolume via displacement field of the
   *     fractured face.
   */
  this->UpdateDeformationForCoupling( domain );

  /* TJ: initialize displacement increment to zero;
   *     set beginningOfStepStress stress_n = stress;
   */
  m_solidSolver->ImplicitStepSetup( time_n, dt, domain,
                                    m_solidSolver->getDofManager(),
                                    m_solidSolver->getSystemMatrix(),
                                    m_solidSolver->getSystemRhs(),
                                    m_solidSolver->getSystemSolution() );

  /* TJ: ResetViews()
   *     initialize m_deltaPressure and m_deltaVolume to zero;
   *     set m_elementAperture0 = effective aperture;
   *     set m_densityOld = m_density, m_porosityOld = m_porosity;
   *     UpdateState()
   */
  m_flowSolver->ImplicitStepSetup( time_n, dt, domain,
                                   m_flowSolver->getDofManager(),
                                   m_flowSolver->getSystemMatrix(),
                                   m_flowSolver->getSystemRhs(),
                                   m_flowSolver->getSystemSolution() );

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  domain->getMeshBody( 0 )->getMeshLevel( 0 )->getElemManager()->forElementRegions< FaceElementRegion >( [&]( FaceElementRegion * const faceElemRegion )
  {
    faceElemRegion->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion * const subRegion )
    {
      arrayView1d< real64 > const &
      separationCoeff0 = subRegion->getReference< array1d< real64 > >( viewKeyStruct::separationCoeff0String );
      arrayView1d< real64 const > const &
      separationCoeff = subRegion->getSeparationCoefficient();
      for( localIndex k=0; k<separationCoeff0.size(); ++k )
      {
        separationCoeff0[k] = separationCoeff[k];
      }
    } );
  } );
#endif

}

void HydrofractureSolver::ImplicitStepComplete( real64 const & time_n,
                                                real64 const & dt,
                                                DomainPartition * const domain )
{
  m_flowSolver->ImplicitStepComplete( time_n, dt, domain );
  m_solidSolver->ImplicitStepComplete( time_n, dt, domain );
}

void HydrofractureSolver::PostProcessInput()
{
  string ctOption = this->getReference< string >( viewKeyStruct::couplingTypeOptionStringString );

  if( ctOption == "SIM_FixedStress" )
  {
    this->m_couplingTypeOption = couplingTypeOption::SIM_FixedStress;
  }
  else if( ctOption == "FIM" )
  {
    this->m_couplingTypeOption = couplingTypeOption::FIM;
  }
  else
  {
    GEOSX_ERROR( "invalid coupling type option: " + ctOption );
  }

  m_solidSolver = this->getParent()->GetGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );
  GEOSX_ERROR_IF( m_solidSolver == nullptr, this->getName() << ": invalid solid solver name: " << m_solidSolverName );

  m_flowSolver = this->getParent()->GetGroup< FlowSolverBase >( m_flowSolverName );
  GEOSX_ERROR_IF( m_flowSolver == nullptr, this->getName() << ": invalid flow solver name: " << m_flowSolverName );
}

void HydrofractureSolver::InitializePostInitialConditions_PreSubGroups( Group * const GEOSX_UNUSED_PARAM( problemManager ) )
{}

HydrofractureSolver::~HydrofractureSolver()
{
  // TODO Auto-generated destructor stub
}

void HydrofractureSolver::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  m_flowSolver->ResetStateToBeginningOfStep( domain );
  m_solidSolver->ResetStateToBeginningOfStep( domain );
}

real64 HydrofractureSolver::SolverStep( real64 const & time_n,
                                        real64 const & dt,
                                        int const cycleNumber,
                                        DomainPartition * const domain )
{
  //TJ: This version of solver is a combination of version 1 and version 2.
  //    It uses two-level iterations ONLY when the tip propagates into a new element.
  //    It uses the initial guess method when the tip propagates internally in the tip
  //    element
  real64 dtReturn = dt;

  m_totalTime = time_n + dt;
  //assume the thickness of the KGD problem is 1.0
  real64 const KGDthickness = 1.0;

  {
    //m_meshSize = 1.0;  //hard-coded mesh size, the only place needs to be changed

    SurfaceGenerator * const mySurface = this->getParent()->GetGroup< SurfaceGenerator >( "SurfaceGen" );
    SortedArray< localIndex > const trailingFaces = mySurface->getTrailingFaces();
    Group * elementSubRegions = domain->GetGroup("MeshBodies")
                                      ->GetGroup<MeshBody>("mesh1")
                                      ->GetGroup<MeshLevel>("Level0")
                                      ->GetGroup<ElementRegionManager>("ElementRegions")
                                      ->GetRegion< FaceElementRegion >( "Fracture" )
                                      ->GetGroup("elementSubRegions");

    FaceElementSubRegion * subRegion = elementSubRegions->GetGroup< FaceElementSubRegion >( "default" );
    FaceElementSubRegion::FaceMapType & faceMap = subRegion->faceList();

    for(auto const & trailingFace : trailingFaces)
    {
      bool found = false;
      // loop over all the face element
      for(localIndex i=0; i<faceMap.size(0); i++)
      {
        // loop over all the (TWO) faces in a face element
        for(localIndex j=0; j<faceMap.size(1); j++)
        {
          // if the trailingFace is one of the two faces in a face element,
          // we find it
          if (faceMap[i][j] == trailingFace)
          {
            m_tipElement = i;
            found = true;
            break;
          }
        } // for localIndex j
        if (found)
          break;
      } // for localIndex i
      GEOSX_ASSERT_MSG( found == true,
                      "Trailing face is not found among the fracture face elements" );
    }

    if (subRegion->size() > 0)
    {
      real64 const tipElmtArea = subRegion->getElementArea()[m_tipElement];
      m_meshSize = tipElmtArea/KGDthickness;
      std::cout << "Mesh size = " << m_meshSize << std::endl;
    }
  }

  SolverBase * const surfaceGenerator =  this->getParent()->GetGroup< SolverBase >( "SurfaceGen" );

  if( m_couplingTypeOption == couplingTypeOption::SIM_FixedStress )
  {
    dtReturn = SplitOperatorStep( time_n, dt, cycleNumber, domain->group_cast< DomainPartition * >() );
  }
  else if( m_couplingTypeOption == couplingTypeOption::FIM )
  {

    ImplicitStepSetup( time_n,
                       dt,
                       domain,
                       m_dofManager,
                       m_matrix,
                       m_rhs,
                       m_solution );

    int const maxIter = m_maxNumResolves + 1;
    m_numResolves[1] = m_numResolves[0];
    int solveIter;

    //TJ
    int const rank = MpiWrapper::Comm_rank( MPI_COMM_WORLD );

    //TJ the hard coded elmt length
    real64 meshSize = m_meshSize;  // this value needs to be changed for a mesh-refinement
    //GEOSX_LOG_RANK_0( "Mesh size = " << meshSize );

    //TJ: We use this global flag to ensure that the MPI processes are either ALL in
    //    the initial guess solver (globalSolverIterFlag = 0), or ALL in the tip iteration
    //    solver (globalSolverIterFlag = 1). The reason we use 0/1 instead of true/false is
    //    because we are going to use MPI::allReduce to maintain the consistency of
    //    globalSolverIterFlag across all the processes. The globalSolverIterFlag needs to be
    //    initialized outside the solverIter loop, so that at the beginning of each time step,
    //    we always use the initial guess solver. Also, globalSolverIterFlag should not be rewritten
    //    at the beginning of each solverIter.
    int globalSolverIterFlag = 0;

    for( solveIter=0; solveIter<maxIter; ++solveIter )
    {
      int locallyFractured = 0;
      int globallyFractured = 0;

      //TJ: each copy of hydrofactureSolver (each mpi process) has its own localSolverFlag variable,
      //    even for those process where no surface element (fracture) is defined.
      int localSolverIterFlag = 0;

      /* TJ: create the sparsity patterns for the diagonal blocks
       *     and the off-diagonal blocks
       */
      SetupSystem( domain,
                   m_dofManager,
                   m_matrix,
                   m_rhs,
                   m_solution );

      //TJ: some variables used for print statements
      Group * elementSubRegions = domain->GetGroup("MeshBodies")
					->GetGroup<MeshBody>("mesh1")
					->GetGroup<MeshLevel>("Level0")
					->GetGroup<ElementRegionManager>("ElementRegions")
					->GetRegion< FaceElementRegion >( "Fracture" )
					->GetGroup("elementSubRegions");

      FaceElementSubRegion * subRegion = elementSubRegions->GetGroup< FaceElementSubRegion >( "default" );
      FaceElementSubRegion::NodeMapType & nodeMap = subRegion->nodeList();
//      FaceElementSubRegion::EdgeMapType & edgeMap = subRegion->edgeList();
      FaceElementSubRegion::FaceMapType & faceMap = subRegion->faceList();

      //TJ: original solve without tip iteration
      //    relying on the fact that the initial value of m_tipIterationFlag = false
      //TJ: For MPI, each copy of hydrofractureSolver has its own PRIVATE member of
      //    m_tipIterationFlag. The value of m_tipIterationFlag needs to maintain a consistent
      //    value across all the copies of hydrofractureSolver.
      //TJ: Since we have globalSolverIterFlag now, we can get rid of the private member
      //    m_tipIterationFlag. Also, at the first solver (solverIter = 0), globalSolverIterFlag is
      //    always ZERO

      if (globalSolverIterFlag == 0)
      {
/*
        std::cout << "Rank " << rank << ": I am in initial guess solver." << std::endl;
*/
	real64 localConvergedTipLoc = -1.0;

        SurfaceGenerator * const mySurface = this->getParent()->GetGroup< SurfaceGenerator >( "SurfaceGen" );
	SortedArray< localIndex > & nodesWithAssignedDisp =
	  mySurface->getReference< SortedArray< localIndex > >("nodesWithAssignedDisp");

	//TJ: We have to clear all the essential B.C. assigned at the split nodes
	nodesWithAssignedDisp.clear();

	if( solveIter>0 )
	{
	  /* TJ: set stress = stress_n to reset the solid stress
	   *     to the beginning of the step. On the other hand,
	   *     the displacement field and the displacement increment
	   *     field are NOT reset. That is, the displacement field and
	   *     the displacement increment field used in the first newton
	   *     iteration of the resolve (rewind) are based on the converged
	   *     solution from the previous solve. Of course, the stress is
	   *     reset to the beginning of the step (as mentioned before).
	   *
	   */
	  m_solidSolver->ResetStressToBeginningOfStep( domain );
	}

        MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
        NodeManager & nodeManager = *mesh.getNodeManager();

        array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > & disp = nodeManager.totalDisplacement();
        array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > & dispIncre = nodeManager.incrementalDisplacement();
//        std::cout << dispIncre.size() << std::endl;

	//TJ: calculate some material parameters
	real64 const shearModulus = domain->GetGroup("Constitutive")
			                  ->GetGroup("rock")
			                  ->getReference<real64>("defaultShearModulus");
	real64 const bulkModulus = domain->GetGroup("Constitutive")
		                         ->GetGroup("rock")
				         ->getReference<real64>("defaultBulkModulus");
        real64 const toughness = mySurface->getReference<real64>("rockToughness");
        real64 const viscosity = domain->GetGroup("Constitutive")
                                       ->GetGroup("water")
				       ->getReference<real64>("defaultViscosity");


        // The unit of injectionRate is kg per second
        real64 const injectionRate = domain->getParent()
                                           ->GetGroup<FieldSpecificationManager>("FieldSpecifications")
                                           ->GetGroup<SourceFluxBoundaryCondition>("sourceTerm")
					   ->getReference<real64>("scale");

        // The injectionRate is only for half domain of the KGD problem,
        // to retrieve the full injection rate, we need to multiply it by 2.0
        real64 const q0 = 2.0 * std::abs(injectionRate) /1.0e3;
        real64 const total_time = dt + time_n;

	real64 const nu = ( 1.5 * bulkModulus - shearModulus ) / ( 3.0 * bulkModulus + shearModulus );
	real64 const E = ( 9.0 * bulkModulus * shearModulus )/ ( 3.0 * bulkModulus + shearModulus );
        real64 const Eprime = E/(1.0-nu*nu);
        real64 const PI = 2 * acos(0.0);
        real64 const Kprime = 4.0*sqrt(2.0/PI)*toughness;
        real64 const mup = 12.0 * viscosity;

	//TJ: find the location of the tip boundary
	SortedArray< localIndex > const trailingFaces = mySurface->getTrailingFaces();
	for(auto const & trailingFace : trailingFaces)
	{
	  bool found = false;
	  // loop over all the face element
	  for(localIndex i=0; i<faceMap.size(0); i++)
	  {
	    // loop over all the (TWO) faces in a face element
	    for(localIndex j=0; j<faceMap.size(1); j++)
	    {
	      // if the trailingFace is one of the two faces in a face element,
	      // we find it
	      if (faceMap[i][j] == trailingFace)
	      {
		m_tipElement = i;
		found = true;
		break;
	      }
	    } // for localIndex j
	    if (found)
	      break;
	  } // for localIndex i
	  GEOSX_ASSERT_MSG( found == true,
			"Trailing face is not found among the fracture face elements" );
/*
	  std::cout << "Rank " << rank
	            << ": Before newton solve (globalSolverIterFlag == 0): " << std::endl;
	  std::cout << "Rank " << rank
	            << ": m_tipElement = "     << m_tipElement
		    << ", converged tipLoc = " << m_convergedTipLoc
		    << std::endl;
*/

	  // First time step
	  if (time_n < 1.0e-6 && solveIter == 0)
	  {
	    if (subRegion->size() > 0)
	    {
//	      std::cout << "Rank " << rank << ": Tip element " << m_tipElement << " in the first time step"  << std::endl;
	      real64 refDispFirstStep = 0.0;
	      real64 volume = std::abs(injectionRate) / 1.0e3 * dt;
	      FaceManager & faceManager = *mesh.getFaceManager();
	      r1_array & faceNormal = faceManager.getReference< r1_array >( FaceManager::viewKeyStruct::faceNormalString );
	      SortedArray< localIndex > const tipNodes = mySurface->getTipNodes();

	      if (viscosity < 2.0e-3) // Toughness-dominated case
	      {
		refDispFirstStep = 0.5 * pow(3.0/2.0*
					     pow(Kprime/Eprime,2.0)*
					     volume/1.0,  // assume that the tip thickness is 1
					     1.0/3.0);
	      }
	      else  // Viscosity-dominated case
	      {
		real64 Lm = pow( Eprime*pow(q0,3.0)*pow(total_time,4.0)/mup, 1.0/6.0 );
		real64 gamma_m0 = 0.616;
		real64 velocity = 2.0/3.0 * Lm * gamma_m0 / total_time;
		real64 Betam = pow(2.0, 1.0/3.0) * pow(3.0, 5.0/6.0);
		//real64 lengthFirstStep = pow( 5.0/3.0*volume/Betam*
		//			      pow(Eprime/(mup*velocity),1.0/3.0)
		//			     ,3.0/5.0);
		refDispFirstStep = 0.5 * pow(5.0/3.0*
					     volume/1.0*
					     pow(Betam, 3.0/2.0)*
					     pow(mup*velocity/Eprime, 0.5),
					     2.0/5.0);
	      }
	      //refDispFirstStep = 0.00019;

	      for (auto const & node : nodeMap[m_tipElement])
	      {
		if ( std::find( tipNodes.begin(), tipNodes.end(), node ) == tipNodes.end() )
		{
		  // Insert node to the set for B.C. manipulation
    //		  nodesWithAssignedDisp.insert(node);
		  for (localIndex i=0; i<faceMap[m_tipElement].size(); i++)
		  {
		    auto const & face = faceMap[m_tipElement][i];
		    for (localIndex j=0; j<faceManager.nodeList()[face].size(); j++)
		    {
		      auto const & nodeOnFace = faceManager.nodeList()(face,j);
		      if (node == nodeOnFace)
		      {
			disp(node, 0) = faceNormal(face)[0] > 0 ? -refDispFirstStep : refDispFirstStep;
			dispIncre(node,0) = disp(node,0) - 0.0;
//			std::cout << "Rank " << rank << ": Node " << node << ": " << disp(node,0) << std::endl;
		      }
		    } // for localIndex j
		  } // for localIndex i
		} // if (not found)
	      }  // for auto node
	    } // if (subRegion->size() > 0)
	  } // if (time_n < 1.0e-6 && solveIter == 0)
	}  // 	for(auto const & trailingFace : trailingFaces)

	if (time_n < 1.0e-6 && solveIter == 0)
	{
	  std::map< string, string_array > fieldNames;
	  fieldNames["node"].push_back( keys::IncrementalDisplacement );
	  fieldNames["node"].push_back( keys::TotalDisplacement );
	  fieldNames["elems"].push_back( FlowSolverBase::viewKeyStruct::pressureString );
	  fieldNames["elems"].push_back( "elementAperture" );

	  CommunicationTools::SynchronizeFields( fieldNames,
						 domain->getMeshBody( 0 )->getMeshLevel( 0 ),
						 domain->getNeighbors() );

	  //TJ: update the elmt aperture due to the change of disp field at the newly splitted nodes
	  //    in the first time step
	  this->UpdateDeformationForCoupling( domain );

	  m_flowSolver->ResetViews( domain );
	} //	if (time_n < 1.0e-6 && solveIter == 0)

	//TJ: print out surface aperture before newton solve

	{
/*
	  std::cout << "Rank " << rank << ": Inside initial guess solver, "
		       "before newton solve: "
		    << std::endl;
	  std::cout << "Rank " << rank << ": Face element subRegion size: " << subRegion->size() << std::endl;
	  if (subRegion->size() == 0)
	    std::cout << "Rank " << rank << ": No face element yet" << "\n";
	  else
	  {
	    std::cout << "Rank " << rank << ": Face element info: " << std::endl;
	    for(localIndex i=0; i<nodeMap.size(); i++)
	    {
	      std::cout << "Rank " << rank << ": Face (node) " << i << ": ";
	      for(localIndex j=0; j<nodeMap[i].size(); j++)
		std::cout << nodeMap[i][j] << " ";
	      std::cout << "\n";
	    }
	    for(localIndex i=0; i<edgeMap.size(); i++)
	    {
	      std::cout << "Rank " << rank << ": Face (edge) " << i << ": ";
	      for(localIndex j=0; j<edgeMap[i].size(); j++)
		std::cout << edgeMap[i][j] << " ";
	      std::cout << "\n";
	    }
	    for(localIndex i=0; i<faceMap.size(0); i++)
	    {
	      std::cout << "Rank " << rank << ": Face (face) " << i << ": ";
	      for(localIndex j=0; j<faceMap.size(1); j++)
		std::cout << faceMap[i][j] << " ";
	      std::cout << "\n";
	    }
	    for(localIndex i=0; i<subRegion->getElementCenter().size(); i++)
	    {
	      std::cout << "Rank " << rank << ": Face (center) " << i << ": ";
	      for(auto & item : subRegion->getElementCenter()[i])
		std::cout << item << " ";
	      std::cout << "\n";
	    }
	    for(localIndex i=0; i<subRegion->getElementAperture().size(); i++)
	    {
	      std::cout << "Rank " << rank << ": Face (aperture) " << i << ": "
			<< subRegion->getElementAperture()[i]
			<< std::endl;
	    }
	  }
*/
	} // print statements


	// currently the only method is implicit time integration
	dtReturn = this->NonlinearImplicitStep( time_n,
						dt,
						cycleNumber,
						domain,
						m_dofManager,
						m_matrix,
						m_rhs,
						m_solution );

	//TJ: print out surface aperture after newton solve

	{
/*
	  std::cout << "Rank " << rank << ": Inside initial guess solver, "
		       "after newton solve: "
		    << std::endl;
	  std::cout << "Rank " << rank << ": Face element subRegion size: " << subRegion->size() << std::endl;
	  if (subRegion->size() == 0)
	    std::cout << "Rank " << rank << ": No face element yet" << "\n";
	  else
	  {
	    std::cout << "Rank " << rank << ": Face element info: " << std::endl;
	    for(localIndex i=0; i<nodeMap.size(); i++)
	    {
	      std::cout << "Rank " << rank << ": Face (node) " << i << ": ";
	      for(localIndex j=0; j<nodeMap[i].size(); j++)
		std::cout << nodeMap[i][j] << " ";
	      std::cout << "\n";
	    }
	    for(localIndex i=0; i<edgeMap.size(); i++)
	    {
	      std::cout << "Rank " << rank << ": Face (edge) " << i << ": ";
	      for(localIndex j=0; j<edgeMap[i].size(); j++)
		std::cout << edgeMap[i][j] << " ";
	      std::cout << "\n";
	    }
	    for(localIndex i=0; i<subRegion->getElementCenter().size(); i++)
	    {
	      std::cout << "Rank " << rank << ": Face (center) " << i << ": ";
	      for(auto & item : subRegion->getElementCenter()[i])
		std::cout << item << " ";
	      std::cout << "\n";
	    }
	    for(localIndex i=0; i<subRegion->getElementAperture().size(); i++)
	    {
	      std::cout << "Rank " << rank << ": Face (aperture) " << i << ": "
			<< subRegion->getElementAperture()[i]
			<< std::endl;
	    }
	  }
*/
	} // print statements

/*
        if (time_n < 1.0e-6) //first step
          m_tipElement = 0;
*/
        if (subRegion->size() > 0)
        {
	  real64 const tipElmtArea = subRegion->getElementArea()[m_tipElement];
	  real64 tipElmtSize = sqrt(tipElmtArea);
	  // hard coded
	  tipElmtSize = meshSize;
	  R1Tensor const tipElmtCenter = subRegion->getElementCenter()[m_tipElement];
	  // Tip propagates in the y-direction (hard coded)
	  real64 const tipBCLocation = tipElmtCenter[1] - 0.5 * tipElmtSize;

	  localIndex_array myChildIndex = nodeManager.getReference<localIndex_array>("childIndex");
	  localIndex refNodeIndex=-1;
	  for(localIndex i=0; i<nodeMap[m_tipElement].size(); i++)
	  {
	    localIndex node = nodeMap[m_tipElement][i];
	    if ( myChildIndex[node] >= 0)
	    {
	      refNodeIndex = node;
	      break;
	    }
	  }

	  //TJ: use the displacement gap at the newly split node pair for the tip asymptotic relation
	  real64 refDisp = std::abs( disp(refNodeIndex,0) - disp(myChildIndex[refNodeIndex],0) );
	  GEOSX_ASSERT_MSG( disp(refNodeIndex,0) < 0.0,
			    "Node crosses the symmetric plane." );
/*
	  std::cout << "Rank " << rank << ": refNodeIndex = " << refNodeIndex << ", childIndex = "
					 << myChildIndex[refNodeIndex] << std::endl;
	  std::cout << "Rank " << rank << ": disp " << refNodeIndex               << " = " << disp(refNodeIndex,0)
				       << ", disp " << myChildIndex[refNodeIndex] << " = " << disp(myChildIndex[refNodeIndex],0)
				       << std::endl;
*/
	  real64 tipX = 0.0;
	  if (viscosity < 2.0e-3) // Toughness-dominated case
	  {
	    //TJ: the tip asymptote w = Kprime / Eprime * x^(1/2)
	    tipX = (1.0*refDisp * Eprime / Kprime) * (1.0*refDisp * Eprime / Kprime);
	  }
	  else // Viscosity-dominated case
	  {
/*	    real64 em = pow( mup/(Eprime*total_time), 1.0/3.0 );
	    real64 Lm = pow( Eprime*pow(q0,3.0)*pow(total_time,4.0)/mup, 1.0/6.0 );
            real64 gamma_m0 = 0.616;
	    real64 coeff_viscous = em * Lm * gamma_m0 * sqrt(3.0) * pow(2.0, 2.0/3.0);
            tipX = tipBCLocation / ( pow(coeff_viscous/refDisp, 3.0/2.0) - 1.0 );
*/
	    real64 Lm = pow( Eprime*pow(q0,3.0)*pow(total_time,4.0)/mup, 1.0/6.0 );
            real64 gamma_m0 = 0.616;
	    real64 velocity = 2.0/3.0 * Lm * gamma_m0 / total_time;
	    real64 Betam = pow(2.0, 1.0/3.0) * pow(3.0, 5.0/6.0);
            tipX = sqrt( Eprime/(mup*velocity) * pow( refDisp/Betam ,3.0) );
	  }
	  m_newTipLocation = tipX + tipBCLocation;
	  localConvergedTipLoc = m_newTipLocation;

	  //TJ: when there is no trailingFaces in the MPI process,
	  //    there is no tip element in this process. So we don't
	  //    calculate the tip location in this process.
          if (trailingFaces.empty())
            localConvergedTipLoc = -1.0;
/*
	  std::cout << "Rank " << rank
	            << ": After newton solve (globalSolverIterFlag == 0): " << std::endl;
	  std::cout << "Rank " << rank
		    << ": Tip element = " << m_tipElement
		    << ", "
		    << "converged tip loc = " << localConvergedTipLoc
		    << std::endl;
*/
        }
        else
        {
//	  int const rank = MpiWrapper::Comm_rank( MPI_COMM_WORLD );
//	  std::cout << "Rank " << rank << ": No face element yet" << "\n";
        } // if (subRegion->size() > 0)

        real64 globalConvergedTipLoc = -1.0;
	MpiWrapper::allReduce( &localConvergedTipLoc,
			       &globalConvergedTipLoc,
			       1,
			       MPI_MAX,
			       MPI_COMM_GEOSX );
	m_convergedTipLoc = globalConvergedTipLoc;
/*
	std::cout << "Rank " << rank
		  << ": converged tip loc = " << m_convergedTipLoc
		  << std::endl;
*/
      }  // if (globalSolverIterFlag == 0)
      else  // solve with tip iteration (globalSolverIterFlag == 1)
      {
	int localTipIterConverged = 0;
	int globalTipIterConverged = 0;
	real64 localConvergedTipLoc = -1.0;
/*
        std::cout << "Rank " << rank << ": I am in tip iteration solver." << std::endl;
        std::cout << "Solve with tip iteration, ";
*/
        integer const maxTipIteration = 50;
/*
        std::cout << "maximum " << maxTipIteration << " allowed."
                  << std::endl;
*/
	m_solidSolver->ResetStressToBeginningOfStep( domain );

        SurfaceGenerator * const mySurface = this->getParent()->GetGroup< SurfaceGenerator >( "SurfaceGen" );
        SortedArray< localIndex > const & nodesWithAssignedDisp =
  	    mySurface->getReference< SortedArray< localIndex > >("nodesWithAssignedDisp");
/*
        std::cout << "Rank " << rank << ": size of nodesWithAssignedDisp: "
    	      << nodesWithAssignedDisp.size() << std::endl;
        for(auto & item : nodesWithAssignedDisp)
          std::cout << item << ", ";
        std::cout << std::endl;
*/
        MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
        NodeManager & nodeManager = *mesh.getNodeManager();
        array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > & disp = nodeManager.totalDisplacement();
        array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > & dispIncre = nodeManager.incrementalDisplacement();

	//TJ: calculate some material parameters
	real64 const shearModulus = domain->GetGroup("Constitutive")
			                  ->GetGroup("rock")
			                  ->getReference<real64>("defaultShearModulus");
	real64 const bulkModulus = domain->GetGroup("Constitutive")
		                         ->GetGroup("rock")
				         ->getReference<real64>("defaultBulkModulus");
        real64 const toughness = mySurface->getReference<real64>("rockToughness");
        real64 const viscosity = domain->GetGroup("Constitutive")
                                       ->GetGroup("water")
				       ->getReference<real64>("defaultViscosity");


        // The unit of injectionRate is kg per second
        real64 const injectionRate = domain->getParent()
                                           ->GetGroup<FieldSpecificationManager>("FieldSpecifications")
                                           ->GetGroup<SourceFluxBoundaryCondition>("sourceTerm")
					   ->getReference<real64>("scale");

        // The injectionRate is only for half domain of the KGD problem,
        // to retrieve the full injection rate, we need to multiply it by 2.0
        real64 const q0 = 2.0 * std::abs(injectionRate) /1.0e3;
        real64 const total_time = dt + time_n;

	real64 const nu = ( 1.5 * bulkModulus - shearModulus ) / ( 3.0 * bulkModulus + shearModulus );
	real64 const E = ( 9.0 * bulkModulus * shearModulus )/ ( 3.0 * bulkModulus + shearModulus );
        real64 const Eprime = E/(1.0-nu*nu);
        real64 const PI = 2 * acos(0.0);
        real64 const Kprime = 4.0*sqrt(2.0/PI)*toughness;
        real64 const mup = 12.0 * viscosity;


	//TJ: tolerance for the tip iteration,
	//    convergence is achieved when || m_newTipLocation - m_oldTipLocation || < tipTol
	real64 const tipTol = 1.0e-6;

	//TJ: temporary variables for binary search
	real64 minTipLocation = 0.0;
	real64 maxTipLocation = 0.0;
	localIndex_array myChildIndex = nodeManager.getReference<localIndex_array>("childIndex");

	//TJ: find the location of the channel boundary
	R1Tensor channelElmtCenter;
	real64 channelElmtArea;
	real64 channelElmtSize;
	// Tip propagates in the y-direction (hard coded)
	real64 channelBCLocation;

	localIndex refNodeIndex;
	if (!nodesWithAssignedDisp.empty())
	{
	  channelElmtCenter = subRegion->getElementCenter()[m_channelElement];
	  channelElmtArea = subRegion->getElementArea()[m_channelElement];
	  channelElmtSize = sqrt(channelElmtArea);
	  channelElmtSize = meshSize;  // hard coded
	  // Tip propagates in the y-direction (hard coded)
	  channelBCLocation = channelElmtCenter[1] + 0.5 * channelElmtSize;


	  for(localIndex i=0; i<nodeMap[m_channelElement].size(); i++)
	  {
	    localIndex node = nodeMap[m_channelElement][i];
	    if (   std::find( nodesWithAssignedDisp.begin(), nodesWithAssignedDisp.end(), node )
		== nodesWithAssignedDisp.end()
		&& myChildIndex[node] >= 0)
	    {
	      refNodeIndex = node;
	      break;
	    }
	  }
        } // if (!nodesWithAssignedDisp.empty())
	m_tipLocationHistory.clear();

        //TJ: tip iterations
        for(localIndex tipIterCount = 0; tipIterCount < maxTipIteration; tipIterCount++)
        {

          GEOSX_LOG_RANK_0( "Tip iteration = " << tipIterCount+1);

          //TJ: print out surface aperture before newton solve
	  {
/*
	    std::cout << "Rank " << rank << ": Inside tip iteration, "
			 "before newton solve: "
		      << std::endl;
	    std::cout << "Rank " << rank << ": Face element subRegion size: " << subRegion->size() << std::endl;
	    if (subRegion->size() == 0)
	      std::cout << "Rank " << rank << ": No face element yet" << "\n";
	    else
	    {
	      std::cout << "Rank " << rank << ": Face element info: " << std::endl;
	      for(localIndex i=0; i<nodeMap.size(); i++)
	      {
		std::cout << "Rank " << rank << ": Face (node) " << i << ": ";
		for(localIndex j=0; j<nodeMap[i].size(); j++)
		  std::cout << nodeMap[i][j] << " ";
		std::cout << "\n";
	      }
	      for(localIndex i=0; i<edgeMap.size(); i++)
	      {
		std::cout << "Rank " << rank << ": Face (edge) " << i << ": ";
		for(localIndex j=0; j<edgeMap[i].size(); j++)
		  std::cout << edgeMap[i][j] << " ";
		std::cout << "\n";
	      }
	      for(localIndex i=0; i<faceMap.size(0); i++)
	      {
		std::cout << "Rank " << rank << ": Face (face) " << i << ": ";
		for(localIndex j=0; j<faceMap.size(1); j++)
		  std::cout << faceMap[i][j] << " ";
		std::cout << "\n";
	      }
	      for(localIndex i=0; i<subRegion->getElementCenter().size(); i++)
	      {
		std::cout << "Rank " << rank << ": Face (center) " << i << ": ";
		for(auto & item : subRegion->getElementCenter()[i])
		  std::cout << item << " ";
		std::cout << "\n";
	      }
	      for(localIndex i=0; i<subRegion->getElementAperture().size(); i++)
	      {
		std::cout << "Rank " << rank << ": Face (aperture) " << i << ": "
			  << subRegion->getElementAperture()[i]
			  << std::endl;
	      }
	    }
*/
	  } // print statements

          //TJ: nonlinear solver
	  dtReturn = this->NonlinearImplicitStep( time_n,
						  dt,
						  cycleNumber,
						  domain,
						  m_dofManager,
						  m_matrix,
						  m_rhs,
						  m_solution );
	  //TJ: print out surface aperture after newton solve
	  {
/*
	    std::cout << "Rank " << rank << ": Inside tip iteration, "
			 "after newton solve: "
		      << std::endl;
	    std::cout << "Rank " << rank << ": Face element subRegion size: " << subRegion->size() << std::endl;
	    if (subRegion->size() == 0)
	      std::cout << "Rank " << rank << ": No face element yet" << "\n";
	    else
	    {
	      std::cout << "Rank " << rank << ": Face element info: " << std::endl;
	      for(localIndex i=0; i<nodeMap.size(); i++)
	      {
		std::cout << "Rank " << rank << ": Face (node) " << i << ": ";
		for(localIndex j=0; j<nodeMap[i].size(); j++)
		  std::cout << nodeMap[i][j] << " ";
		std::cout << "\n";
	      }
	      for(localIndex i=0; i<edgeMap.size(); i++)
	      {
		std::cout << "Rank " << rank << ": Face (edge) " << i << ": ";
		for(localIndex j=0; j<edgeMap[i].size(); j++)
		  std::cout << edgeMap[i][j] << " ";
		std::cout << "\n";
	      }
	      for(localIndex i=0; i<subRegion->getElementCenter().size(); i++)
	      {
		std::cout << "Rank " << rank << ": Face (center) " << i << ": ";
		for(auto & item : subRegion->getElementCenter()[i])
		  std::cout << item << " ";
		std::cout << "\n";
	      }
	      for(localIndex i=0; i<subRegion->getElementAperture().size(); i++)
	      {
		std::cout << "Rank " << rank << ": Face (aperture) " << i << ": "
			  << subRegion->getElementAperture()[i]
			  << std::endl;
	      }
	    }
*/
	  } // print statements

	  if (!nodesWithAssignedDisp.empty())
	  {

	    //TJ: obtain new tip location
/*
	    std::cout << "Rank " << rank << ": the face element on the boundary of the channel is "
		      << m_channelElement << std::endl;
	    real64 const refAper = subRegion->getElementAperture()[m_channelElement];
	    std::cout << "Rank " << rank << ": the aperture of the face element at the channel boundary is "
		      << refAper << std::endl;
*/
	    //TJ: find the upper bound of the tip location inside the partially opened element
	    //    the upper bound is obtained at the first tip iteration, since the gap at the
	    //    newly split node pair is almost zero
	    if (tipIterCount == 0)
	    {
	      m_tipLocationHistory.insert(tipIterCount, m_oldTipLocation);
	      minTipLocation = m_oldTipLocation;   // initialize the lower bound
  /*
	      if (m_newTipLocation <= m_oldTipLocation)
	      {
		std::cout << "Accept that new tip location is behind the channel B.C.!" << std::endl;
		break;
	      }
	      //Safe guard for the corner case
	      GEOSX_ASSERT_MSG( m_newTipLocation > m_oldTipLocation,
			    "Corner case: the initial gap (0.0001) is too large!" );
  */
	    }

	    //TJ: use the displacement gap at the node pair other than the newly split one
	    //    on the face element at the channel boundary for the tip asymptotic relation
	    real64 refDisp = std::abs( disp(refNodeIndex,0) - disp(myChildIndex[refNodeIndex],0) );
/*
	    std::cout << "Rank " << rank << ": refNodeIndex = " << refNodeIndex
		      << ", disp "         << disp(refNodeIndex,0)
		      << ", childIndex = " << myChildIndex[refNodeIndex]
		      << ", disp "         << disp(myChildIndex[refNodeIndex],0)
		      << std::endl;
*/

  //        GEOSX_ASSERT_MSG( disp(refNodeIndex,0) < 0.0,
  //			  "Node crosses the symmetric plane." );

	    // TJ: When nodes cross symmetric plane, instead of terminating the simulation,
	    //     we can use the tiploc as a new upper bound
	    if (disp(refNodeIndex,0) >= 0.0)
	    {
/*
	      std::cout << "Rank " << rank << ": Node crosses the symmetric plane. Half the tip location." << std::endl;
*/
	      if (tipIterCount == 0) // initialize the upper bound
		maxTipLocation = m_convergedTipLoc;
	      else
		maxTipLocation = m_oldTipLocation;
	    }
	    else
	    {
	      real64 tipX = 0.0;
	      if (viscosity < 2.0e-3) // Toughness-dominated case
	      {
		//TJ: the tip asymptote w = Kprime / Eprime * x^(1/2)
		tipX = (1.0*refDisp * Eprime / Kprime) * (1.0*refDisp * Eprime / Kprime);
	      }
	      else // Viscosity-dominated case
	      {
    /*	    real64 em = pow( mup/(Eprime*total_time), 1.0/3.0 );
		real64 Lm = pow( Eprime*pow(q0,3.0)*pow(total_time,4.0)/mup, 1.0/6.0 );
		real64 gamma_m0 = 0.616;
		real64 coeff_viscous = em * Lm * gamma_m0 * sqrt(3.0) * pow(2.0, 2.0/3.0);
		tipX = tipBCLocation / ( pow(coeff_viscous/refDisp, 3.0/2.0) - 1.0 );
    */
		real64 Lm = pow( Eprime*pow(q0,3.0)*pow(total_time,4.0)/mup, 1.0/6.0 );
		real64 gamma_m0 = 0.616;
		real64 velocity = 2.0/3.0 * Lm * gamma_m0 / total_time;
		real64 Betam = pow(2.0, 1.0/3.0) * pow(3.0, 5.0/6.0);
		tipX = sqrt( Eprime/(mup*velocity) * pow( refDisp/Betam ,3.0) );
	      }
	      m_newTipLocation = tipX + channelElmtCenter[1] - 0.5 * channelElmtSize;

	      if (tipIterCount == 0)   // initialize the upper bound
	      {
		maxTipLocation =
		  m_newTipLocation > m_convergedTipLoc ? m_newTipLocation : m_convergedTipLoc;
	      }
	      else
	      {
		if (m_newTipLocation <= channelBCLocation)
		{
		  maxTipLocation = m_oldTipLocation;
		}
		else
		{
		  //TJ: update the lower and upper bounds for the binary search
		  if (m_newTipLocation > m_oldTipLocation)
		    minTipLocation = m_oldTipLocation;
		  else
		    maxTipLocation = m_oldTipLocation;
		}
	      }
	    }   // if (disp(refNodeIndex,0) >= 0.0)

	    // safe guard the binary search
  //          GEOSX_ASSERT_MSG( maxTipLocation > minTipLocation,
  //               "maxTipLocation < minTipLocation" );

	    m_newTipLocation = 0.5 * (minTipLocation + maxTipLocation);

	    //TJ: compare old tip location with the new one
	    if ( std::abs(m_newTipLocation - m_oldTipLocation) < tipTol && tipIterCount > 0)
	    {
	      //TJ: change the flag to false to switch to the initial guess based method
	      //TJ: since we have the globalSolverIterFlag now, we don't need m_tipIterationFlag anymore.
              // m_tipIterationFlag = false;
	      m_tipLocationHistory.insert(tipIterCount+1, m_newTipLocation);
	      m_convergedTipLoc = m_newTipLocation;
	      localConvergedTipLoc = m_convergedTipLoc;
	      m_oldTipLocation = m_newTipLocation;

	      for(localIndex i=0; i<m_tipLocationHistory.size(); i++)
		std::cout << "Rank " << rank << ": Tip location (iter " << i << ") = " << m_tipLocationHistory(i) << std::endl;
	      std::cout << "Rank " << rank << ": tipLoc = " << m_convergedTipLoc << ", "
			<< "Tip iteration converges in " << tipIterCount+1 << " steps."
			<< std::endl;

	      localTipIterConverged = 1;
	    }
	  } // if (!nodesWithAssignedDisp.empty())

	  //TJ: each copy of the HF solvers need to break the tip iteration loop
	  //    at the same time
	  MpiWrapper::allReduce( &localTipIterConverged,
				 &globalTipIterConverged,
				 1,
				 MPI_MAX,
				 MPI_COMM_GEOSX );
	  if (globalTipIterConverged == 1)
	  {
	    real64 globalConvergedTipLoc = -1.0;
	    MpiWrapper::allReduce( &localConvergedTipLoc,
				   &globalConvergedTipLoc,
				   1,
				   MPI_MAX,
				   MPI_COMM_GEOSX );
	    m_convergedTipLoc = globalConvergedTipLoc;
	    break;
	  }

	  if (!nodesWithAssignedDisp.empty())
	  {
	    m_tipLocationHistory.insert(tipIterCount+1, m_newTipLocation);
/*
	    for(localIndex i=0; i<m_tipLocationHistory.size(); i++)
	      std::cout << "Rank " << rank << ": Tip location (iter " << i << ") = " << m_tipLocationHistory(i) << std::endl;
*/
	    std::cout << "Rank " << rank
		      << ": tip location (iter " << tipIterCount+1
		      << ") = " << m_newTipLocation
		      << std::endl;
	    std::cout << std::endl;

	    m_oldTipLocation = m_newTipLocation;

	    //TJ: reset essential B.C. at the newly split nodes
	    real64 relativeDist = m_newTipLocation - channelBCLocation;
	    GEOSX_ASSERT_MSG( relativeDist > 0.0,
			    "Tip location falls behind the edge of the newly generated face element!" );

	    real64 refValue = 0.0;
	    if (viscosity < 2.0e-3) // Toughness-dominated case
	    {
	      //TJ: the tip asymptote w = Kprime / Eprime * x^(1/2)
	      // 0.5 is used to get the magnitude of the displacement from the node to the fracture plan
	      // based on symmetry
	      refValue = 0.5 * Kprime/(1.0*Eprime) * sqrt(relativeDist);
	    }
	    else // Viscosity-dominated case
	    {
  /*	    real64 em = pow( mup/(Eprime*total_time), 1.0/3.0 );
	      real64 Lm = pow( Eprime*pow(q0,3.0)*pow(total_time,4.0)/mup, 1.0/6.0 );
	      real64 gamma_m0 = 0.616;
	      real64 coeff_viscous = em * Lm * gamma_m0 * sqrt(3.0) * pow(2.0, 2.0/3.0);
	      // 0.5 is used to get the magnitude of the displacement from the node to the fracture plan
	      // based on symmetry
	      refValue = 0.5 * coeff_viscous * pow(relativeDist/m_convergedTipLoc, 2.0/3.0);
  */
	      real64 Lm = pow( Eprime*pow(q0,3.0)*pow(total_time,4.0)/mup, 1.0/6.0 );
	      real64 gamma_m0 = 0.616;
	      real64 velocity = 2.0/3.0 * Lm * gamma_m0 / total_time;
	      real64 Betam = pow(2.0, 1.0/3.0) * pow(3.0, 5.0/6.0);
	      refValue = 0.5 * Betam
			     * pow( mup*velocity*relativeDist*relativeDist/Eprime, 1.0/3.0);
	    }

	    for(auto node : nodesWithAssignedDisp)
	    {
	      // the sign of the displacement at the newly split nodes are always the same
	      disp(node, 0) = disp(node,0) > 0 ? refValue : -refValue;
	      dispIncre(node, 0) = disp(node, 0) - 0.0;
/*
	      std::cout << "Rank " << rank << ": Node " << node << ": " << disp(node,0) << std::endl;
*/
	    }
	  } // if (!nodesWithAssignedDisp.empty())

	  std::map< string, string_array > fieldNames;
	  fieldNames["node"].push_back( keys::IncrementalDisplacement );
	  fieldNames["node"].push_back( keys::TotalDisplacement );
	  fieldNames["elems"].push_back( FlowSolverBase::viewKeyStruct::pressureString );
	  fieldNames["elems"].push_back( "elementAperture" );

	  CommunicationTools::SynchronizeFields( fieldNames,
						 domain->getMeshBody( 0 )->getMeshLevel( 0 ),
						 domain->getNeighbors() );

          this->UpdateDeformationForCoupling( domain );

          m_flowSolver->ResetViews( domain );

          GEOSX_ASSERT_MSG( tipIterCount < maxTipIteration-1,
  			  "Tip iteration does not converge!" );
        } // for tipIterCount

      } // if (globalSolverIterFlag == 0 or 1)

      /* TJ: update the stress in solid via displacement increment
       *     stress = stress + materialStiffness*disp_increment
       */
      m_solidSolver->updateStress( domain );
//      std::cout << "Rank " << rank << ": after updateStress." << std::endl;

      if( surfaceGenerator!=nullptr )
      {
        if( surfaceGenerator->SolverStep( time_n, dt, cycleNumber, domain ) > 0 )
        {
          locallyFractured = 1;

          // when there is more than one fracture face element, we use tip-based method
          if (subRegion->size() >= 0)
          {
	    /* TJ: This is where we can prescribe the displacement boundary conditions
	     *     for the newly generated nodes due to split. We use the quantities
	     *     m_tipNodes, m_tipEdges, m_tipFaces, and m_trailingFaces defined in the
	     *     SurfaceGenerator class.
	     */
	    SurfaceGenerator * const mySurface = this->getParent()->GetGroup< SurfaceGenerator >( "SurfaceGen" );
	    SortedArray< localIndex > const tipNodes = mySurface->getTipNodes();
/*	    {
	      std::cout << "A new surface is just generated, "
			   "we can manipulate the node displacements "
			   "at the newly generated nodes via split. "
			<< std::endl;
	      std::cout << "Rank " << rank << " m_tipNodes: ";
	      for(auto & item : tipNodes)
		std::cout << item << " ";
	      std::cout << std::endl;
	    }
*/
	    SortedArray< localIndex > const trailingFaces = mySurface->getTrailingFaces();
/*	    {
	      std::cout << "Rank " << rank << " m_trailingFaces: ";
	      for(auto & item : trailingFaces)
		std::cout << item << " ";
	      std::cout << std::endl;
	    }
*/
	    /* TJ: We manipulate the displacement fields at the newly generated
	     *     nodes via element split
	     */
	    MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
	    NodeManager & nodeManager = *mesh.getNodeManager();
	    FaceManager & faceManager = *mesh.getFaceManager();
	    r1_array & faceNormal = faceManager.getReference< r1_array >( FaceManager::viewKeyStruct::faceNormalString );
    /*
	    std::cout << "Face to Node list: " << std::endl;
	    for(int i = 0; i<faceManager.nodeList().size(); i++)
	    {
	      std::cout << "Face "<< i << ": ";
	      for(int j=0; j<faceManager.nodeList()[i].size(); j++)
		std::cout << faceManager.nodeList()(i,j) << " ";
	      std::cout << std::endl;
	    }
	    std::cout << "Face normal: " << std::endl;
	    for(int i=0; i < faceNormal.size(); i++)
	    {
	      std::cout << "Face " << i << ": " ;
	      for(auto & item : faceNormal(i))
		std::cout << item << " ";
	      std::cout << std::endl;
	    }
	    for(localIndex i=0; i<faceMap.size(0); i++)
	    {
	      std::cout << "Face (face) " << i << ": ";
	      for(localIndex j=0; j<faceMap.size(1); j++)
		std::cout << faceMap[i][j] << " ";
	      std::cout << "\n";
	    }
    */
	    array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > & disp = nodeManager.totalDisplacement();
	    array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > & dispIncre = nodeManager.incrementalDisplacement();
    /*
	    std::cout << "Node 10: "
		      << disp(10,0) << ", "
		      << disp(10,1) << ", "
		      << disp(10,2) << std::endl;
	    std::cout << "Node 26: "
		      << disp(26,0) << ", "
		      << disp(26,1) << ", "
		      << disp(26,2) << std::endl;
	    std::cout << "Node 8 (disp): "
		      << disp(8,0) << ", "
		      << disp(8,1) << ", "
		      << disp(8,2) << std::endl;
	    std::cout << "Node 8 (disp_incre): "
		      << dispIncre(8,0) << ", "
		      << dispIncre(8,1) << ", "
		      << dispIncre(8,2) << std::endl;
    */

	    /* TJ: We still need to finish the following task
	     * 0. Assume we only have ONE trailingFace
	     * 1. Assign the displacement field and disp_increment field at newly splitted nodes;
	     * 2. Create a set to include these splitted nodes as essential B.C.;
	     * 3. Pass info from 1. and 2. as essential B.C. values (via getReference);
	     * 4. Update the solidSolver field such as stress (NOT necessary);
	     * 5. Update the fluidSolver field such as aperture UpdateDeformationForCoupling().
	     * 6. Should we worry about other side effects in the flow solver? Shall we call UpdateState()?
	     */
	    //1. manipulate the displacement field on the newly splitted nodes
	    //   it is important to set the displacement increment properly since
	    //   the rhs assembly relys on the the displacement increment.
	    SortedArray< localIndex > & nodesWithAssignedDisp =
	      mySurface->getReference< SortedArray< localIndex > >("nodesWithAssignedDisp");
	    nodesWithAssignedDisp.clear();

	    real64 const shearModulus = domain->GetGroup("Constitutive")
					      ->GetGroup("rock")
					      ->getReference<real64>("defaultShearModulus");
	    real64 const bulkModulus = domain->GetGroup("Constitutive")
					     ->GetGroup("rock")
					     ->getReference<real64>("defaultBulkModulus");
	    real64 const toughness = mySurface->getReference<real64>("rockToughness");

	    real64 const viscosity = domain->GetGroup("Constitutive")
					   ->GetGroup("water")
					   ->getReference<real64>("defaultViscosity");

	    // The unit of injectionRate is kg per second
	    real64 const injectionRate = domain->getParent()
					       ->GetGroup<FieldSpecificationManager>("FieldSpecifications")
					       ->GetGroup<SourceFluxBoundaryCondition>("sourceTerm")
					       ->getReference<real64>("scale");

	    // The injectionRate is only for half domain of the KGD problem,
	    // to retrieve the full injection rate, we need to multiply it by 2.0
	    real64 const q0 = 2.0 * std::abs(injectionRate) /1.0e3;
	    real64 total_time = dt + time_n;

	    real64 const nu = ( 1.5 * bulkModulus - shearModulus ) / ( 3.0 * bulkModulus + shearModulus );
	    real64 const E = ( 9.0 * bulkModulus * shearModulus )/ ( 3.0 * bulkModulus + shearModulus );
	    real64 const Eprime = E/(1.0-nu*nu);
	    real64 const PI = 2 * acos(0.0);
	    real64 const Kprime = 4.0*sqrt(2.0/PI)*toughness;
	    real64 const mup = 12.0 * viscosity;

	    // initial magnitude of the essential B.C. at the newly splited nodes
	    // this initial value can not be zero, since the sign of the displacement
	    // will be used in the later update.
//	    real64 refValue = 0.0025;

/*	    if (subRegion->size() <= 2)
	      refValue = 1.0e-6;
	    else
	      refValue = 1.0e-4;
*/
	    for(auto const & trailingFace : trailingFaces)
	    {
	      bool found = false;
	      // loop over all the face element
	      for(localIndex i=0; i<faceMap.size(0); i++)
	      {
		// loop over all the (TWO) faces in a face element
		for(localIndex j=0; j<faceMap.size(1); j++)
		{
		  // if the trailingFace is one of the two faces in a face element,
		  // we find it
		  if (faceMap[i][j] == trailingFace)
		  {
		    m_tipElement = i;
		    found = true;
		    break;
		  }
		} // for localIndex j
		if (found)
		  break;
	      } // for localIndex i
	      GEOSX_ASSERT_MSG( found == true,
			    "Trailing face is not found among the fracture face elements" );
/*
	      std::cout << "Rank " << rank << ": m_tipElement = " << m_tipElement << std::endl;
*/
	      real64 const tipElmtArea = subRegion->getElementArea()[m_tipElement];
	      //real64 tipElmtSize = sqrt(tipElmtArea);
	      real64 tipElmtSize = tipElmtArea/KGDthickness;
	      //tipElmtSize = meshSize;  // hard coded
	      m_meshSize = tipElmtSize;
	      meshSize = m_meshSize;
	      std::cout << "Mesh size = " << m_meshSize << std::endl;

	      R1Tensor const tipElmtCenter = subRegion->getElementCenter()[m_tipElement];
	      // Tip propagates in the y-direction (hard coded)
	      real64 const tipBCLocation = tipElmtCenter[1] - 0.5 * tipElmtSize;

	      real64 relativeDist = m_convergedTipLoc - tipBCLocation;
//	      GEOSX_ASSERT_MSG( relativeDist < tipElmtSize,
//			    "Tip propagates more than one element, reduce time step." );

	      GEOSX_ASSERT_MSG( relativeDist > 0.0,
			      "Tip location falls behind the edge of the newly generated face element!" );


	      real64 refValue = 0.0;
	      if (viscosity < 2.0e-3) // Toughness-dominated case
	      {
		//TJ: the tip asymptote w = Kprime / Eprime * x^(1/2)
		// 0.5 is used to get the magnitude of the displacement from the node to the fracture plan
		// based on symmetry
		refValue = 0.5 * Kprime/(1.0*Eprime) * sqrt(relativeDist);
	      }
	      else // Viscosity-dominated case
	      {
/*		real64 em = pow( mup/(Eprime*total_time), 1.0/3.0 );
		real64 Lm = pow( Eprime*pow(q0,3.0)*pow(total_time,4.0)/mup, 1.0/6.0 );
		real64 gamma_m0 = 0.616;
		real64 coeff_viscous = em * Lm * gamma_m0 * sqrt(3.0) * pow(2.0, 2.0/3.0);
		// 0.5 is used to get the magnitude of the displacement from the node to the fracture plan
		// based on symmetry
		refValue = 0.5 * coeff_viscous * pow(relativeDist/m_convergedTipLoc, 2.0/3.0);
*/
		real64 Lm = pow( Eprime*pow(q0,3.0)*pow(total_time,4.0)/mup, 1.0/6.0 );
		real64 gamma_m0 = 0.616;
		real64 velocity = 2.0/3.0 * Lm * gamma_m0 / total_time;
		real64 Betam = pow(2.0, 1.0/3.0) * pow(3.0, 5.0/6.0);
		refValue = 0.5 * Betam
		               * pow( mup*velocity*relativeDist*relativeDist/Eprime, 1.0/3.0);
	      }
	      //TJ: We can try to give different magnitude of refValue as initial guess
	      //    to test whether different initial guess will lead to different converged
	      //    solution
//	      refValue = 1.0e-5;

	      //Find which nodes' displacement need to be manipulated
	      for (auto const & node : nodeMap[m_tipElement])
	      {
		if ( std::find( tipNodes.begin(), tipNodes.end(), node ) == tipNodes.end() )
		{
		  // Insert node to the set for B.C. manipulation
		  nodesWithAssignedDisp.insert(node);
		  for (localIndex i=0; i<faceMap[m_tipElement].size(); i++)
		  {
		    auto const & face = faceMap[m_tipElement][i];
		    for (localIndex j=0; j<faceManager.nodeList()[face].size(); j++)
		    {
		      auto const & nodeOnFace = faceManager.nodeList()(face,j);
		      if (node == nodeOnFace)
		      {

			disp(node, 0) = faceNormal(face)[0] > 0 ? -refValue : refValue;
			dispIncre(node,0) = disp(node,0) - 0.0;
/*
			std::cout << "Rank " << rank << ": Node " << node << ": " << disp(node,0) << std::endl;
*/
		      }
		    } // for localIndex j
		  } // for localIndex i
		} // if (not found)
	      }  // for auto node
	    } // for auto trailingFace

	    //TJ: find the element on the boundary of the channel region

	    m_channelElement = -1;
	    for(localIndex i=0; i<nodeMap.size(); i++)
	    {
	      integer nodeCount = 0;
	      for(localIndex j=0; j<nodeMap[i].size(); j++)
	      {
		if ( std::find( nodesWithAssignedDisp.begin(),
			       nodesWithAssignedDisp.end(),
			       nodeMap[i][j] )
		     != nodesWithAssignedDisp.end() )
		  ++nodeCount;
	      }
	      if (nodeCount==4 && i != m_tipElement)
	      {
		m_channelElement = i;
		break;
	      }
	    }
	    GEOSX_ASSERT_MSG( m_channelElement != -1,
		  "Face elmt on the channel boundary is not found!" );
/*
	    std::cout << "Rank " << rank << ": m_channelElement = " << m_channelElement << std::endl;
*/


	    //TJ: assume that the tip location overlap with the tip element boundary

	    real64 const channelElmtArea = subRegion->getElementArea()[m_channelElement];
	    real64 channelElmtSize = sqrt(channelElmtArea);
	    channelElmtSize = meshSize; // hard coded
	    R1Tensor const channelElmtCenter = subRegion->getElementCenter()[m_channelElement];

	    // Tip propagates in the y-direction (hard coded)
	    m_oldTipLocation = channelElmtCenter[1] + 0.5 * channelElmtSize;


    /*
	    // Pair one
	    disp(10,0) = -refValue;
	    dispIncre(10,0) = -refValue - 0.0;
	    disp(26,0) =  refValue;
	    dispIncre(26,0) = refValue - 0.0;
	    // Pair two
	    disp(11,0) = -refValue;
	    dispIncre(11,0) = -refValue - 0.0;
	    disp(27,0) =  refValue;
	    dispIncre(27,0) = refValue - 0.0;

	    //2. create a set to enforce extra essential
	    //    boundary conditions at the newly splitted nodes
	    nodesWithAssignedDisp.insert(10);
	    nodesWithAssignedDisp.insert(11);
	    nodesWithAssignedDisp.insert(26);
	    nodesWithAssignedDisp.insert(27);
    */

	    //TJ: set the tip iteration flag to get ready for finding
	    //    the tip location through fixed-point iteration
            //TJ: since we have globalSolverIterFlag now, we don't need m_tipIterationFlag anymore.
//	    m_tipIterationFlag = true;
	    localSolverIterFlag = 0;
/*
	    std::cout << "Rank " << rank << ": End of disp manipulation" << std::endl;
*/
          } // if subRegion->size() >= 0
        } // if( surfaceGenerator->SolverStep( time_n, dt, cycleNumber, domain ) > 0 )
        MpiWrapper::allReduce( &locallyFractured,
                               &globallyFractured,
                               1,
                               MPI_MAX,
                               MPI_COMM_GEOSX );

        //TJ: we need to maintain a consistent profile of the globalSolverIterFlag
        MpiWrapper::allReduce( &localSolverIterFlag,
                               &globalSolverIterFlag,
                               1,
                               MPI_MAX,
                               MPI_COMM_GEOSX );
      } // if( surfaceGenerator!=nullptr )
      if( globallyFractured == 0 )
      {
        break;
      }
      else
      {
        std::map< string, string_array > fieldNames;
        fieldNames["node"].push_back( keys::IncrementalDisplacement );
        fieldNames["node"].push_back( keys::TotalDisplacement );
        fieldNames["elems"].push_back( FlowSolverBase::viewKeyStruct::pressureString );
        fieldNames["elems"].push_back( "elementAperture" );

        CommunicationTools::SynchronizeFields( fieldNames,
                                               domain->getMeshBody( 0 )->getMeshLevel( 0 ),
                                               domain->getNeighbors() );

        //TJ: check the face elmt aperture due to the nodal displacement change
/*
        {
	  int const rank = MpiWrapper::Comm_rank( MPI_COMM_WORLD );
	  std::cout << "Rank " << rank
	            << ": Face element info (after manipulating disp, before "
		       "UpdateDeformationForCoupling): " << std::endl;
	  for(localIndex i=0; i<nodeMap.size(); i++)
	  {
	    std::cout << "Face (node) " << i << ": ";
	    for(localIndex j=0; j<nodeMap[i].size(); j++)
	      std::cout << nodeMap[i][j] << " ";
	    std::cout << "\n";
	  }
	  for(localIndex i=0; i<subRegion->getElementCenter().size(); i++)
	  {
	    std::cout << "Face (center) " << i << ": ";
	    for(auto & item : subRegion->getElementCenter()[i])
	      std::cout << item << " ";
	    std::cout << "\n";
	  }
	  for(localIndex i=0; i<subRegion->getElementAperture().size(); i++)
	  {
	    std::cout << "Face (aperture) " << i << ": "
		      << subRegion->getElementAperture()[i]
		      << std::endl;
	  }
        }
*/
        //TJ: 5. update the elmt aperture due to the change of disp field at the newly splitted nodes
        this->UpdateDeformationForCoupling( domain );

        //TJ: check the face elmt aperture due to the nodal displacement change
/*
        {
	  int const rank = MpiWrapper::Comm_rank( MPI_COMM_WORLD );
	  std::cout << "Rank "<< rank
	            << " Face element info (after manipulating disp, after "
		       "UpdateDeformationForCoupling): " << std::endl;
	  for(localIndex i=0; i<nodeMap.size(); i++)
	  {
	    std::cout << "Face (node) " << i << ": ";
	    for(localIndex j=0; j<nodeMap[i].size(); j++)
	      std::cout << nodeMap[i][j] << " ";
	    std::cout << "\n";
	  }
	  for(localIndex i=0; i<subRegion->getElementCenter().size(); i++)
	  {
	    std::cout << "Face (center) " << i << ": ";
	    for(auto & item : subRegion->getElementCenter()[i])
	      std::cout << item << " ";
	    std::cout << "\n";
	  }
	  for(localIndex i=0; i<subRegion->getElementAperture().size(); i++)
	  {
	    std::cout << "Face (aperture) " << i << ": "
		      << subRegion->getElementAperture()[i]
		      << std::endl;
	  }
        }
*/

        if( getLogLevel() >= 1 )
        {
          GEOSX_LOG_RANK_0( "++ Fracture propagation. Re-entering Newton Solve." );
        }
        m_flowSolver->ResetViews( domain );

        //TJ: Is this step necessary?
/*
        MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

        forTargetSubRegions( mesh, [&] ( localIndex const targetIndex, ElementSubRegionBase & subRegion1 )
        {
          m_flowSolver->UpdateState( subRegion1, targetIndex );
        } );
*/
      } //  if( globallyFractured == 0 or 1 )
    }  //for( solveIter=0; solveIter<maxIter; ++solveIter )

/*
    if (solveIter == maxIter)
      std::cout << "maxIter = " << maxIter << std::endl;
*/

    GEOSX_ASSERT_MSG( solveIter < maxIter,
	              "Number of max resolve is achieved, but the fracture still propagates"
	              ". Increase maxNumResolves." );

    // final step for completion of timestep. typically secondary variable updates and cleanup.
    /* TJ: in the flowSolver, pressure += delta_pressure, volume += delta_volume
     *                        what is creationMass?
     *     in the solidSolver, update velocity and (acceleration) through disp_increment
     *                         for instance, velocity = disp_incr/dt
     */
    ImplicitStepComplete( time_n, dtReturn, domain );
    m_numResolves[1] = solveIter;

    //TJ: print the converged tip location
    GEOSX_LOG_RANK_0(                 "At time = " << time_n + dt
		       << ", converged tip loc = " << m_convergedTipLoc );

    //TJ: write the pair of time and tip location to a file
    if (rank == 0)
    {
      real64 const viscosityFile = domain->GetGroup("Constitutive")
                                         ->GetGroup("water")
	                                 ->getReference<real64>("defaultViscosity");
      string fileName = "KGD_timeTip_visco_"
	               + std::to_string( viscosityFile )
                       //+ "_meshSize_"
		       //+ std::to_string( meshSize )
                       + ".hist";
      std::ofstream myFile;
      myFile.open(fileName, std::ios::out | std::ios::app);
      myFile << time_n + dt << "\t" << m_convergedTipLoc << "\n";
      myFile.close();
    }

    //TJ: write face element aperture and pressure


  }

  return dtReturn;
}

void HydrofractureSolver::UpdateDeformationForCoupling( DomainPartition * const domain )
{
  //TJ the hard coded elmt length
  real64 const meshSize = m_meshSize;  // this value needs to be changed for a mesh-refinement
  //GEOSX_LOG_RANK_0( "Mesh size = " << meshSize );

  MeshLevel * const meshLevel = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = meshLevel->getElemManager();
  NodeManager const * const nodeManager = meshLevel->getNodeManager();
  FaceManager * const faceManager = meshLevel->getFaceManager();

  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & u = nodeManager->totalDisplacement();
  arrayView1d< R1Tensor const > const & faceNormal = faceManager->faceNormal();
  // arrayView1d<real64 const> const & faceArea = faceManager->faceArea();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager->nodeList().toViewConst();

  ConstitutiveManager const * const
  constitutiveManager = domain->GetGroup< ConstitutiveManager >( keys::ConstitutiveManager );

  ContactRelationBase const * const
  contactRelation = constitutiveManager->GetGroup< ContactRelationBase >( m_contactRelationName );

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    arrayView1d< real64 > const & aperture = subRegion.getElementAperture();
    arrayView1d< real64 > const & effectiveAperture = subRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::effectiveApertureString );
    arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
    arrayView1d< real64 > const & deltaVolume = subRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::deltaVolumeString );
    arrayView1d< real64 const > const & area = subRegion.getElementArea();
    arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
    arrayView1d< real64 const > const &
    apertureF = subRegion.getReference< array1d< real64 > >( viewKeyStruct::apertureAtFailureString );

    arrayView1d< real64 > const &
    separationCoeff = subRegion.getSeparationCoefficient();

    arrayView1d< real64 > const &
    dSeparationCoeff_dAper = subRegion.getReference< array1d< real64 > >( FaceElementSubRegion::viewKeyStruct::dSeparationCoeffdAperString );
    arrayView1d< real64 const > const &
    separationCoeff0 = subRegion.getReference< array1d< real64 > >( viewKeyStruct::separationCoeff0String );
#endif

    for( localIndex kfe=0; kfe<subRegion.size(); ++kfe )
    {
      localIndex const kf0 = elemsToFaces[kfe][0];
      localIndex const kf1 = elemsToFaces[kfe][1];
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );
      R1Tensor temp;
      for( localIndex a=0; a<numNodesPerFace; ++a )
      {
        temp += u[faceToNodeMap( kf0, a )];
        temp -= u[faceToNodeMap( kf1, a )];
      }

      // TODO this needs a proper contact based strategy for aperture
      //TJ: How do we know aperture is positive?
      aperture[kfe] = -Dot( temp, faceNormal[kf0] ) / numNodesPerFace;
      //TJ: Do we use effectiveAperture as the actual aperture?
      effectiveAperture[kfe] = contactRelation->effectiveAperture( aperture[kfe] );
/*
      if (effectiveAperture[kfe] < 0.0)
	effectiveAperture[kfe] = contactRelation->effectiveAperture( aperture[kfe] );
      else
	effectiveAperture[kfe] = aperture[kfe];
*/


#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
      real64 const s = aperture[kfe] / apertureF[kfe];
      if( separationCoeff0[kfe]<1.0 && s>separationCoeff0[kfe] )
      {
        if( s >= 1.0 )
        {
          separationCoeff[kfe] = 1.0;
          dSeparationCoeff_dAper[kfe] = 0.0;
        }
        else
        {
          separationCoeff[kfe] = s;
          dSeparationCoeff_dAper[kfe] = 1.0/apertureF[kfe];
        }
      }
#endif
      deltaVolume[kfe] = effectiveAperture[kfe] * area[kfe] - volume[kfe];
    }

    //TJ: modify the aperture and deltaVolume at the tip element
    //    according to the nonlinear relationship between aperture
    //    and nodal displacement at the tip
    SurfaceGenerator * const mySurface = this->getParent()->GetGroup< SurfaceGenerator >( "SurfaceGen" );
    SortedArray< localIndex > const trailingFaces = mySurface->getTrailingFaces();
    FaceElementSubRegion::FaceMapType & faceMap = subRegion.faceList();

    for(auto const & trailingFace : trailingFaces)
    {
      bool found = false;
      // loop over all the face element
      for(localIndex i=0; i<faceMap.size(0); i++)
      {
	// loop over all the (TWO) faces in a face element
	for(localIndex j=0; j<faceMap.size(1); j++)
	{
	  // if the trailingFace is one of the two faces in a face element,
	  // we find it
	  if (faceMap[i][j] == trailingFace)
	  {
	    m_tipElement = i;
	    found = true;
	    break;
	  }
	} // for localIndex j
	if (found)
	  break;
      } // for localIndex i
      GEOSX_ASSERT_MSG( found == true,
		    "Trailing face is not found among the fracture face elements" );
//      std::cout << "UpdateDeformationForCoupling :: m_tipElement = " << m_tipElement << std::endl;

      //MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
//      NodeManager & const myNodeManager = meshLevel->getNodeManager();

      MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
      NodeManager & myNodeManager = *mesh.getNodeManager();

      array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > & disp = myNodeManager.totalDisplacement();
      array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > & dispIncre = myNodeManager.incrementalDisplacement();

      FaceElementSubRegion::NodeMapType & nodeMap = subRegion.nodeList();
      localIndex_array myChildIndex = myNodeManager.getReference<localIndex_array>("childIndex");
      localIndex refNodeIndex=-1;
      for(localIndex i=0; i<nodeMap[m_tipElement].size(); i++)
      {
	localIndex node = nodeMap[m_tipElement][i];
	if ( myChildIndex[node] >= 0)
	{
	  refNodeIndex = node;
	  break;
	}
      }



      //TJ: should check refDisp, instead of individual displacement value
//      GEOSX_ASSERT_MSG( u(refNodeIndex,0) < 1.0e-12,
//      		"Node crosses the symmetric plane." );

      //TJ: if nodes cross they symmetric plane, do not terminate.
      //    instead, change the displacement to zero
      if (u(refNodeIndex,0) > 0.0)
      {
	for(localIndex i=0; i<nodeMap[m_tipElement].size(); i++)
	{
	  localIndex node = nodeMap[m_tipElement][i];
	  if ( myChildIndex[node] >= 0)
	  {
	    std::cout << "Nodes cross symmetric plane, zero the displacement of node "
		      << node << " and " << myChildIndex[node] << std::endl;
	    disp(node,0)=-1.0e-4;
	    dispIncre(node,0)=-1.0e-4;
	    disp(myChildIndex[node],0)=1.0e-4;
	    dispIncre(myChildIndex[node],0)=1.0e-4;
	  }
	}
      }
      //TJ: use the displacement gap at the newly split node pair for the tip deltaVolume
      real64 refDisp = std::abs( u(refNodeIndex,0) - u(myChildIndex[refNodeIndex],0) );

      std::cout  << " disp " << refNodeIndex               << " = " << u(refNodeIndex,0)
	         << ", disp " << myChildIndex[refNodeIndex] << " = " << u(myChildIndex[refNodeIndex],0)
	                      << std::endl;
      std::cout << "refDisp = " << refDisp << std::endl;

      real64 const shearModulus = domain->GetGroup("Constitutive")
			                  ->GetGroup("rock")
			                  ->getReference<real64>("defaultShearModulus");
      real64 const bulkModulus = domain->GetGroup("Constitutive")
		                         ->GetGroup("rock")
				         ->getReference<real64>("defaultBulkModulus");
      real64 const toughness = mySurface->getReference<real64>("rockToughness");
      real64 const viscosity = domain->GetGroup("Constitutive")
                                     ->GetGroup("water")
				       ->getReference<real64>("defaultViscosity");


      // The unit of injectionRate is kg per second
      real64 const injectionRate = domain->getParent()
                                         ->GetGroup<FieldSpecificationManager>("FieldSpecifications")
                                         ->GetGroup<SourceFluxBoundaryCondition>("sourceTerm")
					   ->getReference<real64>("scale");

      // The injectionRate is only for half domain of the KGD problem,
      // to retrieve the full injection rate, we need to multiply it by 2.0
      real64 const q0 = 2.0 * std::abs(injectionRate) /1.0e3;
      real64 const total_time = m_totalTime;

      real64 const nu = ( 1.5 * bulkModulus - shearModulus ) / ( 3.0 * bulkModulus + shearModulus );
      real64 const E = ( 9.0 * bulkModulus * shearModulus )/ ( 3.0 * bulkModulus + shearModulus );
      real64 const Eprime = E/(1.0-nu*nu);
      real64 const PI = 2 * acos(0.0);
      real64 const Kprime = 4.0*sqrt(2.0/PI)*toughness;
      real64 const mup = 12.0 * viscosity;

      real64 vTip;
      real64 aperTip;

      if (viscosity < 2.0e-3) // Toughness-dominated case
      {
	//TJ: the tip asymptote w = Kprime / Eprime * x^(1/2)
	vTip = 2.0/3.0 * pow(Eprime/Kprime, 2.0) * pow(refDisp, 3.0);
	// average aperture
	aperTip = vTip/meshSize;
	aperTip = 2.0/3.0 * refDisp;
      }
      else // Viscosity-dominated case
      {
	real64 Lm = pow( Eprime*pow(q0,3.0)*pow(total_time,4.0)/mup, 1.0/6.0 );
	real64 gamma_m0 = 0.616;
	real64 velocity = 2.0/3.0 * Lm * gamma_m0 / total_time;
	real64 Betam = pow(2.0, 1.0/3.0) * pow(3.0, 5.0/6.0);
	vTip = 3.0/5.0 * pow(Betam, -3.0/2.0) * pow( Eprime/mup/velocity, 1.0/2.0)
	               * pow(refDisp, 5.0/2.0);
	// average aperture
	aperTip = vTip/meshSize;
	aperTip = 3.0/5.0 * refDisp;
      }


      aperture[m_tipElement] = aperTip;
      effectiveAperture[m_tipElement] = contactRelation->effectiveAperture( aperture[m_tipElement] );

      deltaVolume[m_tipElement] = vTip/meshSize * area[m_tipElement]
				- volume[m_tipElement];
//      std::cout << "Tip element " << m_tipElement
//	        << " : delta volume = " << deltaVolume[m_tipElement] << std::endl;

//      std::cout << "End of tip volume manipulation." << std::endl;
    }


  } );
}

real64 HydrofractureSolver::SplitOperatorStep( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                               real64 const & dt,
                                               integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                               DomainPartition * const GEOSX_UNUSED_PARAM( domain ) )
{
  real64 dtReturn = dt;
//  real64 dtReturnTemporary = dtReturn;
//
//  m_flowSolver->ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );
//  m_solidSolver->ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );
//  this->ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );
//
//
//
//  fluidSolver.ImplicitStepSetup( time_n, dt, domain,
//                                 fluidSolver.getDofManager(),
//                                 fluidSolver.getSystemMatrix(),
//                                 fluidSolver.getSystemRhs(),
//                                 fluidSolver.getSystemSolution() );
//
//  solidSolver.ImplicitStepSetup( time_n, dt, domain,
//                                 solidSolver.getDofManager(),
//                                 solidSolver.getSystemMatrix(),
//                                 solidSolver.getSystemRhs(),
//                                 solidSolver.getSystemSolution() );
//
//  this->UpdateDeformationForCoupling(domain);
//
//  int iter = 0;
//  while (iter < solverParams->maxIterNewton() )
//  {
//    if (iter == 0)
//    {
//      // reset the states of all slave solvers if any of them has been reset
//      m_flowSolver->ResetStateToBeginningOfStep( domain );
//      m_solidSolver->ResetStateToBeginningOfStep( domain );
//      ResetStateToBeginningOfStep( domain );
//    }
//    LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1  << ", FlowSolver: " );
//
//    // call assemble to fill the matrix and the rhs
//    m_flowSolver->AssembleSystem( domain, getLinearSystemRepository(), time_n+dt, dt );
//
//    // apply boundary conditions to system
//    m_flowSolver->ApplyBoundaryConditions( domain, getLinearSystemRepository(), time_n, dt );
//
//    // call the default linear solver on the system
//    m_flowSolver->SolveSystem( getLinearSystemRepository(),
//                 getSystemSolverParameters() );
//
//    // apply the system solution to the fields/variables
//    m_flowSolver->ApplySystemSolution( getLinearSystemRepository(), 1.0, domain );
//
//    if (dtReturnTemporary < dtReturn)
//    {
//      iter = 0;
//      dtReturn = dtReturnTemporary;
//      continue;
//    }
//
////    if (m_fluidSolver->getSystemSolverParameters()->numNewtonIterations() == 0 && iter > 0 && getLogLevel() >= 1)
////    {
////      GEOSX_LOG_RANK_0( "***** The iterative coupling has converged in " << iter  << " iterations! *****\n" );
////      break;
////    }
//
//    if (getLogLevel() >= 1)
//    {
//      GEOSX_LOG_RANK_0( "\tIteration: " << iter+1  << ", MechanicsSolver: " );
//    }
//
//    // call assemble to fill the matrix and the rhs
//    m_solidSolver->AssembleSystem( domain, getLinearSystemRepository(), time_n+dt, dt );
//
//
//    ApplyFractureFluidCoupling( domain, *getLinearSystemRepository() );
//
//    // apply boundary conditions to system
//    m_solidSolver->ApplyBoundaryConditions( domain, getLinearSystemRepository(), time_n, dt );
//
//    // call the default linear solver on the system
//    m_solidSolver->SolveSystem( getLinearSystemRepository(),
//                 getSystemSolverParameters() );
//
//    // apply the system solution to the fields/variables
//    m_solidSolver->ApplySystemSolution( getLinearSystemRepository(), 1.0, domain );
//
//    if( m_flowSolver->CalculateResidualNorm( getLinearSystemRepository(), domain ) < solverParams->newtonTol() &&
//        m_solidSolver->CalculateResidualNorm( getLinearSystemRepository(), domain ) < solverParams->newtonTol() )
//    {
//      GEOSX_LOG_RANK_0( "***** The iterative coupling has converged in " << iter  << " iterations! *****\n" );
//      break;
//    }
//
//    if (dtReturnTemporary < dtReturn)
//    {
//      iter = 0;
//      dtReturn = dtReturnTemporary;
//      continue;
//    }
////    if (m_solidSolver->getSystemSolverParameters()->numNewtonIterations() > 0)
//    {
//      this->UpdateDeformationForCoupling(domain);
////      m_fluidSolver->UpdateState(domain);
//    }
//    ++iter;
//  }
//
//  this->ImplicitStepComplete( time_n, dt, domain );

  return dtReturn;
}

real64 HydrofractureSolver::ExplicitStep( real64 const & time_n,
                                          real64 const & dt,
                                          const int cycleNumber,
                                          DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->ExplicitStep( time_n, dt, cycleNumber, domain );
  m_flowSolver->SolverStep( time_n, dt, cycleNumber, domain );

  return dt;
}


void HydrofractureSolver::SetupDofs( DomainPartition const * const domain,
                                     DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->SetupDofs( domain, dofManager );
  m_flowSolver->SetupDofs( domain, dofManager );

  // restrict coupling to fracture regions only (as done originally in SetupSystem)
  ElementRegionManager const * const elemManager = domain->getMeshBody( 0 )->getMeshLevel( 0 )->getElemManager();
  string_array fractureRegions;
  elemManager->forElementRegions< FaceElementRegion >( [&]( FaceElementRegion const & elementRegion )
  {
    fractureRegions.push_back( elementRegion.getName() );
  } );

  dofManager.addCoupling( keys::TotalDisplacement,
                          FlowSolverBase::viewKeyStruct::pressureString,
                          DofManager::Connector::Elem,
                          fractureRegions );
}

void HydrofractureSolver::SetupSystem( DomainPartition * const domain,
                                       DofManager & dofManager,
                                       ParallelMatrix & matrix,
                                       ParallelVector & rhs,
                                       ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;
  m_flowSolver->ResetViews( domain );
  //TJ: setup solidSolver sparsity pattern
  m_solidSolver->SetupSystem( domain,
                              m_solidSolver->getDofManager(),
                              m_solidSolver->getSystemMatrix(),
                              m_solidSolver->getSystemRhs(),
                              m_solidSolver->getSystemSolution() );
  //TJ: setup sparsity patterns for fluidSolver and m_derivativeFluxResidual_dAperture
  m_flowSolver->SetupSystem( domain,
                             m_flowSolver->getDofManager(),
                             m_flowSolver->getSystemMatrix(),
                             m_flowSolver->getSystemRhs(),
                             m_flowSolver->getSystemSolution() );

  // setup coupled DofManager
  m_dofManager.setMesh( domain, 0, 0 );
  SetupDofs( domain, dofManager );

  // This is needed for NonlinearImplicitStep to function, even though
  // we're currently not assembling the monolithic matrix/rhs
  localIndex const numLocalDof = dofManager.numLocalDofs();
  matrix.createWithLocalSize( numLocalDof, numLocalDof, 0, MPI_COMM_GEOSX );
  rhs.createWithLocalSize( numLocalDof, MPI_COMM_GEOSX );
  solution.createWithLocalSize( numLocalDof, MPI_COMM_GEOSX );

  // Emulate assembly
  matrix.open();
  matrix.close();

  // By not calling dofManager.reorderByRank(), we keep separate dof numbering for each field,
  // which allows constructing separate sparsity patterns for off-diagonal blocks of the matrix.
  // Once the solver moves to monolithic matrix, we can remove this method and just use SolverBase::SetupSystem.

  /* TJ: m_matrix01 and m_matrix10 are off-diagonal,
   *     and the rest of the code is to set up the
   *     sparsity patterns for m_matrix01 and m_matrix10
   */
  m_matrix01.createWithLocalSize( m_solidSolver->getSystemMatrix().numLocalRows(),
                                  m_flowSolver->getSystemMatrix().numLocalCols(),
                                  9,
                                  MPI_COMM_GEOSX );
  /* TJ: bug, m_matrix10 should be m_flowSolver  ->getSystemMatrix().numLocalRows()
   *                            by m_solidSolver ->getSystemMatrix().numLocalCols()
   */
  m_matrix10.createWithLocalSize( m_flowSolver->getSystemMatrix().numLocalCols(),
                                  m_solidSolver->getSystemMatrix().numLocalRows(),
                                  24,
                                  MPI_COMM_GEOSX );

#if 0
  dofManager.setSparsityPattern( m_matrix01, keys::TotalDisplacement, FlowSolverBase::viewKeyStruct::pressureString );
  dofManager.setSparsityPattern( m_matrix10, FlowSolverBase::viewKeyStruct::pressureString, keys::TotalDisplacement );
#else
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  NodeManager * const nodeManager = mesh->getNodeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

  string const presDofKey = m_flowSolver->getDofManager().getKey( FlowSolverBase::viewKeyStruct::pressureString );
  string const dispDofKey = m_solidSolver->getDofManager().getKey( keys::TotalDisplacement );

  arrayView1d< globalIndex > const &
  dispDofNumber =  nodeManager->getReference< globalIndex_array >( dispDofKey );

  m_matrix01.open();
  m_matrix10.open();

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const & elementSubRegion )
  {
    localIndex const numElems = elementSubRegion.size();
    array1d< array1d< localIndex > > const & elemsToNodes = elementSubRegion.nodeList();
    arrayView1d< globalIndex const > const &
    faceElementDofNumber = elementSubRegion.getReference< array1d< globalIndex > >( presDofKey );

    for( localIndex k=0; k<numElems; ++k )
    {
      globalIndex const activeFlowDOF = faceElementDofNumber[k];
      localIndex const numNodesPerElement = elemsToNodes[k].size();
      array1d< globalIndex > activeDisplacementDOF( 3 * numNodesPerElement );
      array1d< real64 > values( 3*numNodesPerElement );
      values = 1;

      for( localIndex a=0; a<numNodesPerElement; ++a )
      {
        for( int d=0; d<3; ++d )
        {
          activeDisplacementDOF[a * 3 + d] = dispDofNumber[elemsToNodes[k][a]] + d;
        }
      }

      m_matrix01.insert( activeDisplacementDOF.data(),
                         &activeFlowDOF,
                         values.data(),
                         activeDisplacementDOF.size(),
                         1 );

      m_matrix10.insert( &activeFlowDOF,
                         activeDisplacementDOF.data(),
                         values.data(),
                         1,
                         activeDisplacementDOF.size() );
    }
  } );

  NumericalMethodsManager const * numericalMethodManager =
    domain->getParent()->GetGroup< NumericalMethodsManager >( keys::numericalMethodsManager );

  FiniteVolumeManager const * fvManager =
    numericalMethodManager->GetGroup< FiniteVolumeManager >( keys::finiteVolumeManager );

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation( m_flowSolver->getDiscretization() );


  fluxApprox->forStencils< FaceElementStencil >( [&]( FaceElementStencil const & stencil )
  {
    for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );
      typename FaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
      typename FaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
      typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

      FaceElementSubRegion const * const
      elementSubRegion = elemManager->GetRegion( seri[iconn][0] )->GetSubRegion< FaceElementSubRegion >( sesri[iconn][0] );

      array1d< array1d< localIndex > > const & elemsToNodes = elementSubRegion->nodeList();

      arrayView1d< globalIndex const > const &
      faceElementDofNumber = elementSubRegion->getReference< array1d< globalIndex > >( presDofKey );
      for( localIndex k0=0; k0<numFluxElems; ++k0 )
      {
        globalIndex const activeFlowDOF = faceElementDofNumber[sei[iconn][k0]];

        for( localIndex k1=0; k1<numFluxElems; ++k1 )
        {
          localIndex const numNodesPerElement = elemsToNodes[sei[iconn][k1]].size();
          array1d< globalIndex > activeDisplacementDOF( 3 * numNodesPerElement );
          array1d< real64 > values( 3*numNodesPerElement );
          values = 1;

          for( localIndex a=0; a<numNodesPerElement; ++a )
          {
            for( int d=0; d<3; ++d )
            {
              activeDisplacementDOF[a * 3 + d] = dispDofNumber[elemsToNodes[sei[iconn][k1]][a]] + d;
            }
          }

          m_matrix10.insert( &activeFlowDOF,
                             activeDisplacementDOF.data(),
                             values.data(),
                             1,
                             activeDisplacementDOF.size() );
        }
      }
    }//);
  } );

  m_matrix01.close();
  m_matrix10.close();
#endif
}

void HydrofractureSolver::AssembleSystem( real64 const time,
                                          real64 const dt,
                                          DomainPartition * const domain,
                                          DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                          ParallelMatrix & GEOSX_UNUSED_PARAM( matrix ),
                                          ParallelVector & GEOSX_UNUSED_PARAM( rhs ) )
{
  GEOSX_MARK_FUNCTION;

  m_solidSolver->getSystemMatrix().zero();
  m_solidSolver->getSystemRhs().zero();


  /* TJ: Inside the solid solver, the system stiffness matrix dRdU is assembled
   *     in the normal way. When the RHS are assembled, the RHS uses the stress
   *     at the beginning of the step stress_n, plus the dRdU * disp_increment,
   *     that is, RHS = stress_n + dRdU * incrementalDisplacement. Therefore,
   *     incrementalDisplacement is important.
   */
  m_solidSolver->AssembleSystem( time,
                                 dt,
                                 domain,
                                 m_solidSolver->getDofManager(),
                                 m_solidSolver->getSystemMatrix(),
                                 m_solidSolver->getSystemRhs() );

  m_flowSolver->getSystemMatrix().zero();
  m_flowSolver->getSystemRhs().zero();

  m_flowSolver->AssembleSystem( time,
                                dt,
                                domain,
                                m_flowSolver->getDofManager(),
                                m_flowSolver->getSystemMatrix(),
                                m_flowSolver->getSystemRhs() );


  m_matrix01.zero();
  AssembleForceResidualDerivativeWrtPressure( domain, &m_matrix01, &(m_solidSolver->getSystemRhs()) );

  m_matrix10.zero();
  AssembleFluidMassResidualDerivativeWrtDisplacement( domain, &m_matrix10, &(m_flowSolver->getSystemRhs()) );
}

void HydrofractureSolver::ApplyBoundaryConditions( real64 const time,
                                                   real64 const dt,
                                                   DomainPartition * const domain,
                                                   DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                                   ParallelMatrix & GEOSX_UNUSED_PARAM( matrix ),
                                                   ParallelVector & GEOSX_UNUSED_PARAM( rhs ) )
{

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );

  FieldSpecificationManager const & fsManager = FieldSpecificationManager::get();
  string const dispDofKey = m_solidSolver->getDofManager().getKey( keys::TotalDisplacement );
  NodeManager const * const nodeManager = mesh->getNodeManager();

  //TJ: test whether m_nodesWithAssignedDisp is passed here correctly
  {
/*
    int const rank = MpiWrapper::Comm_rank( MPI_COMM_WORLD );
    SurfaceGenerator * const mySurface = this->getParent()->GetGroup< SurfaceGenerator >( "SurfaceGen" );
    SortedArray< localIndex > const & nodesWithAssignedDisp =
    mySurface->getReference< SortedArray< localIndex > >("nodesWithAssignedDisp");
    std::cout << "Rank " << rank << ": size of nodesWithAssignedDisp: "
	      << nodesWithAssignedDisp.size() << std::endl;
    for(auto & item : nodesWithAssignedDisp)
      std::cout << item << ", ";
    std::cout << std::endl;
*/
  }

  //TJ: test whether the displacement field at newly splitted nodes are passed correctly
  {
/*
    int const rank = MpiWrapper::Comm_rank( MPI_COMM_WORLD );
    SurfaceGenerator * const mySurface = this->getParent()->GetGroup< SurfaceGenerator >( "SurfaceGen" );
    SortedArray< localIndex > const & nodesWithAssignedDisp =
    mySurface->getReference< SortedArray< localIndex > >("nodesWithAssignedDisp");
    if (!nodesWithAssignedDisp.empty())
    {
      arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const &
      disp = nodeManager->totalDisplacement();
      arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const &
      dispIncre = nodeManager->incrementalDisplacement();

      for(auto node : nodesWithAssignedDisp)
      {
	std::cout << "Rank " << rank
	          << ": Node " << node << " (disp): "
	          << disp(node, 0) << ", "
	          << disp(node, 1) << ", "
	          << disp(node, 2) << ", "
                  << std::endl;
	std::cout << "Rank " << rank
	          << ": Node " << node << " (dispIncre): "
	          << dispIncre(node, 0) << ", "
	          << dispIncre(node, 1) << ", "
	          << dispIncre(node, 2) << ", "
                  << std::endl;
      }
    }
*/
  }

  GEOSX_MARK_FUNCTION;
  m_solidSolver->ApplyBoundaryConditions( time,
                                          dt,
                                          domain,
                                          m_solidSolver->getDofManager(),
                                          m_solidSolver->getSystemMatrix(),
                                          m_solidSolver->getSystemRhs() );

  //TJ: Apply the assigned displacement field at newly generated nodes as
  //    essential B.C.
  SurfaceGenerator * const mySurface = this->getParent()->GetGroup< SurfaceGenerator >( "SurfaceGen" );
  SortedArray< localIndex > const & nodesWithAssignedDisp =
          mySurface->getReference< SortedArray< localIndex > >("nodesWithAssignedDisp");

  //TJ: matrix open and close have to be outside the   if (!nodesWithAssignedDisp.empty())
  m_solidSolver->getSystemMatrix().open();
  m_solidSolver->getSystemRhs().open();
  if (!nodesWithAssignedDisp.empty())
  {
/*    integer newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;

    {
      string filename = "before_matrix00_" + std::to_string( time ) + "_" + std::to_string( newtonIter ) + ".mtx";
      m_solidSolver->getSystemMatrix().write( filename, LAIOutputFormat::MATRIX_MARKET );
      GEOSX_LOG_RANK_0( "matrix00: written to " << filename );
    }
    {
      string filename = "before_residual0_" + std::to_string( time ) + "_" + std::to_string( newtonIter ) + ".mtx";
      m_solidSolver->getSystemRhs().write( filename, LAIOutputFormat::MATRIX_MARKET );
      GEOSX_LOG_RANK_0( "residual0: written to " << filename );
    }
*/

    // Hard code for displacement in the x-direction
    integer const component = 0;
    arrayView1d< globalIndex const > const & dofMap
                                        = nodeManager->getReference< array1d< globalIndex > >( dispDofKey );
    real64_array rhsContribution( nodesWithAssignedDisp.size() );
    globalIndex_array dof( nodesWithAssignedDisp.size() );

    arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & disp
                                        = nodeManager->totalDisplacement();
    //std::cout << "We are before zero out matrix00 entries." << std::endl;

    integer counter = 0;
    for(auto node : nodesWithAssignedDisp)
    {
      //std::cout << "We are inside Node: " << node << std::endl;

      dof( counter ) = dofMap[node]+component;
      real64 bcValue = disp(node, component);
      real64 fieldValue = bcValue;

      FieldSpecificationEqual::template SpecifyFieldValue<LAInterface>( dof( counter ),
									m_solidSolver->getSystemMatrix(),
				                                        rhsContribution( counter ),
				                                        bcValue,
				                                        fieldValue );

      ++counter;
    }

    //std::cout << "We passed zero out matrix00 entries." << std::endl;
    //std::cout << "We are before prescribe RHS due to zero out matrix00 entries." << std::endl;

    FieldSpecificationEqual::template PrescribeRhsValues< LAInterface >( m_solidSolver->getSystemRhs(),
									 counter,
									 dof.data(),
									 rhsContribution.data() );

    //std::cout << "We are after prescribe RHS due to zero out matrix00 entries." << std::endl;



/*    {
      string filename = "after_matrix00_" + std::to_string( time ) + "_" + std::to_string( newtonIter ) + ".mtx";
      m_solidSolver->getSystemMatrix().write( filename, LAIOutputFormat::MATRIX_MARKET );
      GEOSX_LOG_RANK_0( "matrix00: written to " << filename );
    }
    {
      string filename = "after_residual0_" + std::to_string( time ) + "_" + std::to_string( newtonIter ) + ".mtx";
      m_solidSolver->getSystemRhs().write( filename, LAIOutputFormat::MATRIX_MARKET );
      GEOSX_LOG_RANK_0( "residual0: written to " << filename );
    }
*/
  }  //   if (!nodesWithAssignedDisp.empty())

  //TJ: matrix open and close have to be outside the   if (!nodesWithAssignedDisp.empty())
  m_solidSolver->getSystemMatrix().close();
  {
/*
    int const rank = MpiWrapper::Comm_rank( MPI_COMM_WORLD );
    std::cout << "Rank " << rank << " passed close system matrix." << std::endl;
*/
  }

  m_solidSolver->getSystemRhs().close();
  {
/*
    int const rank = MpiWrapper::Comm_rank( MPI_COMM_WORLD );
    std::cout << "Rank " << rank << " passed close system RHS." << std::endl;
*/
  }

  {
/*
    int const rank = MpiWrapper::Comm_rank( MPI_COMM_WORLD );
    std::cout << "Rank " << rank << " passed extra B.C.s applied at matrix_00." << std::endl;
*/
  }


  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager->getReference< globalIndex_array >( dispDofKey );
  arrayView1d< integer const > const & nodeGhostRank = nodeManager->ghostRank();

  m_matrix01.open();
  fsManager.Apply( time + dt,
                   domain,
                   "nodeManager",
                   keys::TotalDisplacement,
                   [&]( FieldSpecificationBase const * const bc,
                        string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group * const,
                        string const & )
  {
    SortedArray< localIndex > localSet;
    for( auto const & a : targetSet )
    {
      if( nodeGhostRank[a]<0 )
      {
        localSet.insert( a );
      }
    }
    bc->ZeroSystemRowsForBoundaryCondition< LAInterface >( localSet.toViewConst(),
                                                           dispDofNumber,
                                                           m_matrix01 );
  } );
  m_matrix01.close();


  //TJ: matrix open and close have to be outside the   if (!nodesWithAssignedDisp.empty())
  m_matrix01.open();
  if (!nodesWithAssignedDisp.empty())
  {
/*    integer newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;

    {
      string filename = "before_matrix01_" + std::to_string( time ) + "_" + std::to_string( newtonIter ) + ".mtx";
      m_matrix01.write( filename, LAIOutputFormat::MATRIX_MARKET );
      GEOSX_LOG_RANK_0( "matrix01: written to " << filename );
    }
*/
    integer const component = 0;
    arrayView1d< globalIndex const > const & dofMap
                                        = nodeManager->getReference< array1d< globalIndex > >( dispDofKey );
    for(auto a : nodesWithAssignedDisp)
    {
      globalIndex const dof = dofMap[a]+component;
      m_matrix01.clearRow(dof);
    }

    {
/*      string filename = "after_matrix01_" + std::to_string( time ) + "_" + std::to_string( newtonIter ) + ".mtx";
      m_matrix01.write( filename, LAIOutputFormat::MATRIX_MARKET );
      GEOSX_LOG_RANK_0( "matrix01: written to " << filename );
*/
    }
  }  //    if (!nodesWithAssignedDisp.empty())
  m_matrix01.close();


  {
/*
    int const rank = MpiWrapper::Comm_rank( MPI_COMM_WORLD );
    std::cout << "Rank " << rank << " passed extra B.C.s applied at matrix_01." << std::endl;
*/
  }

  m_flowSolver->ApplyBoundaryConditions( time,
                                         dt,
                                         domain,
                                         m_flowSolver->getDofManager(),
                                         m_flowSolver->getSystemMatrix(),
                                         m_flowSolver->getSystemRhs() );

  string const presDofKey = m_flowSolver->getDofManager().getKey( FlowSolverBase::viewKeyStruct::pressureString );

  m_matrix10.open();
  fsManager.Apply( time + dt,
                   domain,
                   "ElementRegions",
                   FlowSolverBase::viewKeyStruct::pressureString,
                   [&]( FieldSpecificationBase const * const fs,
                        string const &,
                        SortedArrayView< localIndex const > const & lset,
                        Group * subRegion,
                        string const & ) -> void
  {
    arrayView1d< globalIndex const > const &
    dofNumber = subRegion->getReference< array1d< globalIndex > >( presDofKey );
    arrayView1d< integer const > const & ghostRank = subRegion->group_cast< ObjectManagerBase * >()->ghostRank();

    SortedArray< localIndex > localSet;
    for( auto const & a : lset )
    {
      if( ghostRank[a]<0 )
      {
        localSet.insert( a );
      }
    }

    fs->ZeroSystemRowsForBoundaryCondition< LAInterface >( localSet.toViewConst(),
                                                           dofNumber,
                                                           m_matrix10 );
  } );
  m_matrix10.close();

  // debugging info.  can be trimmed once everything is working.
  if( getLogLevel()==2 )
  {
    // Before outputting anything generate permuation matrix and permute.
//    ElementRegionManager * const elemManager = mesh->getElemManager();

//    LAIHelperFunctions::CreatePermutationMatrix(nodeManager,
//                                                m_solidSolver->getSystemMatrix().numGlobalRows(),
//                                                m_solidSolver->getSystemMatrix().numGlobalCols(),
//                                                3,
//                                                m_solidSolver->getDofManager().getKey( keys::TotalDisplacement ),
//                                                m_permutationMatrix0);
//
//    LAIHelperFunctions::CreatePermutationMatrix(elemManager,
//                                                m_flowSolver->getSystemMatrix().numGlobalRows(),
//                                                m_flowSolver->getSystemMatrix().numGlobalCols(),
//                                                1,
//                                                m_flowSolver->getDofManager().getKey(
// FlowSolverBase::viewKeyStruct::pressureString ),
//                                                m_permutationMatrix1);

//    GEOSX_LOG_RANK_0("***********************************************************");
//    GEOSX_LOG_RANK_0("matrix00");
//    GEOSX_LOG_RANK_0("***********************************************************");
//    LAIHelperFunctions::PrintPermutedMatrix(m_solidSolver->getSystemMatrix(), m_permutationMatrix0, std::cout);
//    m_solidSolver->getSystemMatrix().print(std::cout);
    MpiWrapper::Barrier();

    GEOSX_LOG_RANK_0( "***********************************************************" );
    GEOSX_LOG_RANK_0( "matrix01" );
    GEOSX_LOG_RANK_0( "***********************************************************" );
//    LAIHelperFunctions::PrintPermutedMatrix(m_matrix01, m_permutationMatrix0, m_permutationMatrix1, std::cout);
    m_matrix01.print( std::cout );
    MpiWrapper::Barrier();

    GEOSX_LOG_RANK_0( "***********************************************************" );
    GEOSX_LOG_RANK_0( "matrix10" );
    GEOSX_LOG_RANK_0( "***********************************************************" );
//    LAIHelperFunctions::PrintPermutedMatrix(m_matrix10, m_permutationMatrix1, m_permutationMatrix0, std::cout);
    m_matrix10.print( std::cout );
    MpiWrapper::Barrier();

    GEOSX_LOG_RANK_0( "***********************************************************" );
    GEOSX_LOG_RANK_0( "matrix11" );
    GEOSX_LOG_RANK_0( "***********************************************************" );
//    LAIHelperFunctions::PrintPermutedMatrix(m_flowSolver->getSystemMatrix(), m_permutationMatrix1, std::cout);
    m_flowSolver->getSystemMatrix().print( std::cout );
    MpiWrapper::Barrier();

//    GEOSX_LOG_RANK_0("***********************************************************");
//    GEOSX_LOG_RANK_0("residual0");
//    GEOSX_LOG_RANK_0("***********************************************************");
//    LAIHelperFunctions::PrintPermutedVector(m_solidSolver->getSystemRhs(), m_permutationMatrix0, std::cout);
//    m_solidSolver->getSystemRhs().print(std::cout);
    MpiWrapper::Barrier();

    GEOSX_LOG_RANK_0( "***********************************************************" );
    GEOSX_LOG_RANK_0( "residual1" );
    GEOSX_LOG_RANK_0( "***********************************************************" );
//    LAIHelperFunctions::PrintPermutedVector(m_flowSolver->getSystemRhs(), m_permutationMatrix1, std::cout);
    m_flowSolver->getSystemRhs().print( std::cout );
    MpiWrapper::Barrier();
  }

  if( getLogLevel() >= 10 )
  {
    integer newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;

    {
      string filename = "matrix00_" + std::to_string( time ) + "_" + std::to_string( newtonIter ) + ".mtx";
      m_solidSolver->getSystemMatrix().write( filename, LAIOutputFormat::MATRIX_MARKET );
      GEOSX_LOG_RANK_0( "matrix00: written to " << filename );
    }
    {
      string filename = "matrix01_" + std::to_string( time ) + "_" + std::to_string( newtonIter ) + ".mtx";
      m_matrix01.write( filename, LAIOutputFormat::MATRIX_MARKET );
      GEOSX_LOG_RANK_0( "matrix01: written to " << filename );
    }
    {
      string filename = "matrix10_" + std::to_string( time ) + "_" + std::to_string( newtonIter ) + ".mtx";
      m_matrix10.write( filename, LAIOutputFormat::MATRIX_MARKET );
      GEOSX_LOG_RANK_0( "matrix10: written to " << filename );
    }
    {
      string filename = "matrix11_" + std::to_string( time ) + "_" + std::to_string( newtonIter ) + ".mtx";
      m_flowSolver->getSystemMatrix().write( filename, LAIOutputFormat::MATRIX_MARKET );
      GEOSX_LOG_RANK_0( "matrix11: written to " << filename );
    }
    {
      string filename = "residual0_" + std::to_string( time ) + "_" + std::to_string( newtonIter ) + ".mtx";
      m_solidSolver->getSystemRhs().write( filename, LAIOutputFormat::MATRIX_MARKET );
      GEOSX_LOG_RANK_0( "residual0: written to " << filename );
    }
    {
      string filename = "residual1_" + std::to_string( time ) + "_" + std::to_string( newtonIter ) + ".mtx";
      m_flowSolver->getSystemRhs().write( filename, LAIOutputFormat::MATRIX_MARKET );
      GEOSX_LOG_RANK_0( "residual1: written to " << filename );
    }
  }
}

real64
HydrofractureSolver::
  CalculateResidualNorm( DomainPartition const * const domain,
                         DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                         ParallelVector const & GEOSX_UNUSED_PARAM( rhs ) )
{
  GEOSX_MARK_FUNCTION;

  /*
     real64 const fluidResidual = m_flowSolver->getSystemRhs().norm2();
     real64 const solidResidual = m_solidSolver->getSystemRhs().norm2();
   */

  real64 const fluidResidual = m_flowSolver->CalculateResidualNorm( domain,
                                                                    m_flowSolver->getDofManager(),
                                                                    m_flowSolver->getSystemRhs() );

  real64 const solidResidual = m_solidSolver->CalculateResidualNorm( domain,
                                                                     m_solidSolver->getDofManager(),
                                                                     m_solidSolver->getSystemRhs() );

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    char output[200] = {0};
    sprintf( output,
             "( Rfluid, Rsolid ) = (%4.2e, %4.2e) ; ",
             fluidResidual,
             solidResidual );
    std::cout<<output;
  }

  return fluidResidual + solidResidual;
}



void
HydrofractureSolver::
  AssembleForceResidualDerivativeWrtPressure( DomainPartition * const domain,
                                              ParallelMatrix * const matrix01,
                                              ParallelVector * const rhs0 )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );

  FaceManager const * const faceManager = mesh->getFaceManager();
  NodeManager * const nodeManager = mesh->getNodeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

  arrayView1d< R1Tensor const > const & faceNormal = faceManager->faceNormal();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager->nodeList().toViewConst();

  arrayView1d< R1Tensor > const &
  fext = nodeManager->getReference< array1d< R1Tensor > >( SolidMechanicsLagrangianFEM::viewKeyStruct::forceExternal );
  fext = {0, 0, 0};

  string const presDofKey = m_flowSolver->getDofManager().getKey( FlowSolverBase::viewKeyStruct::pressureString );
  string const dispDofKey = m_solidSolver->getDofManager().getKey( keys::TotalDisplacement );

  arrayView1d< globalIndex > const &
  dispDofNumber =  nodeManager->getReference< globalIndex_array >( dispDofKey );

  matrix01->open();
  rhs0->open();

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {

    arrayView1d< globalIndex > const &
    faceElementDofNumber = subRegion.getReference< array1d< globalIndex > >( presDofKey );

    if( subRegion.hasWrapper( "pressure" ) )
    {
      arrayView1d< real64 const > const & fluidPressure = subRegion.getReference< array1d< real64 > >( "pressure" );
      arrayView1d< real64 const > const & deltaFluidPressure = subRegion.getReference< array1d< real64 > >( "deltaPressure" );
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView1d< real64 > const & area = subRegion.getElementArea();
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();

      forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        R1Tensor Nbar = faceNormal[elemsToFaces[kfe][0]];
        Nbar -= faceNormal[elemsToFaces[kfe][1]];
        Nbar.Normalize();

        localIndex const kf0 = elemsToFaces[kfe][0];
        localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );

        globalIndex rowDOF[24];
        real64 nodeRHS[24];
        stackArray2d< real64, 12*12 > dRdP( numNodesPerFace*3, 1 );
        globalIndex colDOF = faceElementDofNumber[kfe];

        real64 const Ja = area[kfe] / numNodesPerFace;

        //          std::cout<<"fluidPressure["<<kfe<<"] = "<<fluidPressure[kfe]+deltaFluidPressure[kfe]<<std::endl;
        real64 nodalForceMag = ( fluidPressure[kfe]+deltaFluidPressure[kfe] ) * Ja;
        R1Tensor nodalForce( Nbar );
        nodalForce *= nodalForceMag;

        //          std::cout << "    rank " << MpiWrapper::Comm_rank(MPI_COMM_GEOSX) << ", faceElement " << kfe <<
        // std::endl;
        //          std::cout << "    fluid pressure " << fluidPressure[kfe]+deltaFluidPressure[kfe] << std::endl;
        //          std::cout << "    nodalForce " << nodalForce << std::endl;
        for( localIndex kf=0; kf<2; ++kf )
        {
          localIndex const faceIndex = elemsToFaces[kfe][kf];


          for( localIndex a=0; a<numNodesPerFace; ++a )
          {

            for( int i=0; i<3; ++i )
            {
              rowDOF[3*a+i] = dispDofNumber[faceToNodeMap( faceIndex, a )] + i;
              nodeRHS[3*a+i] = -nodalForce[i] * pow( -1, kf );
              fext[faceToNodeMap( faceIndex, a )][i] += -nodalForce[i] * pow( -1, kf );

              dRdP( 3*a+i, 0 ) = -Ja * Nbar[i] * pow( -1, kf );
              // this is for debugging
              //                if (dispDofNumber[faceToNodeMap(faceIndex, a)] == 0 ||
              // dispDofNumber[faceToNodeMap(faceIndex, a)] == 6 || dispDofNumber[faceToNodeMap(faceIndex, a)] == 12 ||
              // dispDofNumber[faceToNodeMap(faceIndex, a)] == 18)
              //                  std::cout << "rank " << MpiWrapper::Comm_rank(MPI_COMM_GEOSX) << "DOF index " <<
              // dispDofNumber[faceToNodeMap(faceIndex, a)] + i << " contribution " << nodeRHS[3*a+i] << std::endl;

            }
          }
          if( ghostRank[kfe] < 0 )
          {

            rhs0->add( rowDOF,
                       nodeRHS,
                       numNodesPerFace*3 );


            matrix01->add( rowDOF,
                           &colDOF,
                           dRdP.data(),
                           numNodesPerFace * 3,
                           1 );
          }
        }
      } );
    }
  } );

  matrix01->close();
  rhs0->close();
}

void
HydrofractureSolver::
  AssembleFluidMassResidualDerivativeWrtDisplacement( DomainPartition const * const domain,
                                                      ParallelMatrix * const matrix10,
                                                      ParallelVector * const GEOSX_UNUSED_PARAM( rhs0 ) )
{
  GEOSX_MARK_FUNCTION;

  //TJ the hard coded elmt length
  real64 const meshSize = m_meshSize;  // this value needs to be changed for a mesh-refinement
  //GEOSX_LOG_RANK_0( "Mesh size = " << meshSize );

  MeshLevel const & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  FaceManager const * const faceManager = mesh.getFaceManager();
  NodeManager const * const nodeManager = mesh.getNodeManager();
  ConstitutiveManager const * const constitutiveManager = domain->getConstitutiveManager();

  string const presDofKey = m_flowSolver->getDofManager().getKey( FlowSolverBase::viewKeyStruct::pressureString );
  string const dispDofKey = m_solidSolver->getDofManager().getKey( keys::TotalDisplacement );

  CRSMatrixView< real64 const, localIndex const > const &
  dFluxResidual_dAperture = m_flowSolver->getDerivativeFluxResidual_dAperture().toViewConst();

  ContactRelationBase const * const
  contactRelation = constitutiveManager->GetGroup< ContactRelationBase >( m_contactRelationName );

  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & disp = nodeManager->totalDisplacement();
  SurfaceGenerator * const mySurface = this->getParent()->GetGroup< SurfaceGenerator >( "SurfaceGen" );

  real64 const shearModulus = domain->GetGroup("Constitutive")
			                  ->GetGroup("rock")
			                  ->getReference<real64>("defaultShearModulus");
  real64 const bulkModulus = domain->GetGroup("Constitutive")
		                         ->GetGroup("rock")
				         ->getReference<real64>("defaultBulkModulus");
  real64 const toughness = mySurface->getReference<real64>("rockToughness");
  real64 const viscosity = domain->GetGroup("Constitutive")
                                 ->GetGroup("water")
				       ->getReference<real64>("defaultViscosity");


  // The unit of injectionRate is kg per second
  real64 const injectionRate = domain->getParent()
                                     ->GetGroup<FieldSpecificationManager>("FieldSpecifications")
                                     ->GetGroup<SourceFluxBoundaryCondition>("sourceTerm")
					   ->getReference<real64>("scale");

  // The injectionRate is only for half domain of the KGD problem,
  // to retrieve the full injection rate, we need to multiply it by 2.0
  real64 const q0 = 2.0 * std::abs(injectionRate) /1.0e3;
  real64 const total_time = m_totalTime;

  real64 const nu = ( 1.5 * bulkModulus - shearModulus ) / ( 3.0 * bulkModulus + shearModulus );
  real64 const E = ( 9.0 * bulkModulus * shearModulus )/ ( 3.0 * bulkModulus + shearModulus );
  real64 const Eprime = E/(1.0-nu*nu);
  real64 const PI = 2 * acos(0.0);
  real64 const Kprime = 4.0*sqrt(2.0/PI)*toughness;
  real64 const mup = 12.0 * viscosity;




  matrix10->open();

  forTargetSubRegionsComplete< FaceElementSubRegion >( mesh,
                                                       [&]( localIndex const,
                                                            localIndex const,
                                                            localIndex const,
                                                            ElementRegionBase const & region,
                                                            FaceElementSubRegion const & subRegion )
  {
    string const
    & fluidName = m_flowSolver->fluidModelNames()[m_flowSolver->targetRegionIndex( region.getName() )];
    SingleFluidBase const & fluid = GetConstitutiveModel< SingleFluidBase >( subRegion, fluidName );

    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
    arrayView1d< globalIndex const > const
    & presDofNumber = subRegion.getReference< array1d< globalIndex > >( presDofKey );
    arrayView1d< globalIndex const > const
    & dispDofNumber = nodeManager->getReference< array1d< globalIndex > >( dispDofKey );

    arrayView2d< real64 const > const & dens = fluid.density();

    arrayView1d< real64 const > const & aperture = subRegion.getElementAperture();
    arrayView1d< real64 const > const & area = subRegion.getElementArea();

    arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();
    ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager->nodeList().toViewConst();

    arrayView1d< R1Tensor const > const & faceNormal = faceManager->faceNormal();



//    arrayView1d< real64 const > const & separationCoeff = subRegion.getSeparationCoefficient();
//    arrayView1d<real64 const> const &
//    dseparationCoeff_dAper  =
// subRegion.getReference<array1d<real64>>(FaceElementSubRegion::viewKeyStruct::dSeparationCoeffdAperString);

    //TJ
//    FaceElementSubRegion::NodeMapType & nodeMap = subRegion.nodeList();


    forAll< serialPolicy >( subRegion.size(), [=]( localIndex ei )
    {
      //if (elemGhostRank[ei] < 0)
      {
        globalIndex const elemDOF = presDofNumber[ei];
        localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[ei][0] );
        real64 const
        dAccumulationResidualdAperture = dens[ei][0] * area[ei];                                                        //*
                                                                                                                        // (
                                                                                                                        // separationCoeff[ei]
                                                                                                                        // +
        //   aperture[ei] * dseparationCoeff_dAper[ei] );


        globalIndex nodeDOF[8 * 3];

        R1Tensor Nbar = faceNormal[elemsToFaces[ei][0]];
        Nbar -= faceNormal[elemsToFaces[ei][1]];
        Nbar.Normalize();

        stackArray1d< real64, 24 > dRdU( 2 * numNodesPerFace * 3 );

        // Accumulation derivative
        if( elemGhostRank[ei] < 0 )
        {
          //GEOS_LOG_RANK( "dAccumulationResidualdAperture("<<ei<<") = "<<dAccumulationResidualdAperture );
          for( localIndex kf = 0; kf < 2; ++kf )
          {
            for( localIndex a = 0; a < numNodesPerFace; ++a )
            {
              for( int i = 0; i < 3; ++i )
              {
                nodeDOF[kf * 3 * numNodesPerFace + 3 * a + i] =
                  dispDofNumber[faceToNodeMap( elemsToFaces[ei][kf], a )] + i;
                real64 const dGap_dU = -pow( -1, kf ) * Nbar[i] / numNodesPerFace;
                real64 const dAper_dU = contactRelation->dEffectiveAperture_dAperture( aperture[ei] ) * dGap_dU;
                dRdU( kf * 3 * numNodesPerFace + 3 * a + i ) = dAccumulationResidualdAperture * dAper_dU;
              }
            }
          }

          //TJ: for the tip element, the relationship between aperture and displacement
          //    is nonlinear. Therefore, we need to rewrite dAper_dU and dRdU
          SortedArray< localIndex > const trailingFaces = mySurface->getTrailingFaces();
	  SortedArray< localIndex > const tipNodes = mySurface->getTipNodes();

          for(auto const & trailingFace : trailingFaces)
          {
            //TJ: elmt ei is a tip element, do the modification on dRdU
            if ( (elemsToFaces[ei][0] == trailingFace) || (elemsToFaces[ei][1] == trailingFace) )
            {
//              std::cout << "elmt " << ei << " is a tip element." << std::endl;
              //TJ: calculate 1/2 * ( (u11-u10) + (u21-u20) ) dot normal_direction
              real64 averageGap = 0.0;
              for (localIndex kf = 0; kf < 2; ++kf)
              {
                for( localIndex a = 0; a < numNodesPerFace; ++a )
                {
                  localIndex node = faceToNodeMap( elemsToFaces[ei][kf], a );
  	          if ( std::find( tipNodes.begin(), tipNodes.end(), node ) == tipNodes.end() )
  	          {
 //                   std::cout << "Face " << kf << ": node " << node << std::endl;
                    R1Tensor temp = disp[node];
                    averageGap += (-pow(-1,kf)) * Dot( temp, Nbar)/2 ;
  	          }
                }
              }
//              std::cout << averageGap << std::endl;

              real64 coeff;
              real64 dGap_dU_tip;
              //real64 dAper_dU_tip;
              if (viscosity < 2.0e-3) // Toughness-dominated case
              {
        	//TJ: the tip asymptote w = Kprime / Eprime * x^(1/2)
        	coeff = 2.0/3.0 * pow(Eprime/Kprime, 2.0)/meshSize;
                for (localIndex kf = 0; kf < 2; ++kf)
                {
                  for( localIndex a = 0; a < numNodesPerFace; ++a )
                  {
                    localIndex node = faceToNodeMap( elemsToFaces[ei][kf], a );
  	            if ( std::find( tipNodes.begin(), tipNodes.end(), node ) == tipNodes.end() )
  	            {
  	              for (int i = 0; i < 3; ++i)
  	              {
  	                dGap_dU_tip = 3.0/2.0 * coeff * pow(averageGap, 2.0) * (-pow(-1, kf)) *Nbar[i];
  	                //dAper_dU_tip = contactRelation->dEffectiveAperture_dAperture( aperture[ei] ) * dGap_dU_tip;
  	                dRdU( kf * 3 * numNodesPerFace + 3 * a + i ) = dAccumulationResidualdAperture * dGap_dU_tip;
  	              }
  	            }
  	            else
  	            {
  	              for (int i = 0; i < 3; ++i)
  	              {
  	                dRdU( kf * 3 * numNodesPerFace + 3 * a + i ) = 0.0;
  	              }
  	            } // if else
                  } // for( localIndex a = 0; a < numNodesPerFace; ++a )
                } // for (localIndex kf = 0; kf < 2; ++kf)
              }
              else // Viscosity-dominated case
              {
        	real64 Lm = pow( Eprime*pow(q0,3.0)*pow(total_time,4.0)/mup, 1.0/6.0 );
        	real64 gamma_m0 = 0.616;
        	real64 velocity = 2.0/3.0 * Lm * gamma_m0 / total_time;
        	real64 Betam = pow(2.0, 1.0/3.0) * pow(3.0, 5.0/6.0);
        	coeff = 3.0/5.0 * pow(Betam, -3.0/2.0) * pow(Eprime/mup/velocity, 0.5)/meshSize;
                for (localIndex kf = 0; kf < 2; ++kf)
                {
                  for( localIndex a = 0; a < numNodesPerFace; ++a )
                  {
                    localIndex node = faceToNodeMap( elemsToFaces[ei][kf], a );
  	            if ( std::find( tipNodes.begin(), tipNodes.end(), node ) == tipNodes.end() )
  	            {
  	              for (int i = 0; i < 3; ++i)
  	              {
  	                dGap_dU_tip = 5.0/4.0 * coeff * pow(averageGap, 3.0/2.0) * (-pow(-1, kf)) * Nbar[i];
  	                //dAper_dU_tip = contactRelation->dEffectiveAperture_dAperture( aperture[ei] ) * dGap_dU_tip;
  	                dRdU( kf * 3 * numNodesPerFace + 3 * a + i ) = dAccumulationResidualdAperture * dGap_dU_tip;
  	              }
  	            }
  	            else
  	            {
  	              for (int i = 0; i < 3; ++i)
  	              {
  	                dRdU( kf * 3 * numNodesPerFace + 3 * a + i ) = 0.0;
  	              }
  	            } // if else
                  } // for( localIndex a = 0; a < numNodesPerFace; ++a )
                } // for (localIndex kf = 0; kf < 2; ++kf)
              } // toughness or viscosity dominated case
            }  // if ei is a tip element
          } // loop over all the trailing faces

          matrix10->add( elemDOF,
                         nodeDOF,
                         dRdU.data(),
                         2 * numNodesPerFace * 3 );
        }

        // flux derivative
        localIndex const numColumns = dFluxResidual_dAperture.numNonZeros( ei );
        arraySlice1d< localIndex const > const & columns = dFluxResidual_dAperture.getColumns( ei );
        arraySlice1d< real64 const > const & values = dFluxResidual_dAperture.getEntries( ei );

        for( localIndex kfe2 = 0; kfe2 < numColumns; ++kfe2 )
        {
          real64 dRdAper = values[kfe2];
          localIndex const ei2 = columns[kfe2];
//          GEOS_LOG_RANK( "dRdAper("<<ei<<", "<<ei2<<") = "<<dRdAper );

          for( localIndex kf = 0; kf < 2; ++kf )
          {
            for( localIndex a = 0; a < numNodesPerFace; ++a )
            {
              for( int i = 0; i < 3; ++i )
              {
                nodeDOF[kf * 3 * numNodesPerFace + 3 * a + i] =
                  dispDofNumber[faceToNodeMap( elemsToFaces[ei2][kf], a )] + i;
                real64 const dGap_dU = -pow( -1, kf ) * Nbar[i] / numNodesPerFace;
                real64 const
                dAper_dU = contactRelation->dEffectiveAperture_dAperture( aperture[ei2] ) * dGap_dU;
                dRdU( kf * 3 * numNodesPerFace + 3 * a + i ) = dRdAper * dAper_dU;
              }
            }
          }

          //TJ: for the tip element, the relationship between aperture and displacement
          //    is nonlinear. Therefore, we need to rewrite dAper_dU and dRdU
          SortedArray< localIndex > const trailingFaces = mySurface->getTrailingFaces();
	  SortedArray< localIndex > const tipNodes = mySurface->getTipNodes();

          for(auto const & trailingFace : trailingFaces)
          {
            //TJ: elmt ei is a tip element, do the modification on dRdU
            if ( (elemsToFaces[ei2][0] == trailingFace) || (elemsToFaces[ei2][1] == trailingFace) )
            {
//              std::cout << "elmt " << ei2 << " is a tip element." << std::endl;
              //TJ: calculate 1/2 * ( (u11-u10) + (u21-u20) ) dot normal_direction
              real64 averageGap = 0.0;
              for (localIndex kf = 0; kf < 2; ++kf)
              {
                for( localIndex a = 0; a < numNodesPerFace; ++a )
                {
                  localIndex node = faceToNodeMap( elemsToFaces[ei2][kf], a );
  	          if ( std::find( tipNodes.begin(), tipNodes.end(), node ) == tipNodes.end() )
  	          {
//                    std::cout << "Face " << kf << ": node " << node << std::endl;
                    R1Tensor temp = disp[node];
                    averageGap += (-pow(-1,kf)) * Dot( temp, Nbar)/2 ;
  	          }
                }
              }
//              std::cout << averageGap << std::endl;

              real64 coeff;
              real64 dGap_dU_tip;
              //real64 dAper_dU_tip;
              if (viscosity < 2.0e-3) // Toughness-dominated case
              {
        	//TJ: the tip asymptote w = Kprime / Eprime * x^(1/2)
        	coeff = 2.0/3.0 * pow(Eprime/Kprime, 2.0)/meshSize;
                for (localIndex kf = 0; kf < 2; ++kf)
                {
                  for( localIndex a = 0; a < numNodesPerFace; ++a )
                  {
                    localIndex node = faceToNodeMap( elemsToFaces[ei2][kf], a );
  	            if ( std::find( tipNodes.begin(), tipNodes.end(), node ) == tipNodes.end() )
  	            {
  	              for (int i = 0; i < 3; ++i)
  	              {
  	                dGap_dU_tip = 3.0/2.0 * coeff * pow(averageGap, 2.0) * (-pow(-1, kf)) *Nbar[i];
  	                dGap_dU_tip = 2.0/3.0 * (-pow(-1, kf)) *Nbar[i]/2;
  	                //dAper_dU_tip = contactRelation->dEffectiveAperture_dAperture( aperture[ei2] ) * dGap_dU_tip;
  	                dRdU( kf * 3 * numNodesPerFace + 3 * a + i ) = dRdAper * dGap_dU_tip;
  	              }
  	            }
  	            else
  	            {
  	              for (int i = 0; i < 3; ++i)
  	              {
  	                dRdU( kf * 3 * numNodesPerFace + 3 * a + i ) = 0.0;
  	              }
  	            } // if else
                  } // for( localIndex a = 0; a < numNodesPerFace; ++a )
                } // for (localIndex kf = 0; kf < 2; ++kf)
              }
              else // Viscosity-dominated case
              {
        	real64 Lm = pow( Eprime*pow(q0,3.0)*pow(total_time,4.0)/mup, 1.0/6.0 );
        	real64 gamma_m0 = 0.616;
        	real64 velocity = 2.0/3.0 * Lm * gamma_m0 / total_time;
        	real64 Betam = pow(2.0, 1.0/3.0) * pow(3.0, 5.0/6.0);
        	coeff = 3.0/5.0 * pow(Betam, -3.0/2.0) * pow(Eprime/mup/velocity, 0.5)/meshSize;
                for (localIndex kf = 0; kf < 2; ++kf)
                {
                  for( localIndex a = 0; a < numNodesPerFace; ++a )
                  {
                    localIndex node = faceToNodeMap( elemsToFaces[ei2][kf], a );
  	            if ( std::find( tipNodes.begin(), tipNodes.end(), node ) == tipNodes.end() )
  	            {
  	              for (int i = 0; i < 3; ++i)
  	              {
  	                dGap_dU_tip = 5.0/4.0 * coeff * pow(averageGap, 3.0/2.0) * (-pow(-1, kf)) * Nbar[i];
  	                dGap_dU_tip = 3.0/5.0 * (-pow(-1, kf)) *Nbar[i]/2;
  	                //dAper_dU_tip = contactRelation->dEffectiveAperture_dAperture( aperture[ei2] ) * dGap_dU_tip;
  	                dRdU( kf * 3 * numNodesPerFace + 3 * a + i ) = dRdAper * dGap_dU_tip;
  	              }
  	            }
  	            else
  	            {
  	              for (int i = 0; i < 3; ++i)
  	              {
  	                dRdU( kf * 3 * numNodesPerFace + 3 * a + i ) = 0.0;
  	              }
  	            } // if else
                  } // for( localIndex a = 0; a < numNodesPerFace; ++a )
                } // for (localIndex kf = 0; kf < 2; ++kf)
              } // toughness or viscosity dominated case
            }  // if ei2 is a tip element
          } // loop over all the trailing faces

          matrix10->add( elemDOF,
                         nodeDOF,
                         dRdU.data(),
                         2 * numNodesPerFace * 3 );

        }
      }
    } );
  } );

  matrix10->close();
}

void
HydrofractureSolver::
  ApplySystemSolution( DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                       ParallelVector const & GEOSX_UNUSED_PARAM( solution ),
                       real64 const scalingFactor,
                       DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->ApplySystemSolution( m_solidSolver->getDofManager(),
                                      m_solidSolver->getSystemSolution(),
                                      scalingFactor,
                                      domain );
  m_flowSolver->ApplySystemSolution( m_flowSolver->getDofManager(),
                                     m_flowSolver->getSystemSolution(),
                                     -scalingFactor,
                                     domain );

  this->UpdateDeformationForCoupling( domain );
}

}

#ifdef USING_TRILINOS

#include "Epetra_FEVector.h"
#include "Epetra_FECrsMatrix.h"
#include "EpetraExt_MatrixMatrix.h"
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
#include "Teuchos_Time.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

#endif

namespace geosx
{

void HydrofractureSolver::SolveSystem( DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                       ParallelMatrix &,
                                       ParallelVector &,
                                       ParallelVector & )
{
  GEOSX_MARK_FUNCTION;

#ifdef USING_TRILINOS

  /*
     globalIndex numU = m_solidSolver->getSystemRhs().globalSize();
     globalIndex numP = m_flowSolver->getSystemRhs().globalSize();
     GEOSX_LOG_RANK_0("size = " << numU << " + " << numP);
   */

  integer const newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;

  using namespace Teuchos;
  using namespace Thyra;

  Teuchos::Time clock( "solveClock" );

  GEOSX_MARK_BEGIN( Setup );
  Epetra_FECrsMatrix * p_matrix[2][2];
  Epetra_FEVector * p_rhs[2];
  Epetra_FEVector * p_solution[2];

  p_rhs[0] = &m_solidSolver->getSystemRhs().unwrapped();
  p_rhs[1] = &m_flowSolver->getSystemRhs().unwrapped();

  p_solution[0] = &m_solidSolver->getSystemSolution().unwrapped();
  p_solution[1] = &m_flowSolver->getSystemSolution().unwrapped();

  p_matrix[0][0] = &m_solidSolver->getSystemMatrix().unwrapped();
  p_matrix[0][1] = &m_matrix01.unwrapped();
  p_matrix[1][0] = &m_matrix10.unwrapped();
  p_matrix[1][1] = &m_flowSolver->getSystemMatrix().unwrapped();

  // scale and symmetrize

  m_densityScaling = 1e-3;
  m_pressureScaling = 1e9;

  p_matrix[0][1]->Scale( m_pressureScaling );
  p_matrix[1][0]->Scale( m_pressureScaling*m_densityScaling );
  p_matrix[1][1]->Scale( m_pressureScaling*m_pressureScaling*m_densityScaling );
  p_rhs[1]->Scale( m_pressureScaling*m_densityScaling );

  // SCHEME CHOICES
  //
  // there are several flags to control solver behavior.
  // these should be compared in a scaling study.
  //
  // -- whether to use a block diagonal or a
  //    block triangular preconditioner.
  // -- whether to use BiCGstab or GMRES for the
  //    krylov solver.  GMRES is generally more robust,
  //    BiCGstab sometimes shows better parallel performance.
  //    false is probably better.

  const bool use_diagonal_prec = true;
  const bool use_bicgstab      = (m_linearSolverParameters.solverType == "bicgstab");

  // set initial guess to zero

  p_solution[0]->PutScalar( 0.0 );
  p_solution[1]->PutScalar( 0.0 );

  // create separate displacement component matrix

  clock.start( true );
  if( newtonIter==0 )
  {
    m_blockDiagUU.reset( new ParallelMatrix());
    LAIHelperFunctions::SeparateComponentFilter< TrilinosInterface >( m_solidSolver->getSystemMatrix(), *m_blockDiagUU, 3 );
  }

  // create schur complement approximation matrix

  Epetra_CrsMatrix * schurApproxPP = NULL; // confirm we delete this at end of function!
  {
    Epetra_Vector diag( p_matrix[0][0]->RowMap());
    Epetra_Vector diagInv( p_matrix[0][0]->RowMap());

    p_matrix[0][0]->ExtractDiagonalCopy( diag );
    diagInv.Reciprocal( diag );

    Epetra_FECrsMatrix DB( *p_matrix[0][1] );
    DB.LeftScale( diagInv );
    DB.FillComplete();

    Epetra_FECrsMatrix BtDB( Epetra_DataAccess::Copy, p_matrix[1][1]->RowMap(), 1 );
    EpetraExt::MatrixMatrix::Multiply( *p_matrix[1][0], false, DB, false, BtDB );
    EpetraExt::MatrixMatrix::Add( BtDB, false, -1.0, *p_matrix[1][1], false, 1.0, schurApproxPP );

    schurApproxPP->FillComplete();
  }
  double auxTime = clock.stop();
  GEOSX_MARK_END( Setup );

  // we want to use thyra to wrap epetra operators and vectors
  // for individual blocks.  this is an ugly conversion, but
  // it is basically just window dressing.
  //
  // note the use of Teuchos::RCP reference counted pointers.
  // The general syntax is usually one of:
  //
  //   RCP<T> Tptr = rcp(new T)
  //   RCP<T> Tptr = nonMemberConstructor();
  //   RCP<T> Tptr (t_ptr,false)
  //
  // where "false" implies the RCP does not own the object and
  // should not attempt to delete it when finished.

  GEOSX_MARK_BEGIN( THYRA_SETUP );

  RCP< const Thyra::LinearOpBase< double > >  matrix_block[2][2];
  RCP< Thyra::MultiVectorBase< double > >     lhs_block[2];
  RCP< Thyra::MultiVectorBase< double > >     rhs_block[2];

  for( unsigned i=0; i<2; ++i )
    for( unsigned j=0; j<2; ++j )
    {
      RCP< Epetra_Operator > mmm ( &*p_matrix[i][j], false );
      matrix_block[i][j] = Thyra::epetraLinearOp( mmm );
    }

  RCP< Epetra_Operator > bbb( &m_blockDiagUU->unwrapped(), false );
  RCP< Epetra_Operator > ppp( schurApproxPP, false );

  RCP< const Thyra::LinearOpBase< double > >  blockDiagOp = Thyra::epetraLinearOp( bbb );
  RCP< const Thyra::LinearOpBase< double > >  schurOp = Thyra::epetraLinearOp( ppp );

  for( unsigned i=0; i<2; ++i )
  {
    RCP< Epetra_MultiVector > lll ( &*p_solution[i], false );
    RCP< Epetra_MultiVector > rrr ( &*p_rhs[i], false );

    lhs_block[i] = Thyra::create_MultiVector( lll, matrix_block[i][i]->domain());
    rhs_block[i] = Thyra::create_MultiVector( rrr, matrix_block[i][i]->range());
  }

  // now use thyra to create an operator representing
  // the full block 2x2 system

  RCP< const Thyra::LinearOpBase< double > > matrix = Thyra::block2x2( matrix_block[0][0],
                                                                       matrix_block[0][1],
                                                                       matrix_block[1][0],
                                                                       matrix_block[1][1] );

  // creating a representation of the blocked
  // rhs and lhs is a little uglier.

  RCP< Thyra::ProductMultiVectorBase< double > > rhs;
  {
    Teuchos::Array< RCP< Thyra::MultiVectorBase< double > > > mva;
    Teuchos::Array< RCP< const Thyra::VectorSpaceBase< double > > > mvs;

    for( unsigned i=0; i<2; ++i )
    {
      mva.push_back( rhs_block[i] );
      mvs.push_back( rhs_block[i]->range());
    }

    RCP< const Thyra::DefaultProductVectorSpace< double > > vs = Thyra::productVectorSpace< double >( mvs );

    rhs = Thyra::defaultProductMultiVector< double >( vs, mva );
  }

  RCP< Thyra::ProductMultiVectorBase< double > > lhs;

  {
    Teuchos::Array< RCP< Thyra::MultiVectorBase< double > > > mva;
    Teuchos::Array< RCP< const Thyra::VectorSpaceBase< double > > > mvs;

    for( unsigned i=0; i<2; ++i )
    {
      mva.push_back( lhs_block[i] );
      mvs.push_back( lhs_block[i]->range());
    }

    RCP< const Thyra::DefaultProductVectorSpace< double > > vs = Thyra::productVectorSpace< double >( mvs );

    lhs = Thyra::defaultProductMultiVector< double >( vs, mva );
  }

  GEOSX_MARK_END( THYRA_SETUP );

  // for the preconditioner, we need two approximate inverses,
  // we store both "sub operators" in a 1x2 array:

  RCP< const Thyra::LinearOpBase< double > > sub_op[2];

  clock.start( true );
  GEOSX_MARK_BEGIN( PRECONDITIONER );

  for( unsigned i=0; i<2; ++i ) // loop over diagonal blocks
  {
    RCP< Teuchos::ParameterList > list = rcp( new Teuchos::ParameterList( "precond_list" ), true );

    if( m_linearSolverParameters.preconditionerType == "amg" )
    {
      list->set( "Preconditioner Type", "ML" );
      list->sublist( "Preconditioner Types" ).sublist( "ML" ).set( "Base Method Defaults", "SA" );
      list->sublist( "Preconditioner Types" ).sublist( "ML" ).sublist( "ML Settings" ).set( "PDE equations", (i==0?3:1));
      list->sublist( "Preconditioner Types" ).sublist( "ML" ).sublist( "ML Settings" ).set( "ML output", 0 );
      list->sublist( "Preconditioner Types" ).sublist( "ML" ).sublist( "ML Settings" ).set( "aggregation: type", "Uncoupled" );
      list->sublist( "Preconditioner Types" ).sublist( "ML" ).sublist( "ML Settings" ).set( "aggregation: threshold", 1e-3 );

      if( i==0 ) // smoother for mechanics block
      {
        list->sublist( "Preconditioner Types" ).sublist( "ML" ).sublist( "ML Settings" ).set( "smoother: type", "Chebyshev" );
        list->sublist( "Preconditioner Types" ).sublist( "ML" ).sublist( "ML Settings" ).set( "smoother: sweeps", 3 );
        list->sublist( "Preconditioner Types" ).sublist( "ML" ).sublist( "ML Settings" ).set( "coarse: type", "Chebyshev" );
        list->sublist( "Preconditioner Types" ).sublist( "ML" ).sublist( "ML Settings" ).set( "coarse: sweeps", 3 );
      }
      else // smoother for flow block
      {
        list->sublist( "Preconditioner Types" ).sublist( "ML" ).sublist( "ML Settings" ).set( "smoother: type", "Chebyshev" );
        list->sublist( "Preconditioner Types" ).sublist( "ML" ).sublist( "ML Settings" ).set( "smoother: sweeps", 3 );
      }

    }
    else // use ILU for both blocks
    {
      list->set( "Preconditioner Type", "Ifpack" );
      list->sublist( "Preconditioner Types" ).sublist( "Ifpack" ).set( "Prec Type", "ILU" );
    }

    Stratimikos::DefaultLinearSolverBuilder builder;
    builder.setParameterList( list );

    RCP< const Thyra::PreconditionerFactoryBase< double > > strategy = createPreconditioningStrategy( builder );
    RCP< Thyra::PreconditionerBase< double > > tmp;

    if( i==0 )
      tmp = prec( *strategy, blockDiagOp );
    else
      tmp = prec( *strategy, schurOp );
    //tmp = prec(*strategy,matrix_block[i][i]);

    sub_op[i] = tmp->getUnspecifiedPrecOp();
  }


  // create zero operators for off diagonal blocks

  RCP< const Thyra::LinearOpBase< double > > zero_01
    = rcp( new Thyra::DefaultZeroLinearOp< double >( matrix_block[0][0]->range(),
                                                     matrix_block[1][1]->domain()));

  RCP< const Thyra::LinearOpBase< double > > zero_10
    = rcp( new Thyra::DefaultZeroLinearOp< double >( matrix_block[1][1]->range(),
                                                     matrix_block[0][0]->domain()));

  // now build the block preconditioner

  RCP< const Thyra::LinearOpBase< double > > preconditioner;

  if( use_diagonal_prec )
  {
    preconditioner = Thyra::block2x2( sub_op[0], zero_01, zero_10, sub_op[1] );
  }
  else
  {
    RCP< const Thyra::LinearOpBase< double > > eye_00
      = Teuchos::rcp( new Thyra::DefaultIdentityLinearOp< double >( matrix_block[0][0]->range()));

    RCP< const Thyra::LinearOpBase< double > > eye_11
      = Teuchos::rcp( new Thyra::DefaultIdentityLinearOp< double >( matrix_block[1][1]->range()));

    RCP< const Thyra::LinearOpBase< double > > mAinvB1, mB2Ainv;

    mAinvB1 = Thyra::scale( -1.0, Thyra::multiply( sub_op[0], matrix_block[0][1] ) );
    mB2Ainv = Thyra::scale( -1.0, Thyra::multiply( matrix_block[1][0], sub_op[0] ) );

    RCP< const Thyra::LinearOpBase< double > > Linv, Dinv, Uinv, Eye;

    Linv = Thyra::block2x2( eye_00, zero_01, mB2Ainv, eye_11 );
    Dinv = Thyra::block2x2( sub_op[0], zero_01, zero_10, sub_op[1] );
    Uinv = Thyra::block2x2( eye_00, mAinvB1, zero_10, eye_11 );

    //preconditioner = Thyra::multiply(Uinv,Dinv);
    //preconditioner = Thyra::multiply(Dinv,Linv);
    preconditioner = Thyra::multiply( Uinv, Dinv, Linv );
  }

  GEOSX_MARK_END( PRECONDITIONER );
  double setupTime = clock.stop();

  // define solver strategy for blocked system. this is
  // similar but slightly different from the sub operator
  // construction, since now we have a user defined preconditioner

  {
    RCP< Teuchos::ParameterList > list = rcp( new Teuchos::ParameterList( "list" ));

    list->set( "Linear Solver Type", "AztecOO" );
    list->set( "Preconditioner Type", "None" ); // will use user-defined P
    list->sublist( "Linear Solver Types" ).sublist( "AztecOO" ).sublist( "Forward Solve" ).set( "Max Iterations",
                                                                                                m_linearSolverParameters.krylov.maxIterations );
    list->sublist( "Linear Solver Types" ).sublist( "AztecOO" ).sublist( "Forward Solve" ).set( "Tolerance", m_linearSolverParameters.krylov.relTolerance );

    if( use_bicgstab )
      list->sublist( "Linear Solver Types" ).sublist( "AztecOO" ).sublist( "Forward Solve" ).sublist( "AztecOO Settings" ).set( "Aztec Solver", "BiCGStab" );
    else
      list->sublist( "Linear Solver Types" ).sublist( "AztecOO" ).sublist( "Forward Solve" ).sublist( "AztecOO Settings" ).set( "Aztec Solver", "GMRES" );

    if( m_linearSolverParameters.logLevel > 1 )
      list->sublist( "Linear Solver Types" ).sublist( "AztecOO" ).sublist( "Forward Solve" ).sublist( "AztecOO Settings" ).set( "Output Frequency", 1 );

    Stratimikos::DefaultLinearSolverBuilder builder;
    builder.setParameterList( list );

    RCP< const Thyra::LinearOpWithSolveFactoryBase< double > > strategy = createLinearSolveStrategy( builder );
    RCP< Thyra::LinearOpWithSolveBase< double > > solver = strategy->createOp();

    Thyra::initializePreconditionedOp< double >( *strategy,
                                                 matrix,
                                                 Thyra::rightPrec< double >( preconditioner ),
                                                 solver.ptr());

    clock.start( true );
    GEOSX_MARK_BEGIN( SOLVER );

    // !!!! Actual Solve !!!!
    Thyra::SolveStatus< double > status = solver->solve( Thyra::NOTRANS, *rhs, lhs.ptr());

    GEOSX_MARK_END( SOLVER );
    double solveTime = clock.stop();

    /* TODO: replace with SolverBase status output */

    integer numKrylovIter = status.extraParameters->get< int >( "Iteration Count" );
    if( getLogLevel()>=2 )
    {
      GEOSX_LOG_RANK_0( "\t\tLinear Solver | Iter = " << numKrylovIter <<
                        " | TargetReduction " << m_linearSolverParameters.krylov.relTolerance <<
                        " | AuxTime " << auxTime <<
                        " | SetupTime " << setupTime <<
                        " | SolveTime " << solveTime );
    }


    p_solution[1]->Scale( m_pressureScaling );
    p_rhs[1]->Scale( 1/(m_pressureScaling*m_densityScaling));
  }

  delete schurApproxPP;

  //TODO: remove all this once everything is working
  if( getLogLevel() == 2 )
  {
    /*
       ParallelVector permutedSol;
       ParallelVector const & solution = m_solidSolver->getSystemSolution();
       permutedSol.createWithLocalSize(m_solidSolver->getSystemMatrix().numLocalRows(), MPI_COMM_GEOSX);
       m_permutationMatrix0.multiply(solution, permutedSol);
       permutedSol.close();
     */

    /*
       GEOSX_LOG_RANK_0("***********************************************************");
       GEOSX_LOG_RANK_0("solution0");
       GEOSX_LOG_RANK_0("***********************************************************");
       solution.print(std::cout);
       std::cout<<std::endl;
       MPI_Barrier(MPI_COMM_GEOSX);

       GEOSX_LOG_RANK_0("***********************************************************");
       GEOSX_LOG_RANK_0("solution0");
       GEOSX_LOG_RANK_0("***********************************************************");
       permutedSol.print(std::cout);
     */
    GEOSX_LOG_RANK_0( "***********************************************************" );
    GEOSX_LOG_RANK_0( "solution1" );
    GEOSX_LOG_RANK_0( "***********************************************************" );
    p_solution[1]->Print( std::cout );

  }
#endif
}

real64
HydrofractureSolver::ScalingForSystemSolution( DomainPartition const * const domain,
                                               DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                               ParallelVector const & GEOSX_UNUSED_PARAM( solution ) )
{
  return m_solidSolver->ScalingForSystemSolution( domain,
                                                  m_solidSolver->getDofManager(),
                                                  m_solidSolver->getSystemSolution() );
}

void HydrofractureSolver::SetNextDt( real64 const & currentDt,
                                     real64 & nextDt )
{

  if( m_numResolves[0] == 0 && m_numResolves[1] == 0 )
  {
    this->SetNextDtBasedOnNewtonIter( currentDt, nextDt );
  }
  else
  {
    SolverBase * const surfaceGenerator =  this->getParent()->GetGroup< SolverBase >( "SurfaceGen" );
    nextDt = surfaceGenerator->GetTimestepRequest() < 1e99 ? surfaceGenerator->GetTimestepRequest() : currentDt;
  }
  GEOSX_LOG_LEVEL_RANK_0( 3, this->getName() << ": nextDt request is "  << nextDt );
}

void HydrofractureSolver::initializeNewFaceElements( DomainPartition const & )
{
//  m_flowSolver->
}

REGISTER_CATALOG_ENTRY( SolverBase, HydrofractureSolver, std::string const &, Group * const )
} /* namespace geosx */
