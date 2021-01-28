/*
 * WaveEquation.cpp
 *
 *  Created on: Jan 12, 2021
 *      Author: settgast
 */

#include "AcousticWaveEquationSEM.hpp"

#include "dataRepository/KeyNames.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "constitutive/solid/SolidBase.hpp"

namespace geosx
{

using namespace dataRepository;

AcousticWaveEquationSEM::AcousticWaveEquationSEM( const std::string & name,
                            Group * const parent ):
  SolverBase( name,
              parent )
{}

AcousticWaveEquationSEM::~AcousticWaveEquationSEM()
{
  // TODO Auto-generated destructor stub
}


void AcousticWaveEquationSEM::RegisterDataOnMesh( Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {

    MeshLevel & meshLevel = *(mesh.second->group_cast< MeshBody * >()->getMeshLevel( 0 ));

    NodeManager & nodes = *(meshLevel.getNodeManager());

    nodes.registerExtrinsicData< extrinsicMeshData::Pressure_nm1,
				 extrinsicMeshData::Pressure_n,
				 extrinsicMeshData::Pressure_np1,
				 extrinsicMeshData::ForcingRHS,
				 extrinsicMeshData::MassVector,
				 extrinsicMeshData::DampingVector,
				 extrinsicMeshData::StiffnessVector>( this->getName() );


    ElementRegionManager & elemManager = *(meshLevel.getElemManager());

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerExtrinsicData< extrinsicMeshData::MediumVelocity >( this->getName() );
    } );


  }
}


void AcousticWaveEquationSEM::InitializePreSubGroups( Group * const rootGroup )
{
  SolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup< DomainPartition >( keys::domain );
  /* No constituve model
  // Validate solid models in target regions
  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel & meshLevel = *Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );
    ValidateModelMapping< constitutive::SolidBase >( *meshLevel.getElemManager(), m_solidMaterialNames );
  }
  */
  NumericalMethodsManager const & numericalMethodManager = domain->getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const * const
  feDiscretization = feDiscretizationManager.GetGroup< FiniteElementDiscretization >( m_discretizationName );
  GEOSX_ERROR_IF( feDiscretization == nullptr, getName() << ": FE discretization not found: " << m_discretizationName );
}


void AcousticWaveEquationSEM::InitializePostInitialConditions_PreSubGroups( Group * const problemManager )
{
  DomainPartition * domain = problemManager->GetGroup< DomainPartition >( keys::domain );
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  NodeManager & nodeManager = *mesh.getNodeManager();
  FaceManager & faceManager = *mesh.getFaceManager();

  /// get the array of indicators: 1 if the node is on the boundary; 0 otherwise
  arrayView1d< integer const > const nodesDomainBoundaryIndicator = nodeManager.getDomainBoundaryIndicator();
  arrayView1d< integer > const & facesDomainBoundaryIndicator = faceManager.getDomainBoundaryIndicator();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();
  
  /// Get table containing all the face normals
  arrayView2d< real64 const > const faceNormal  = faceManager.faceNormal();
  ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();
  
  arrayView1d< real64 > const mass = nodeManager.getExtrinsicData< extrinsicMeshData::MassVector >();
  /// damping matrix to be computed for each dof in the boundary of the mesh
  arrayView1d< real64 > const damping = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector >();

  forTargetRegionsComplete( mesh, [&]( localIndex const,
                                       localIndex const,
                                       ElementRegionBase & elemRegion )
  {
    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const,
                                                                       CellElementSubRegion & elementSubRegion )
    {

      arrayView2d< real64 const > const & detJ = elementSubRegion.detJ();
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

      /// get the map element to faces
      arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();
      
      arrayView1d< real64 > const c = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumVelocity >();
      std::cout << c.size() << std::endl;
      c.setValues< serialPolicy >(1500);
      
      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::dispatch3D( fe,
                                 [&]
                                 ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
        constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
	constexpr localIndex numFacesPerElem = 6; // FE_TYPE::numNodes;
	constexpr localIndex numNodesPerFace = 4; // FE_TYPE::numNodes;
	
        real64 N[numNodesPerElem];
        for( localIndex k=0; k < elemsToNodes.size( 0 ); ++k )
        {
          real64 const invC2 = 1.0 / ( c[k] * c[k] );
	  //real64 const alpha = 1.0/c[k];

          for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
          {
            FE_TYPE::calcN( q, N );

            for( localIndex a=0; a< numNodesPerElem; ++a )
            {
	      /// update mass matrix
              mass[elemsToNodes[k][a]] +=  invC2 * detJ[k][q] * N[a];
            }
          }
        }

	//std::cout << "Coucou c'est bon" << std::endl;
	/// update damping matrix
	for( localIndex k=0; k < elemsToFaces.size( 0 ); ++k )
	  {
	    real64 const alpha = 1.0/c[k];	    
	    
	    for( localIndex kfe=0; kfe< numFacesPerElem; ++kfe )
	      {
		/// Face on the domain boundary
		if(facesDomainBoundaryIndicator[elemsToFaces[k][kfe]]==1)
		  {
		    real64 xLocal[numNodesPerElem][3];
		    for( localIndex a=0; a< numNodesPerElem; ++a )
		      {
			for( localIndex i=0 ; i<3; ++i )
			  {
			    xLocal[a][i] = X( elemsToNodes(k,a), i);
			  }
		      }
		    for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
		      {
			FE_TYPE::calcN( q, N );
	
			///Compute invJ = DF^{-1}
			real64 invJ[3][3]={{0}};
			FE_TYPE::invJacobianTransformation(q,xLocal,invJ);
			
			for(localIndex a=0; a < numNodesPerFace; ++a)
			  {
			    //if(nodesDomainBoundaryIndicator[elemsToNodes[k][a]] == 1)
			    // {
			    /// compute ds=||detJ*invJ*normalFace_{kfe}||
			    real64 tmp[3]={0};
			    real64 ds = 0.0;
			    for(localIndex i=0; i<3; ++i)
			      {
				for(localIndex j = 0; j < 3; ++j)
				  {
				    tmp[i] += invJ[j][i]*faceNormal[elemsToFaces[k][kfe]][j];
				  }
				ds +=tmp[i]*tmp[i];
			      }
			    ds = std::sqrt(ds);
			    //std::cout << kfe << std::endl;
			    //std::cout << facesToNodes[elemsToFaces[k][kfe]][a] << std::endl;
			    //damping[elemsToNodes[k][facesToNodes[elemsToFaces[k][kfe]][a]]] += alpha*detJ[k][q]*ds*N[a];
			    localIndex numFaceGl = elemsToFaces[k][kfe];
			    localIndex numNodeGl = facesToNodes[numFaceGl][a];
			    damping[numNodeGl] += alpha*detJ[k][q]*ds*N[a];
			  }
		      }
		  }
	      }
	  }
	
      } );
    } );
  } );

}


real64 AcousticWaveEquationSEM::SolverStep( real64 const & time_n,
                                            real64 const & dt,
                                            integer const cycleNumber,
                                            DomainPartition & domain )
{
  return ExplicitStep( time_n, dt, cycleNumber, domain );
}


/// Returns the value of a Ricker at time t0 with central Fourier frequency f0 and center time tc
real64 AcousticWaveEquationSEM::EvaluateRicker(real64 const & t0, real64 const & tc, real64 const & f0)
{
    real64 pulse = 0;
    real64 pi = 3.14;
    real64 tmp = f0*(t0-tc);
    real64 f0tm1_2 = 2*(tmp*pi)*(tmp*pi);
    real64 gaussian_term = exp(-f0tm1_2); 
    pulse = -(t0-tc)*gaussian_term;
        
    return pulse;
  }

real64 AcousticWaveEquationSEM::ExplicitStep( real64 const & time_n,
                                              real64 const & dt,
                                              integer const cycleNumber,
                                              DomainPartition & domain )
{

  GEOSX_UNUSED_VAR( time_n, dt, cycleNumber );

  MeshLevel & mesh = *(domain.getMeshBody( 0 )->getMeshLevel( 0 ));

  NodeManager & nodes = *mesh.getNodeManager();
  
  arrayView1d< real64 > const mass = nodes.getExtrinsicData< extrinsicMeshData::MassVector >();
  arrayView1d< real64 > const damping = nodes.getExtrinsicData< extrinsicMeshData::DampingVector >();
  
  arrayView1d< real64 > const p_nm1 = nodes.getExtrinsicData< extrinsicMeshData::Pressure_nm1 >();
  arrayView1d< real64 > const p_n = nodes.getExtrinsicData< extrinsicMeshData::Pressure_n >();
  arrayView1d< real64 > const p_np1 = nodes.getExtrinsicData< extrinsicMeshData::Pressure_np1 >();
  
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodes.referencePosition().toViewConst(); //nodes.referencePosition();

  /// Vector to contain the product of the stiffness matrix R_h and the pressure p_n
  arrayView1d< real64 > const stiffnessVector = nodes.getExtrinsicData< extrinsicMeshData::StiffnessVector >();

  /// Vector to compute rhs
  arrayView1d< real64 > const rhs = nodes.getExtrinsicData< extrinsicMeshData::ForcingRHS >();
  
  forTargetRegionsComplete( mesh, [&]( localIndex const,
                                       localIndex const,
                                       ElementRegionBase & elemRegion )
  {
    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const,
                                                                       CellElementSubRegion & elementSubRegion )
    {
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

      /// get the barycenter of the elem
      arrayView2d< real64 const > const & elemCenterCoord = elementSubRegion.getElementCenter();
      
      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::dispatch3D( fe,
                                 [&]
                                 ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
        constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
	
        real64 N[numNodesPerElem];
        real64 gradN[ numNodesPerElem ][ 3 ];

        for( localIndex k=0; k < elemsToNodes.size( 0 ); ++k )
        {
          real64 xLocal[numNodesPerElem][3];
	  /// Local stiffness matrix for the element k
	  real64 Rh_k[numNodesPerElem][numNodesPerElem];
	  
          for( localIndex a=0; a< numNodesPerElem; ++a )
          {
            for( localIndex i=0 ; i<3; ++i )
            {
              xLocal[a][i] = X( elemsToNodes(k,a), i);
            }
          }


          for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
          {
	    ///Calculate the basis function N at the node q
            FE_TYPE::calcN( q, N );
	    ///Compute gradN = invJ*\hat{\nabla}N at the node q and return the determinant of the transformation matrix J
	    real64 const detJ = finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN);
	    
	    for( localIndex i=0; i<numNodesPerElem; ++i )
            {
	      for(localIndex j=0; j<numNodesPerElem; ++j )
		{
		  Rh_k[i][j] = 0.0;
		  for(localIndex a=0; a<2; ++a)
		    {
		      Rh_k[i][j] +=  detJ * gradN[i][a]*gradN[j][a];
		    }
		}
	    }
	    ///Compute local Rh_k*p_n and save in the global vector  
            for( localIndex l=0; l< numNodesPerElem; ++l )
	      {
		stiffnessVector[elemsToNodes[k][q]] += Rh_k[q][l]*p_n[elemsToNodes[k][l]] ;
	      }
          }
	  //if(elemCenterCoord[k][0]== 500.0 && ())
	  // {  
	  //  }
        }
      } );
    } );
  } );

  // rhs at 500 500 500
  
  /// Calculate your time integrators
  real64 dt2 = dt*dt;
  for( localIndex a=0; a<nodes.size(); ++a )
  {
    // pressure update here
    p_np1[a] = (1.0/(mass[a]+0.5*dt*damping[a]))*(2*mass[a]*p_n[a]-dt2*stiffnessVector[a] - (mass[a] - 0.5*dt*damping[a])*p_nm1[a] + dt2*rhs[a] );
    p_nm1[a]=p_n[a];
    p_n[a] = p_np1[a];

    stiffnessVector[a] = 0.0;
    rhs[a] = 0.0;
  }
  std::cout << "Pressure[0] = " << p_n[2] << std::endl;

  return dt;
}


//void WaveEquation::ApplyBoundaryConditions( real64 const time,
//                                            real64 const dt,
//                                            DomainPartition & domain,
//                                            DofManager const & dofManager,
//                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
//                                            arrayView1d< real64 > const & localRhs )
//{
//
//}


REGISTER_CATALOG_ENTRY( SolverBase, AcousticWaveEquationSEM, string const &, dataRepository::Group * const )

} /* namespace geosx */
