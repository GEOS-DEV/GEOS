/*
 * WaveEquation.cpp
 *
 *  Created on: Jan 12, 2021
 *      Author: settgast
 */

#include "AcousticWaveEquationSEM.hpp"

#include "dataRepository/KeyNames.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"

namespace geosx
{

using namespace dataRepository;

AcousticWaveEquationSEM::AcousticWaveEquationSEM( const std::string & name,
                                                  Group * const parent ):
  SolverBase( name,
              parent )
{

  registerWrapper( viewKeyStruct::sourceCoordinatesString, &m_sourceCoordinates )->
    setInputFlag( InputFlags::REQUIRED )->
    setSizedFromParent( 0 )->
    setDescription( "Coordinates (x,y,z) of the sources" );

  registerWrapper( viewKeyStruct::sourceNodeIdsString, &m_sourceNodeIds )->
    setInputFlag( InputFlags::FALSE )->
    setSizedFromParent( 0 )->
    setDescription( "Indices of the nodes (in the right order) for each source point" );

  registerWrapper( viewKeyStruct::sourceConstantsString, &m_sourceConstants )->
    setInputFlag( InputFlags::FALSE )->
    setSizedFromParent( 0 )->
    setDescription( "Constant part of the source for the nodes listed in m_sourceNodeIds" );

  registerWrapper( viewKeyStruct::sourceIsLocalString, &m_sourceIsLocal )->
    setInputFlag( InputFlags::FALSE )->
    setSizedFromParent( 0 )->
    setDescription( "Flag that indicates whether the source is local to this MPI rank" );


  registerWrapper( viewKeyStruct::timeSourceFrequencyString, &m_timeSourceFrequency )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Central frequency for the time source" );

}

AcousticWaveEquationSEM::~AcousticWaveEquationSEM()
{
  // TODO Auto-generated destructor stub
}

void AcousticWaveEquationSEM::postProcessInput()
{

  // here we check that the correct number of dimensions has been provided
  GEOSX_ERROR_IF( m_sourceCoordinates.size( 1 ) != 3,
                  "Invalid number of physical coordinates for the sources" );

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );
  localIndex const numNodesPerElem = 8;
  m_sourceNodeIds.resizeDimension< 0, 1 >( numSourcesGlobal, numNodesPerElem );
  m_sourceConstants.resizeDimension< 0, 1 >( numSourcesGlobal, numNodesPerElem );
  m_sourceIsLocal.resizeDimension< 0 >( numSourcesGlobal );

}

void AcousticWaveEquationSEM::registerDataOnMesh( Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->getSubGroups() )
  {

    MeshLevel & meshLevel = *(mesh.second->groupCast< MeshBody * >()->getMeshLevel( 0 ));

    NodeManager & nodes = *(meshLevel.getNodeManager());

    nodes.registerExtrinsicData< extrinsicMeshData::Pressure_nm1,
                                 extrinsicMeshData::Pressure_n,
                                 extrinsicMeshData::Pressure_np1,
                                 extrinsicMeshData::ForcingRHS,
                                 extrinsicMeshData::MassVector,
                                 extrinsicMeshData::DampingVector,
                                 extrinsicMeshData::StiffnessVector >( this->getName() );


    ElementRegionManager & elemManager = *(meshLevel.getElemManager());

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerExtrinsicData< extrinsicMeshData::MediumVelocity >( this->getName() );
    } );


  }
}

void AcousticWaveEquationSEM::precomputeSourceTerm( MeshLevel & mesh )
{
  NodeManager & nodeManager = *(mesh.getNodeManager());

  // get the position of the nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

  // get the source information
  arrayView2d< real64 const > const sourceCoordinates = m_sourceCoordinates.toViewConst();
  arrayView2d< localIndex > const sourceNodeIds = m_sourceNodeIds.toView();
  arrayView2d< real64 > const sourceConstants = m_sourceConstants.toView();
  arrayView1d< localIndex > const sourceIsLocal = m_sourceIsLocal.toView();
  sourceNodeIds.setValues< serialPolicy >( -1 );
  sourceConstants.setValues< serialPolicy >( -1 );
  sourceIsLocal.setValues< serialPolicy >( 0 );

  forTargetRegionsComplete( mesh, [&]( localIndex const,
                                       localIndex const,
                                       ElementRegionBase & elemRegion )
  {
    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const,
                                                                       CellElementSubRegion & elementSubRegion )
    {

      // get the face/node information
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::dispatch3D( fe,
                                 [&]
                                   ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        localIndex const numFacesPerElem = elementSubRegion.numFacesPerElement();
        localIndex constexpr numNodesPerElem = FE_TYPE::numNodes;
        array1d< array1d< localIndex > > faceNodes( numFacesPerElem );

        // loop over all the elements in the subRegion
        // this will potentially become a RAJA loop
        for( localIndex k = 0; k < elementSubRegion.size(); ++k )
        {

          // collect the faces for this element
          for( localIndex kf = 0; kf < numFacesPerElem; ++kf )
          {
            elementSubRegion.getFaceNodes( k, kf, faceNodes[kf] );
          }

          // loop over all the sources that haven't been found yet
          // If we don't care about having multiple sources, we can remove this loop
          for( localIndex isrc = 0; isrc < sourceCoordinates.size( 0 ); ++isrc )
          {
            // if the source has not been found yet
            if( sourceIsLocal[isrc] == 0 )
            {
              // get the coordinates of the source
              real64 const coords[3] = { sourceCoordinates[isrc][0],
                                         sourceCoordinates[isrc][1],
                                         sourceCoordinates[isrc][2] };

              // if the point is in the element, we can compute the constant part of the source term
              // The search can be optimized a lot
              if( computationalGeometry::IsPointInsidePolyhedron( X, faceNodes, coords ) )
              {
                sourceIsLocal[isrc] = 1;
                std::cout << "I found the source in element " << k << " at location ("
                          << coords[0] << ", " << coords[1] << ", " << coords[2] << ")" << std::endl;

                /// Get all the node of element k containing the source point
                real64 xLocal[numNodesPerElem][3];
                for( localIndex a=0; a< numNodesPerElem; ++a )
                {
                  std::cout << " For node " << a;
                  for( localIndex i=0; i<3; ++i )
                  {
                    xLocal[a][i] = X( elemsToNodes( k, a ), i );
                    std::cout << " x_"<< i << " = " << xLocal[a][i];
                  }
                  std::cout << " " << std::endl;
                }
                std::cout << "xLocal Ok "<< std::endl;

                /// Transformation of the source coordinate to unit reference element coordinate
                /// coordsOnRefElem = invJ*(coords-coordsNode_0)
                real64 coordsOnRefElem[3]; // 3D coord of the source in Ref element
                //real64 valsOnRefElem = 0.0; // Value of the basis function evaluated at coordsOnRefElem
                localIndex q=0; // index of node 0 in the element

                ///
                //real64 N[ numNodesPerElem ];
                real64 gradN[ numNodesPerElem ][ 3 ];
                //FE_TYPE::calcN( q, N );
                real64 const detJ = finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );

                ///Compute invJ = DF^{-1}
                real64 invJ[3][3]={{0}};
                FE_TYPE::invJacobianTransformation( q, xLocal, invJ );
                std::cout << "invJ Ok "<< std::endl;

                /// compute (coords - coordsNode_0)
                real64 coordsRef[3]={0};
                for( localIndex i=0; i<3; ++i )
                {
                  coordsRef[i] = coords[i] - xLocal[q][i];
                }

                /// Compute coordsOnRefElem = invJ*coordsRef
                for( localIndex i=0; i<3; ++i )
                {
                  // Init at (-1,-1,-1) as the origin of the referential elem
                  coordsOnRefElem[i] =-1.0;
                  for( localIndex j=0; j<3; ++j )
                  {
                    coordsOnRefElem[i] += invJ[i][j]*coordsRef[j];
                  }
                }

                /* std::cout << " V1 coordsOnRefElem " << coordsOnRefElem[0] << " " << coordsOnRefElem[1] << " " << coordsOnRefElem[2] <<
                   std::endl;

                   if(numNodesPerElem ==8)
                   {
                    /// V2 to compute coords On ref Elem
                    coordsOnRefElem[0] = -1.0 + 2.0*(coords[0]-xLocal[0][0])/(xLocal[1][0] - xLocal[0][0]);
                    coordsOnRefElem[1] = -1.0 + 2.0*(coords[1]-xLocal[0][1])/(xLocal[2][1] - xLocal[0][1]);
                    coordsOnRefElem[2] = -1.0 + 2.0*(coords[2]-xLocal[0][2])/(xLocal[4][2] - xLocal[0][2]);

                    std::cout << " V2 coordsOnRefElem " << coordsOnRefElem[0] << " " << coordsOnRefElem[1] << " " << coordsOnRefElem[2] <<
                       std::endl;
                   }
                 */

                /// Evaluate basis functions at coord source on unit ref element
                real64 Ntest[8];
                finiteElement::LagrangeBasis1::TensorProduct3D::value( coordsOnRefElem, Ntest );

                std::cout << "Ntest Ok "<< std::endl;

                // save all the node indices and constant part of source term here
                for( localIndex a=0; a< numNodesPerElem; ++a )
                {

                  sourceNodeIds[isrc][a] = elemsToNodes[k][a];
                  sourceConstants[isrc][a] = detJ*Ntest[a]; // precompute detJ x N[position of source] here

                  std::cout << "For source #" << isrc << " I save node #" << sourceNodeIds[isrc][a] << " and constant value = " << sourceConstants[isrc][a] << std::endl;

                }
              }
            }
          }
        }
      } );
    } );
  } );
  ///GEOSX_ERROR_IF( true, "Stop test " );
}

void AcousticWaveEquationSEM::addSourceToRightHandSide( real64 const & time, arrayView1d< real64 > const rhs )
{
  // get the precomputed source information
  arrayView2d< localIndex const > const sourceNodeIds = m_sourceNodeIds.toViewConst();
  arrayView2d< real64 const > const sourceConstants   = m_sourceConstants.toViewConst();
  arrayView1d< localIndex const > const sourceIsLocal = m_sourceIsLocal.toViewConst();

  real64 const fi = evaluateRicker( time, this->m_timeSourceFrequency );

  // loop over all the sources
  for( localIndex isrc = 0; isrc < sourceConstants.size( 0 ); ++isrc )
  {
    if( sourceIsLocal[isrc] == 1 ) // check if the source is on this MPI rank
    {
      // for all the nodes
      for( localIndex inode = 0; inode < sourceConstants.size( 1 ); ++inode )
      {
        // multiply the precomputed part by the ricker
        rhs[sourceNodeIds[isrc][inode]] += sourceConstants[isrc][inode] * fi;
      }
    }
  }
}

void AcousticWaveEquationSEM::initializePreSubGroups( Group * const rootGroup )
{
  SolverBase::initializePreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->getGroup< DomainPartition >( keys::domain );

  NumericalMethodsManager const & numericalMethodManager = domain->getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const * const
  feDiscretization = feDiscretizationManager.getGroup< FiniteElementDiscretization >( m_discretizationName );
  GEOSX_ERROR_IF( feDiscretization == nullptr, getName() << ": FE discretization not found: " << m_discretizationName );
}


void AcousticWaveEquationSEM::initializePostInitialConditionsPreSubGroups( Group * const problemManager )
{
  DomainPartition * domain = problemManager->getGroup< DomainPartition >( keys::domain );
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  precomputeSourceTerm( mesh );

  NodeManager & nodeManager = *mesh.getNodeManager();
  FaceManager & faceManager = *mesh.getFaceManager();

  /// get the array of indicators: 1 if the node is on the boundary; 0 otherwise
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

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

      /// get the map element to faces
      arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();

      arrayView1d< real64 > const c = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumVelocity >();
      
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
        real64 gradN[ numNodesPerElem ][ 3 ];
        for( localIndex k=0; k < elemsToNodes.size( 0 ); ++k )
        {
          real64 const invC2 = 1.0 / ( c[k] * c[k] );
          real64 xLocal[numNodesPerElem][3];
          for( localIndex a=0; a< numNodesPerElem; ++a )
          {
            for( localIndex i=0; i<3; ++i )
            {
              xLocal[a][i] = X( elemsToNodes( k, a ), i );
            }
          }

          for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
          {
            FE_TYPE::calcN( q, N );
            real64 const detJ = finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );

            for( localIndex a=0; a< numNodesPerElem; ++a )
            {
              /// update mass matrix
              mass[elemsToNodes[k][a]] +=  invC2 * detJ * N[a];
            }
          }
        }

        real64 sumMass = 0.0;
        for( localIndex a=0; a<nodeManager.size(); ++a )
        {
          sumMass +=mass[a];
        }

        /* Unit test
           // Test for mass matrix sumTerm*c2 should be volume of the domaine
           // assuming MediumVelocity c = 1500
           sumMass *=1500*1500;
           std::cout << "Sum mass terms time C2 = " << sumMass << std::endl;

           GEOSX_ERROR_IF( true, " Stop test Mass Ok" );
         */

        /// update damping matrix
        for( localIndex k=0; k < elemsToFaces.size( 0 ); ++k )
        {
          real64 const alpha = 1.0/c[k];

          for( localIndex kfe=0; kfe< numFacesPerElem; ++kfe )
          {
            /// Face on the domain boundary
            if( facesDomainBoundaryIndicator[elemsToFaces[k][kfe]]==1 )
            {
              real64 xLocal[numNodesPerElem][3];
              for( localIndex a=0; a< numNodesPerElem; ++a )
              {
                for( localIndex i=0; i<3; ++i )
                {
                  xLocal[a][i] = X( elemsToNodes( k, a ), i );
                }
              }
              for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
              {
                FE_TYPE::calcN( q, N );
                real64 const detJ = finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );

                ///Compute invJ = DF^{-1}
                real64 invJ[3][3]={{0}};
                FE_TYPE::invJacobianTransformation( q, xLocal, invJ );

                for( localIndex a=0; a < numNodesPerFace; ++a )
                {
                  /// compute ds=||detJ*invJ*normalFace_{kfe}||
                  real64 tmp[3]={0};
                  real64 ds = 0.0;
                  for( localIndex i=0; i<3; ++i )
                  {
                    for( localIndex j = 0; j < 3; ++j )
                    {
                      tmp[i] += invJ[j][i]*faceNormal[elemsToFaces[k][kfe]][j];
                    }
                    ds +=tmp[i]*tmp[i];
                  }
                  ds = std::sqrt( ds );

                  localIndex numFaceGl = elemsToFaces[k][kfe];
                  localIndex numNodeGl = facesToNodes[numFaceGl][a];
                  damping[numNodeGl] += alpha*detJ*ds*N[a];
                }
              }
            }
          }
        }

      } );
    } );
  } );

}


real64 AcousticWaveEquationSEM::solverStep( real64 const & time_n,
                                            real64 const & dt,
                                            integer const cycleNumber,
                                            DomainPartition & domain )
{
  return explicitStep( time_n, dt, cycleNumber, domain );
}


/// Returns the value of a Ricker at time t0 with central Fourier frequency f0
real64 AcousticWaveEquationSEM::evaluateRicker( real64 const & t0, real64 const & f0 )
{
  // Center time
  real64 T0 = 1.0/f0;
  real64 pulse = 0.0;
  if((t0 <= -0.9*T0) || (t0 >= 2.9*T0))
    return pulse;

  real64 pi = 3.14;
  real64 tmp = f0*t0-1.0;
  real64 f0tm1_2 = 2*(tmp*pi)*(tmp*pi);
  real64 gaussian_term = exp( -f0tm1_2 );
  pulse = -(t0-1)*gaussian_term;

  return pulse;
}

real64 AcousticWaveEquationSEM::explicitStep( real64 const & time_n,
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

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodes.referencePosition().toViewConst();

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

        /* Unit test
           // For unit test, intialize p_n with x-coordinate of mesh vertices {(x,y,z), x,y,z \in R}
           for( localIndex a=0; a<nodes.size(); ++a )
           {
            p_n[a] =  X(a,0); //1.0; // X(a,0);
           }
         */

        for( localIndex k=0; k < elemsToNodes.size( 0 ); ++k )
        {
          real64 xLocal[numNodesPerElem][3];
          /// Local stiffness matrix for the element k
          real64 Rh_k[numNodesPerElem][numNodesPerElem] = {{0}};

          for( localIndex a=0; a< numNodesPerElem; ++a )
          {
            for( localIndex i=0; i<3; ++i )
            {
              xLocal[a][i] = X( elemsToNodes( k, a ), i );
            }
          }


          for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
          {
            ///Calculate the basis function N at the node q
            FE_TYPE::calcN( q, N );
            ///Compute gradN = invJ*\hat{\nabla}N at the node q and return the determinant of the transformation matrix J
            real64 const detJ = finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );
            //Rh_k = {{0}};

            for( localIndex i=0; i<numNodesPerElem; ++i )
            {
              for( localIndex j=0; j<numNodesPerElem; ++j )
              {
                Rh_k[i][j] = 0.0;
                for( localIndex a=0; a < 3; ++a )
                {
                  Rh_k[i][j] +=  detJ * gradN[i][a]*gradN[j][a];
                }

                ///Compute local Rh_k*p_n and save in the global vector
                stiffnessVector[elemsToNodes[k][i]] += Rh_k[i][j]*p_n[elemsToNodes[k][j]];
              }

            }

          }

        }
        /* Unit test
           // compute <\nabla p, \nabla p> = <p, Rh p>
           real64 prodScalar = 0.0;
           for( localIndex a=0; a<nodes.size(); ++a )
           {
            prodScalar +=stiffnessVector[a]*p_n[a];
           }
           std::cout << " Scalar product p_n and stiffnessVector = " << prodScalar << std::endl;
           GEOSX_ERROR_IF( true, " Stop test Stiffness matrix Ok" );
         */
      } );
    } );
  } );

  /// Add source term to rhs
  addSourceToRightHandSide( time_n, rhs );

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
  std::cout << "Pressure[505050] = " << p_n[505050] << std::endl;

  return dt;
}

REGISTER_CATALOG_ENTRY( SolverBase, AcousticWaveEquationSEM, string const &, dataRepository::Group * const )

} /* namespace geosx */
