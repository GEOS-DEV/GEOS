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

      arrayView2d< real64 const > const & detJ = elementSubRegion.detJ();
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

      /// get the map element to faces
      arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();

      arrayView1d< real64 > const c = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumVelocity >();
      std::cout << "Size of the medium velocity array c is " << c.size() << std::endl;
      std::cout << "c[c.size()-1] = " << c[c.size()-1] << std::endl;

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
/// Returns the value of the second derivative of a Ricker at time t0 with central Fourier frequency f0
real64 AcousticWaveEquationSEM::evaluateSecondDerivativeRicker( real64 const & t0, real64 const & f0 )
{


  // derivative of the Ricker
  real64 der_pulse = 0.0;
  real64 T0 = 1.0/f0;
  if((t0 <= -0.9*T0) || (t0 >= 2.9*T0))
    return der_pulse;

  real64 pi=3.14;
  real64 tmp = f0*t0-1.0;
  real64 f0tm1_2 = (tmp*pi)*(tmp*pi);
  der_pulse = 2*pi*pi*f0*tmp*(3.0 - 2.0*f0tm1_2)*exp( -f0tm1_2 );
  return der_pulse;
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

      /// get the barycenter of the element
      //arrayView2d< real64 const > const & elemCenterLocation = elementSubRegion.getElementCenter();
      /// The source location
      //real64 sourceLocation[3]={0};
      //localIndex sourceInElemIndex=0;

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

            for( localIndex i=0; i<numNodesPerElem; ++i )
            {
              for( localIndex j=0; j<numNodesPerElem; ++j )
              {
                Rh_k[i][j] = 0.0;
                for( localIndex a=0; a < 3; ++a )
                {
                  Rh_k[i][j] +=  detJ * gradN[i][a]*gradN[j][a];
                }
              }
            }
            ///Compute local Rh_k*p_n and save in the global vector
            for( localIndex l=0; l< numNodesPerElem; ++l )
            {
              stiffnessVector[elemsToNodes[k][q]] += Rh_k[q][l]*p_n[elemsToNodes[k][l]];
            }
          }

  	  /*
          ///Try to get an element in the center of the domain
	  //if( std::abs( elemCenterLocation[k][0]- 505.0 ) <=1.0 && (std::abs( elemCenterLocation[k][1]- 505.0 ) <=1.0 && std::abs( elemCenterLocation[k][2]- 505.0 ) <=1.0)) // 1KM WITH 100 POINTS
          
          //if( std::abs( elemCenterLocation[k][0]- 1515.0 ) <=1.0 && (std::abs( elemCenterLocation[k][1]- 1515.0 ) <=1.0 && std::abs( elemCenterLocation[k][2]- 1515.0 ) <=1.0)) //3KM WITH 100 POINTS
	  if( std::abs( elemCenterLocation[k][0]- 1010.0 ) <=1.0 && (std::abs( elemCenterLocation[k][1]- 1010.0 ) <=1.0 && std::abs( elemCenterLocation[k][2]- 1010.0 ) <=1.0)) // 2KM WITH 100 POINTS
          
          {
            sourceInElemIndex = k;
            sourceLocation[0] = elemCenterLocation[k][0];
            sourceLocation[1] = elemCenterLocation[k][1];
            sourceLocation[2] = elemCenterLocation[k][2];
	    std::cout << sourceLocation[0] << "   " << sourceLocation[1] << "   " << sourceLocation[2] << std::endl;
            real64 frequency = 5.0;
            real64 fi =0.0;
            //if(time_n <=0.4)
            //{
            fi = this->evaluateRicker( time_n, frequency );

            std::cout << "Ricker at t = " << time_n << "s is fi = " << fi << std::endl;
            for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
            {
              FE_TYPE::calcN( q, N );
              real64 const detJ = finiteElement.template getGradN< FE_TYPE >( sourceInElemIndex, q, xLocal, gradN );

              rhs[elemsToNodes[k][q]] +=  fi* detJ * N[q];
              //std::cout << "For index " << elemsToNodes[k][q] << " rhs = " << rhs[elemsToNodes[k][q]] << std::endl;
            }
          }
  	  */
	  
        }
      } );
    } );
  } );

  // apply the source before marching ahead in time
  applyRickerSource( time_n, domain );
  
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

void AcousticWaveEquationSEM::applyRickerSource( real64 const time,
                                                 DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();
  FunctionManager const & functionManager = FunctionManager::instance();

  // get the nodeManager to access the rhs
  NodeManager & nodes = *domain.getMeshBody( 0 )->getMeshLevel( 0 )->getNodeManager();
  // get the rhs vector that we want to fill in this function
  arrayView1d< real64 > const rhs = nodes.getExtrinsicData< extrinsicMeshData::ForcingRHS >();
  rhs.setValues< serialPolicy >( 0.0 );
  // get the position of the nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodes.referencePosition().toViewConst();
  
  fsManager.apply( time,
                   &domain,
                   "ElementRegions",
                   string( "Ricker" ),
                   [&]( FieldSpecificationBase const * const fs,
                        string const &,
                        SortedArrayView< localIndex const > const & targetElemSet,
                        Group * const group,
                        string const & )
  {
    // get the name of the function specified in the XML file 
    string const & functionName = fs->getFunctionName();

    // get the elemsToNodes map for this subRegion
    CellElementSubRegion * subRegion = group->groupCast< CellElementSubRegion * >();
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = subRegion->nodeList();

    // at this point, "value" should contain the value of the Ricker
    // we have to multiply it by N[a] * detJ. This is done below.

    finiteElement::FiniteElementBase const &
      fe = subRegion->getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
    finiteElement::dispatch3D( fe,
                               [&]
                               ( auto const finiteElement )
    {
      using FE_TYPE = TYPEOFREF( finiteElement );

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
      
      real64 N[numNodesPerElem]{};
      real64 gradN[numNodesPerElem][3] = {{ 0.0 }};
      real64 xLocal[numNodesPerElem][3] = {{ 0.0 }};      
      
      if( functionName.empty() || functionManager.getGroupReference< FunctionBase >( functionName ).isFunctionOfTime() == 2 )
      {
	// we can uncomment this and use the data from the table
        /*
        real64 value = fs->getScale(); // this is the "scale" of the XML file (by default, scale=1)

        // if no table is specified, then the source is constant, equal to the value of "scale"
        if( !functionName.empty() )
        {
          // if a table is specified, then interpolate in the table
          FunctionBase const & function = functionManager.getGroupReference< FunctionBase >( functionName );
          value *= function.evaluate( &time );
        }
        */

	/// For now we call the Ricker function to compute value
        real64 const frequency = 5.0;
        real64 const value = evaluateRicker( time, frequency );

        // we check that the targetElemSet only contain one element
        GEOSX_ERROR_IF( targetElemSet.size() > 1, "the target set of applyRickerSource contains more than one element. For now, the box for the source must contain only one element" );

	// we check that the targetElemSet only contain at least one element
        GEOSX_ERROR_IF( targetElemSet.size() < 1, "the target set of applyRickerSource does not contain any element of the mesh. Please Specify a correct box coordinate for the source in the xml file." );
        
        // apply source here
        for( localIndex i = 0; i < targetElemSet.size(); ++i )
        {
          localIndex const k = targetElemSet[ i ];
          for( localIndex a=0; a< numNodesPerElem; ++a )
          {
            for( localIndex inode=0; inode<3; ++inode )
            {
              xLocal[a][inode] = X( elemsToNodes( k, a ), inode );
            }
          }

          for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
          {
            FE_TYPE::calcN( q, N );
            real64 const detJ = finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );
            rhs[elemsToNodes[k][q]] += value * detJ * N[q];
          }
        }
      }
      else
      {
        // we can discuss this later, but for now let's not use this option
        GEOSX_ERROR( "This option is not supported yet" );
      }
    } );
  } );
}
  
REGISTER_CATALOG_ENTRY( SolverBase, AcousticWaveEquationSEM, string const &, dataRepository::Group * const )

} /* namespace geosx */
