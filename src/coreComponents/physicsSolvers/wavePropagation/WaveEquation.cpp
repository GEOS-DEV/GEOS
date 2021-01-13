/*
 * WaveEquation.cpp
 *
 *  Created on: Jan 12, 2021
 *      Author: settgast
 */

#include "WaveEquation.hpp"

#include "dataRepository/KeyNames.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "constitutive/solid/SolidBase.hpp"

namespace geosx
{

using namespace dataRepository;

WaveEquation::WaveEquation( const std::string & name,
                            Group * const parent ):
  SolverBase( name,
              parent )
{}

WaveEquation::~WaveEquation()
{
  // TODO Auto-generated destructor stub
}


void WaveEquation::RegisterDataOnMesh( Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    NodeManager * const nodes = mesh.second->group_cast< MeshBody * >()->getMeshLevel( 0 )->getNodeManager();

    nodes->registerExtrinsicData< extrinsicMeshData::Pressure_nm1,
                                  extrinsicMeshData::Pressure_n,
                                  extrinsicMeshData::Pressure_np1 >( this->getName() );

  }
}


void WaveEquation::InitializePreSubGroups( Group * const rootGroup )
{
  SolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup< DomainPartition >( keys::domain );

  // Validate solid models in target regions
  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel & meshLevel = *Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );
    ValidateModelMapping< constitutive::SolidBase >( *meshLevel.getElemManager(), m_solidMaterialNames );
  }

  NumericalMethodsManager const & numericalMethodManager = domain->getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const * const
  feDiscretization = feDiscretizationManager.GetGroup< FiniteElementDiscretization >( m_discretizationName );
  GEOSX_ERROR_IF( feDiscretization == nullptr, getName() << ": FE discretization not found: " << m_discretizationName );
}


void WaveEquation::InitializePostInitialConditions_PreSubGroups( Group * const problemManager )
{
  DomainPartition * domain = problemManager->GetGroup< DomainPartition >( keys::domain );
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  NodeManager & nodes = *mesh.getNodeManager();

  ElementRegionManager & elementRegionManager = *mesh.getElemManager();

  arrayView1d< real64 > & mass = nodes.getReference< array1d< real64 > >( keys::Mass );

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > >
  rho = elementRegionManager.ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 const > >( "density",
                                                                                                              targetRegionNames(),
                                                                                                              solidMaterialNames() );

  forTargetRegionsComplete( mesh, [&]( localIndex const,
                                       localIndex const er,
                                       ElementRegionBase & elemRegion )
  {
    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const esr, CellElementSubRegion & elementSubRegion )
    {

      arrayView2d< real64 const > const & detJ = elementSubRegion.detJ();
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
        for( localIndex k=0; k < elemsToNodes.size( 0 ); ++k )
        {
          real64 elemMass = 0;
          for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
          {
            elemMass += rho[er][esr][k][q] * detJ[k][q];
            FE_TYPE::calcN( q, N );

            for( localIndex a=0; a< numNodesPerElem; ++a )
            {
              mass[elemsToNodes[k][a]] += rho[er][esr][k][q] * detJ[k][q] * N[a];
            }
          }
        }
      } );
    } );
  } );

}


real64 WaveEquation::SolverStep( real64 const & time_n,
                   real64 const & dt,
                   integer const cycleNumber,
                   DomainPartition & domain )
{
  return ExplicitStep( time_n, dt, cycleNumber, domain );
}

real64 WaveEquation::ExplicitStep( real64 const & time_n,
                                   real64 const & dt,
                                   integer const cycleNumber,
                                   DomainPartition & domain )
{

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

} /* namespace geosx */
