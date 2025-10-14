

#include "LaplaceExample.hpp"
#include "mesh/DomainPartition.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "finiteElement/elementFormulations/FiniteElementBase.hpp"


namespace geos
{
using namespace dataRepository;

LaplaceExample::LaplaceExample( const string & name,
                                Group * const parent ):
  PhysicsSolverBase( name, parent ),
  m_fieldName("pressure")
{
    this->registerWrapper( viewKeyStruct::fieldVarName(), &m_fieldName ).
      setInputFlag( InputFlags::REQUIRED ).
      setDescription( "Name of field variable" );
}


void LaplaceExample::registerDataOnMesh( Group & meshBodies )
{
  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    NodeManager & nodes = meshBody.getBaseDiscretization().getNodeManager();

    nodes.registerWrapper< real64_array >( m_fieldName ).
      setApplyDefaultValue( 0.0 ).
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setDescription( "Primary field variable" );
  } );
}

real64 LaplaceExample::solverStep( real64 const & time_n,
                                  real64 const & dt,
                                  const int cycleNumber,
                                  DomainPartition & domain )
{
  return this->linearImplicitStep( time_n, dt, cycleNumber, domain );
}

void LaplaceExample::setupDofs( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                DofManager & dofManager ) const
{
  dofManager.addField( m_fieldName,
                       FieldLocation::Node,
                       1,
                       getMeshTargets() );

  dofManager.addCoupling( m_fieldName,
                          m_fieldName,
                          DofManager::Connector::Elem );
}

void LaplaceExample::implicitStepSetup( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                       real64 const & GEOS_UNUSED_PARAM( dt ),
                                       DomainPartition & domain )
{
  Timestamp const meshModificationTimestamp = getMeshModificationTimestamp( domain );

  // Only build the sparsity pattern if the mesh has changed
  if( meshModificationTimestamp > getSystemSetupTimestamp() )
  {
    setupSystem( domain, m_dofManager, m_localMatrix, m_rhs, m_solution );
    setSystemSetupTimestamp( meshModificationTimestamp );
  }
}

void LaplaceExample::assembleSystem( real64 const GEOS_UNUSED_PARAM( time_n ),
                                 real64 const,
                                 DomainPartition & domain,
                                 DofManager const & dofManager,
                                 CRSMatrixView< real64, globalIndex const > const & matrix,
                                 arrayView1d< real64 > const & rhs )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & targetRegions )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    string const dofKey = dofManager.getKey( m_fieldName );
    arrayView1d< globalIndex const > const & dofIndex =  nodeManager.getReference< array1d< globalIndex > >( dofKey );
    const localIndex dofRankOffset = dofManager.rankOffset();

    string const finiteElementName = this->getDiscretizationName();


    ElementRegionManager & elementRegionManager = mesh.getElemManager();

    // Loop over all sub-regions in regions of type SUBREGION_TYPE, that are listed in the targetRegions array.
    elementRegionManager.forElementSubRegions< CellElementSubRegion >( targetRegions,
                                                                [&] ( localIndex const, 
                                                                      auto & elementSubRegion )
    {
      localIndex const numElems = elementSubRegion.size();

      finiteElement::FiniteElementBase const & subRegionFE = elementSubRegion.template getReference< finiteElement::FiniteElementBase >( finiteElementName );
      finiteElement::FiniteElementDispatchHandler< BASE_FE_TYPES >::dispatch3D( subRegionFE,
                                                                                    [&] ( auto const finiteElement )
      {
        using FE_TYPE = decltype(finiteElement);
        static constexpr localIndex numSupportPointsPerElem = FE_TYPE::numNodes;
        static constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;



        // get data views off the mesh objects
        auto const m_elemsToNodes = elementSubRegion.nodeList().toViewConst();
        arrayView1d<int const> const elemGhostRank = elementSubRegion.ghostRank();
        arrayView1d<real64 const> const primaryField = nodeManager.template getReference< array1d< real64 > >( m_fieldName );
        arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition();

        // this is the kernel launch
        for( localIndex k=0; k<numElems; ++k )
        {
          // gather data from mesh and store in local stack arrays
          real64 xLocal[numSupportPointsPerElem][3];
          real64 localPrimaryFieldVar[numSupportPointsPerElem];
          globalIndex localRowDofIndex[numSupportPointsPerElem];
          globalIndex localColDofIndex[numSupportPointsPerElem];
          for( localIndex a=0; a<numSupportPointsPerElem; ++a )
          {
            localIndex const localNodeIndex = m_elemsToNodes( k, a );
            xLocal[a][0] = X(localNodeIndex,0);
            xLocal[a][1] = X(localNodeIndex,1);
            xLocal[a][2] = X(localNodeIndex,2);
            localPrimaryFieldVar[ a ] = primaryField[ localNodeIndex ];
            localRowDofIndex[a] = dofIndex[localNodeIndex];
            localColDofIndex[a] = dofIndex[localNodeIndex];
          }

          real64 localJacobian[numSupportPointsPerElem][numSupportPointsPerElem] = {{0}};
          real64 localResidual[numSupportPointsPerElem] = {0};

          // quadrature point loop
          for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
          {
            real64 dNdX[ numSupportPointsPerElem ][ 3 ];
            real64 const detJ = FE_TYPE::calcGradN( q, xLocal, dNdX );
            for( localIndex a = 0; a < numSupportPointsPerElem; ++a )
            {
              for( localIndex b = 0; b < numSupportPointsPerElem; ++b )
              {
                localJacobian[ a ][ b ] += LvArray::tensorOps::AiBi< 3 >( dNdX[a], dNdX[b] ) * detJ;
              }
            }
          } // end quadrature point loop

          // form element residual from the fully formed element Jacobian dotted with
          // the primary field and map the element local Jacobian/Residual to the
          // global matrix/vector.
          for( localIndex a = 0; a < numSupportPointsPerElem; ++a )
          {
            for( localIndex b = 0; b < numSupportPointsPerElem; ++b )
            {
              localResidual[ a ] += localJacobian[ a ][ b ] * localPrimaryFieldVar[ b ];
            }
          }

          // scatter system to global matrix and rhs
          static constexpr localIndex numRows = numSupportPointsPerElem;
          static constexpr localIndex numCols = numSupportPointsPerElem;
          for( int a = 0; a < numRows; ++a )
          {
            localIndex const dof = LvArray::integerConversion< localIndex >( localRowDofIndex[ a ] - dofRankOffset );
            if( dof < 0 || dof >= matrix.numRows() ) continue;
            matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                                    localColDofIndex,
                                                                                    localJacobian[ a ],
                                                                                    numCols );

            RAJA::atomicAdd< parallelDeviceAtomic >( &rhs[ dof ], localResidual[ a ] );
          }


        } // end of kernel


      } );
    } );
  } );
}


void LaplaceExample::applyBoundaryConditions( real64 const time_n,
                                             real64 const dt,
                                             DomainPartition & domain,
                                             DofManager const & dofManager,
                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                             arrayView1d< real64 > const & localRhs )
{
  real64 const time = time_n + dt;
  FieldSpecificationManager const & fsManager = FieldSpecificationManager::getInstance();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & )
  {

    fsManager.apply< NodeManager >( time,
                                    mesh,
                                    m_fieldName,
                                    [&]( FieldSpecificationBase const & bc,
                                         string const &,
                                         SortedArrayView< localIndex const > const & targetSet,
                                         NodeManager & targetGroup,
                                         string const & GEOS_UNUSED_PARAM( fieldName ) )
    {
      bc.applyBoundaryConditionToSystem< FieldSpecificationEqual,
                                         parallelDevicePolicy< > >( targetSet,
                                                                    time,
                                                                    targetGroup,
                                                                    m_fieldName,
                                                                    dofManager.getKey( m_fieldName ),
                                                                    dofManager.rankOffset(),
                                                                    localMatrix,
                                                                    localRhs );
    } );
  } );
}


void LaplaceExample::applySystemSolution( DofManager const & dofManager,
                                         arrayView1d< real64 const > const & localSolution,
                                         real64 const scalingFactor,
                                         real64 const dt,
                                         DomainPartition & domain )
{
  GEOS_UNUSED_VAR( dt );
  dofManager.addVectorToField( localSolution,
                               m_fieldName,
                               m_fieldName,
                               scalingFactor );


  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addFields( FieldLocation::Node, { m_fieldName } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );
}



REGISTER_CATALOG_ENTRY( PhysicsSolverBase, LaplaceExample, string const &, Group * const )


} // namespace geos