/*
 * PhysicsLoopInterface.hpp
 *
 *  Created on: Mar 24, 2020
 *      Author: settgast
 */

#ifndef GEOSX_PHYSICSSOLVERS_PHYSICSLOOPINTERFACE_HPP_
#define GEOSX_PHYSICSSOLVERS_PHYSICSLOOPINTERFACE_HPP_

#include "common/DataTypes.hpp"
#include "constitutive/solid/LinearElasticIsotropic.hpp"
#include "constitutive/solid/solidSelector.hpp"
#include "mesh/ElementRegionManager.hpp"


namespace geosx
{
namespace physicsLoopInterface
{

template< typename LAMBDA >
void
discretizationLaunchSelector( localIndex NUM_NODES_PER_ELEM,
                              localIndex NUM_QUADRATURE_POINTS,
                              LAMBDA && lambda )
{
  if( NUM_NODES_PER_ELEM==8 && NUM_QUADRATURE_POINTS==8 )
  {
    lambda( std::integral_constant< int, 8 >(), std::integral_constant< int, 8 >() );
  }
  else if( NUM_NODES_PER_ELEM==8 && NUM_QUADRATURE_POINTS==1 )
  {
    lambda( std::integral_constant< int, 8 >(), std::integral_constant< int, 1 >() );
  }
  else if( NUM_NODES_PER_ELEM==4 && NUM_QUADRATURE_POINTS==1 )
  {
    lambda( std::integral_constant< int, 4 >(), std::integral_constant< int, 1 >() );
  }
  else
  {
    GEOSX_ERROR( "Valid Branch not found." );
  }
}

class FiniteElementRegionLoopKernelBase
{
public:

  template< int NUM_NODES_PER_ELEM, int NUM_DOF_PER_NODE >
  struct StackVariables
  {
public:
    static constexpr int numNodesPerElem = NUM_NODES_PER_ELEM;
    static constexpr int ndof = NUM_DOF_PER_NODE * NUM_NODES_PER_ELEM;

    StackVariables():
      elementLocalDofIndex( ndof ),
      R( ndof ),
      dRdU( ndof, ndof )
    {
      dRdU = 0.0;
      R = 0.0;
    }

    stackArray1d< globalIndex, ndof >   elementLocalDofIndex;
    stackArray1d< real64, ndof >        R;
    stackArray2d< real64, ndof *ndof >  dRdU;

  };


  FiniteElementRegionLoopKernelBase( arrayView1d< globalIndex const > const & inputDofNumber,
                                     ParallelMatrix & inputMatrix,
                                     ParallelVector & inputRhs ):
    m_dofNumber( inputDofNumber ),
    m_matrix( inputMatrix ),
    m_rhs( inputRhs )
  {}

  arrayView1d< globalIndex const > m_dofNumber;
  ParallelMatrix & m_matrix;
  ParallelVector & m_rhs;
  arrayView2d< localIndex const, cells::NODE_MAP_USD > elemsToNodes;

  void initializeNonElementViews( MeshLevel & GEOSX_UNUSED_PARAM( mesh ) )
  {}

  template< typename ELEMENT_SUBREGION_TYPE >
  void initializeElementSubRegionViews( ELEMENT_SUBREGION_TYPE & elementSubRegion )
  {
    elemsToNodes = elementSubRegion.nodeList();
  }



  template< typename STACK_VARIABLE_TYPE >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void preKernel( localIndex const k,
                  STACK_VARIABLE_TYPE & stack ) const
  {
    for( localIndex a=0; a<STACK_VARIABLE_TYPE::numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = elemsToNodes( k, a );
      for( int i=0; i<3; ++i )
      {
        stack.elementLocalDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
      }
    }

  }

  template< typename STACK_VARIABLE_TYPE, typename CONSTITUTIVE_UPDATE >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void updateKernel( localIndex const GEOSX_UNUSED_PARAM( k ),
                     localIndex const GEOSX_UNUSED_PARAM( q ),
                     STACK_VARIABLE_TYPE & GEOSX_UNUSED_PARAM( stack ),
                     CONSTITUTIVE_UPDATE const & GEOSX_UNUSED_PARAM( constitutiveUpdate )  ) const
  {}

  template< typename STACK_VARIABLE_TYPE >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void stiffnessKernel( localIndex const GEOSX_UNUSED_PARAM( k ),
                        localIndex const GEOSX_UNUSED_PARAM( q ),
                        STACK_VARIABLE_TYPE & GEOSX_UNUSED_PARAM( stack ) ) const
  {}

  template< typename STACK_VARIABLE_TYPE, typename CONSTITUTIVE_UPDATE >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void integrationKernel( localIndex const GEOSX_UNUSED_PARAM( k ),
                          localIndex const GEOSX_UNUSED_PARAM( q ),
                          STACK_VARIABLE_TYPE & GEOSX_UNUSED_PARAM( stack ),
                          CONSTITUTIVE_UPDATE const & GEOSX_UNUSED_PARAM( constitutiveUpdate ) ) const
  {}

  template< typename STACK_VARIABLE_TYPE >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 postKernel( STACK_VARIABLE_TYPE const & stack ) const
  {
    m_matrix.insert( stack.elementLocalDofIndex.data(),
                     stack.elementLocalDofIndex.data(),
                     stack.dRdU.data(),
                     stack.ndof,
                     stack.ndof );
    return 0;
  }

};

template< typename POLICY, typename KERNEL_CLASS = FiniteElementRegionLoopKernelBase >
static
real64 FiniteElementRegionLoop( MeshLevel & mesh,
                                string_array const & targetRegions,
                                string const & solidMaterialName,
                                FiniteElementDiscretization const * const feDiscretization,
                                arrayView1d< globalIndex const > const & inputDofNumber,
                                ParallelMatrix & inputMatrix,
                                ParallelVector & inputRhs )
{
  KERNEL_CLASS kernelClass( inputDofNumber, inputMatrix, inputRhs );

  RAJA::ReduceMax< serialReduce, double > maxForce( 0 );


  kernelClass.initializeNonElementViews( mesh );

  ElementRegionManager & elementRegionManager = *(mesh.getElemManager());

  elementRegionManager.forElementSubRegions< CellElementSubRegion >( targetRegions,
                                                                     [&] ( auto & elementSubRegion )
  {
    localIndex const numElems = elementSubRegion.size();
//      typedef TYPEOFREF( elementSubRegion ) SUBREGIONTYPE;

    arrayView1d< integer const > const & elemGhostRank = elementSubRegion.ghostRank();

    localIndex const
    numQuadraturePointsPerElem = feDiscretization == nullptr ?
                                 1 :
                                 feDiscretization->m_finiteElement->n_quadrature_points();

    discretizationLaunchSelector( elementSubRegion.numNodesPerElement(),
                                  numQuadraturePointsPerElem,
                                  [&]( auto constNNPE,
                                       auto constNQPPE )
    {
      constexpr int NUM_NODES_PER_ELEM = decltype( constNNPE )::value;
      constexpr int NUM_QUADRATURE_POINTS = decltype( constNQPPE )::value;

      constitutive::SolidBase * const
      solidBase = elementSubRegion.GetConstitutiveModels()->template GetGroup< constitutive::SolidBase >( solidMaterialName );

      kernelClass.initializeElementSubRegionViews( elementSubRegion );

      constitutive::constitutiveUpdatePassThru( solidBase, [&]( auto & castedConstitutiveRelation )
      {
        using CONSTITUTIVE_TYPE = TYPEOFREF( castedConstitutiveRelation );

        typename CONSTITUTIVE_TYPE::KernelWrapper constitutiveWrapper = castedConstitutiveRelation.createKernelWrapper();

        forAll< POLICY >( numElems,
                          [=] GEOSX_HOST_DEVICE ( localIndex const k )
            {

              typename KERNEL_CLASS:: template StackVariables< NUM_NODES_PER_ELEM, 3 > stack;

              if( elemGhostRank[k] < 0 )
              {
                kernelClass.preKernel( k, stack );
                for( integer q=0; q<NUM_QUADRATURE_POINTS; ++q )
                {
                  kernelClass.updateKernel( k, q, stack, constitutiveWrapper );

                  kernelClass.stiffnessKernel( k, q, stack );

                  kernelClass.integrationKernel( k, q, stack, constitutiveWrapper );
                }
                maxForce.max( kernelClass.postKernel( stack ) );
              }
            } );
      } );
    } );
  } );

  return maxForce.get();
}
}
}



#endif /* GEOSX_PHYSICSSOLVERS_PHYSICSLOOPINTERFACE_HPP_ */
