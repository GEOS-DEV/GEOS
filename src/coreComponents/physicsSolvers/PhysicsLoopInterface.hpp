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
#include "mesh/ElementRegionManager.hpp"


namespace geosx
{
namespace physicsLoopInterface
{

class FiniteElementRegionLoopKernelBase
{
public:

  template< int NDOF, int NUM_NODES_PER_ELEM >
  struct StackVariables
  {
public:
    static constexpr int ndof = NDOF;

    StackVariables():
      elementLocalDofIndex( NDOF ),
      R( NDOF ),
      dRdU( NDOF, NDOF )
    {
      dRdU = 0.0;
      R = 0.0;
    }

    stackArray1d< globalIndex, NDOF >   elementLocalDofIndex;
    stackArray1d< real64, NDOF >        R;
    stackArray2d< real64, NDOF *NDOF >  dRdU;

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
                  localIndex const numNodesPerElem,
                  STACK_VARIABLE_TYPE & stack ) const
  {
    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = elemsToNodes( k, a );
      for( int i=0; i<3; ++i )
      {
        stack.elementLocalDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
      }
    }

  }

  template< typename STACK_VARIABLE_TYPE, typename CONSTITUTIVE_TYPE >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void updateKernel( localIndex const GEOSX_UNUSED_PARAM( k ),
                     localIndex const GEOSX_UNUSED_PARAM( q ),
                     STACK_VARIABLE_TYPE & GEOSX_UNUSED_PARAM( stack ),
                     CONSTITUTIVE_TYPE & GEOSX_UNUSED_PARAM( constitutive ),
                     localIndex const GEOSX_UNUSED_PARAM( numNodesPerElem ) ) const
  {}

  template< typename STACK_VARIABLE_TYPE >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void stiffnessKernel( localIndex const GEOSX_UNUSED_PARAM( k ),
                        localIndex const GEOSX_UNUSED_PARAM( q ),
                        localIndex const GEOSX_UNUSED_PARAM( a ),
                        localIndex const GEOSX_UNUSED_PARAM( b ),
                        STACK_VARIABLE_TYPE & GEOSX_UNUSED_PARAM( stack ) ) const
  {}

  template< typename STACK_VARIABLE_TYPE, typename CONSTITUTIVE_TYPE >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void integrationKernel( localIndex const GEOSX_UNUSED_PARAM( k ),
                          localIndex const GEOSX_UNUSED_PARAM( q ),
                          localIndex const GEOSX_UNUSED_PARAM( numNodesPerElem ),
                          STACK_VARIABLE_TYPE & GEOSX_UNUSED_PARAM( stack ),
                          CONSTITUTIVE_TYPE const & GEOSX_UNUSED_PARAM( constitutive ) ) const
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

template< typename KERNEL_CLASS >
static
real64 FiniteElementRegionLoop( MeshLevel & mesh,
                                string_array const & targetRegions,
                                string const & solidMaterialName,
                                KERNEL_CLASS & kernelClass )
{
  RAJA::ReduceMax< serialReduce, double > maxForce( 0 );

  kernelClass.initializeNonElementViews( mesh );

  ElementRegionManager & elementRegionManager = *(mesh.getElemManager());

  elementRegionManager.forElementSubRegions< CellElementSubRegion >( targetRegions,
                                                                     [&] ( auto & elementSubRegion )
  {
    localIndex const numElems = elementSubRegion.size();
//      typedef TYPEOFREF( elementSubRegion ) SUBREGIONTYPE;

    localIndex const NUM_NODES_PER_ELEM = 8;  //elementSubRegion.numNodesPerElement();
    localIndex const NUM_QUADRATURE_POINTS = 8;
    localIndex const ndof = 3 * NUM_NODES_PER_ELEM;

    arrayView1d< integer const > const & elemGhostRank = elementSubRegion.ghostRank();

    constitutive::LinearElasticIsotropic * const
    material = elementSubRegion.GetConstitutiveModels()->template GetGroup< constitutive::LinearElasticIsotropic >( solidMaterialName );

    constitutive::LinearElasticIsotropicUpdates constitutiveWrapper = material->createKernelWrapper();

    kernelClass.initializeElementSubRegionViews( elementSubRegion );


    RAJA::forall< serialPolicy >( RAJA::TypedRangeSegment< localIndex >( 0, numElems ),
                                  [=] GEOSX_DEVICE ( localIndex const k )
        {

          typename KERNEL_CLASS:: template StackVariables< ndof, NUM_NODES_PER_ELEM > stack;

          if( elemGhostRank[k] < 0 )
          {
            kernelClass.preKernel( k, NUM_NODES_PER_ELEM, stack );

            for( integer q=0; q<NUM_QUADRATURE_POINTS; ++q )
            {
              kernelClass.updateKernel( k, q, stack, constitutiveWrapper, NUM_NODES_PER_ELEM );

              for( localIndex a=0; a<NUM_NODES_PER_ELEM; ++a )
              {
                for( localIndex b=0; b<NUM_NODES_PER_ELEM; ++b )
                {
                  kernelClass.stiffnessKernel( k, q, a, b, stack );
                }
              }

              kernelClass.integrationKernel( k,
                                             q,
                                             NUM_NODES_PER_ELEM,
                                             stack,
                                             constitutiveWrapper );
            }

            maxForce.max( kernelClass.postKernel( stack ) );

          }
        } );
  } );

  return maxForce.get();
}
}
}



#endif /* GEOSX_PHYSICSSOLVERS_PHYSICSLOOPINTERFACE_HPP_ */
