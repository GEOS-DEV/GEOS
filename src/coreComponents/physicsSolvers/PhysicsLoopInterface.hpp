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



class FiniteElementRegionLoop
{
public:

  struct Parameters
  {};

  template< int NUM_NODES_PER_ELEM, int NUM_DOF_PER_NODE >
  struct StackVariables
  {
  public:
    static constexpr int numNodesPerElem = NUM_NODES_PER_ELEM;
    static constexpr int numDofPerNode= NUM_DOF_PER_NODE;
    static constexpr int ndof = NUM_DOF_PER_NODE * NUM_NODES_PER_ELEM;

    //    GEOSX_HOST_DEVICE
    StackVariables():
      elementLocalDofIndex{0},
      localResidual{0.0},
      localJacobian{{0.0}}
      {}

      globalIndex elementLocalDofIndex[ndof];
      real64      localResidual[ndof];
      real64      localJacobian[ndof][ndof];
  };

  template< typename SUBREGION_TYPE,
  typename CONSTITUTIVE_TYPE >
  class Kernels
  {
  public:

    Kernels( arrayView1d< globalIndex const > const & inputDofNumber,
             ParallelMatrix & inputMatrix,
             ParallelVector & inputRhs,
             NodeManager const & nodeManager,
             SUBREGION_TYPE & elementSubRegion,
             CONSTITUTIVE_TYPE & inputConstitutiveType ):
      Kernels( inputDofNumber,
               inputMatrix,
               inputRhs,
               nodeManager,
               elementSubRegion,
               inputConstitutiveType,
               typename CONSTITUTIVE_TYPE::KernelWrapper() )
    {}

    template< typename STACK_VARIABLE_TYPE >
    //  GEOSX_HOST_DEVICE
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

    template< typename STACK_VARIABLE_TYPE >
    //  GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void updateKernel( localIndex const GEOSX_UNUSED_PARAM( k ),
                       localIndex const GEOSX_UNUSED_PARAM( q ),
                       STACK_VARIABLE_TYPE & GEOSX_UNUSED_PARAM( stack ) ) const
    {}

    template< typename PARAMETERS_TYPE, typename STACK_VARIABLE_TYPE >
    //  GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void stiffnessKernel( localIndex const GEOSX_UNUSED_PARAM( k ),
                          localIndex const GEOSX_UNUSED_PARAM( q ),
                          PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM( parameters ),
                          STACK_VARIABLE_TYPE & GEOSX_UNUSED_PARAM( stack ) ) const
    {}

    template< typename PARAMETERS_TYPE, typename STACK_VARIABLE_TYPE >
    //  GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void integrationKernel( localIndex const GEOSX_UNUSED_PARAM( k ),
                            localIndex const GEOSX_UNUSED_PARAM( q ),
                            PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM( parameters ),
                            STACK_VARIABLE_TYPE & GEOSX_UNUSED_PARAM( stack ) ) const
    {}

    template< typename PARAMETERS_TYPE, typename STACK_VARIABLE_TYPE >
    //  GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    real64 postKernel( PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM( parameters ),
                       STACK_VARIABLE_TYPE & stack ) const
    {
      m_matrix.insert( stack.elementLocalDofIndex,
                       stack.elementLocalDofIndex,
                       &(stack.localJacobian[0][0]),
                       stack.ndof,
                       stack.ndof );
      return 0;
    }

    arrayView1d< globalIndex const > const m_dofNumber;
    ParallelMatrix & m_matrix;
    ParallelVector & m_rhs;
    typename SUBREGION_TYPE::NodeMapType::base_type::ViewTypeConst const elemsToNodes;
    arrayView1d< integer const > const elemGhostRank;
    typename CONSTITUTIVE_TYPE::KernelWrapper const constitutiveUpdate;

protected:
    Kernels( arrayView1d< globalIndex const > const & inputDofNumber,
             ParallelMatrix & inputMatrix,
             ParallelVector & inputRhs,
             NodeManager const & GEOSX_UNUSED_PARAM(nodeManager),
             SUBREGION_TYPE & elementSubRegion,
             CONSTITUTIVE_TYPE & GEOSX_UNUSED_PARAM(inputConstitutiveType),
             typename CONSTITUTIVE_TYPE::KernelWrapper const & inputConstitutiveUpdate ):
     m_dofNumber( inputDofNumber ),
     m_matrix( inputMatrix ),
     m_rhs( inputRhs ),
     elemsToNodes( elementSubRegion.nodeList() ),
     elemGhostRank( elementSubRegion.ghostRank() ),
     constitutiveUpdate( inputConstitutiveUpdate )
  {}

  };

  template< typename POLICY,
  int NUM_NODES_PER_ELEM,
  int NUM_QUADRATURE_POINTS,
  typename STACK_VARIABLES,
  typename PARAMETERS_TYPE,
  typename KERNEL_CLASS >
  static
  real64 Launch( localIndex const numElems,
                 PARAMETERS_TYPE const & parameters,
                 KERNEL_CLASS const & kernelClass )
  {
    RAJA::ReduceMax< serialReduce, real64 > maxResidual( 0 );

    forAll< POLICY >( numElems,
                      [=] ( localIndex const k )
                      {
      STACK_VARIABLES stack;

      if( kernelClass.elemGhostRank[k] < 0 )
      {
        kernelClass.preKernel( k, stack );
        for( integer q=0; q<NUM_QUADRATURE_POINTS; ++q )
        {
          kernelClass.updateKernel( k, q, stack );

          kernelClass.stiffnessKernel( k, q, parameters, stack );

          kernelClass.integrationKernel( k, q, parameters, stack );
        }
        maxResidual.max( kernelClass.postKernel( parameters, stack ) );
      }
                      } );
    return maxResidual.get();
  }


  template< typename POLICY,
            typename UPDATE_CLASS,
            typename ... KERNEL_PARAMS >
  static
  real64 Execute( MeshLevel & mesh,
                  string_array const & targetRegions,
                  string const & solidMaterialName,
                  FiniteElementDiscretization const * const feDiscretization,
                  arrayView1d< globalIndex const > const & inputDofNumber,
                  ParallelMatrix & inputMatrix,
                  ParallelVector & inputRhs,
                  typename UPDATE_CLASS::Parameters const & parameters )
  {

    real64 maxResidual = 0;

    NodeManager const & nodeManager = *(mesh.getNodeManager());
    ElementRegionManager & elementRegionManager = *(mesh.getElemManager());


    elementRegionManager.forElementSubRegions<CellElementSubRegion>( targetRegions,
                                                                     [&] ( auto & elementSubRegion )
    {
      localIndex const numElems = elementSubRegion.size();
      typedef TYPEOFREF( elementSubRegion ) SUBREGIONTYPE;

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
        solidBase = elementSubRegion.template getConstitutiveModel<constitutive::SolidBase>( solidMaterialName );

        constitutive::constitutiveUpdatePassThru( solidBase, [&]( auto & castedConstitutiveRelation )
        {
          using CONSTITUTIVE_TYPE = TYPEOFREF( castedConstitutiveRelation );

          typename UPDATE_CLASS::template Kernels<SUBREGIONTYPE,CONSTITUTIVE_TYPE> kernelClass( inputDofNumber,
                                                                              inputMatrix,
                                                                              inputRhs,
                                                                              nodeManager,
                                                                              elementSubRegion,
                                                                              castedConstitutiveRelation );

          maxResidual = std::max( maxResidual,
                                  Launch< POLICY,
                                          NUM_NODES_PER_ELEM,
                                          NUM_QUADRATURE_POINTS,
                                          typename UPDATE_CLASS::template StackVariables<NUM_NODES_PER_ELEM,3>>( numElems,
                                                                                                                 parameters,
                                                                                                                 kernelClass ) );
        } );
      } );
    } );

    return maxResidual;
  }

  template< typename POLICY, typename ... PARAMS >
  static
  real64 FillSparsity( PARAMS && ... params )
  {
    return Execute<POLICY, FiniteElementRegionLoop>( std::forward<PARAMS>(params)..., FiniteElementRegionLoop::Parameters() );
  }


};

}
}



#endif /* GEOSX_PHYSICSSOLVERS_PHYSICSLOOPINTERFACE_HPP_ */
