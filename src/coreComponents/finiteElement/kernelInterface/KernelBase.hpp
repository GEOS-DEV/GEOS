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
 * @file KernelBase.hpp
 */

#ifndef GEOSX_FINITEELEMENT_KERNELBASE_HPP_
#define GEOSX_FINITEELEMENT_KERNELBASE_HPP_

#include "common/DataTypes.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "mesh/ElementRegionManager.hpp"


// 0 is the parameter classes
// 1 is std::tuple
// 2 is camp
#define CONSTRUCTOR_PARAM_OPTION 0

#if CONSTRUCTOR_PARAM_OPTION==1
namespace std
{

namespace detail
{
template< class T, class Tuple, std::size_t... I >
constexpr T make_from_tuple_impl( Tuple && t, std::index_sequence< I... > )
{
  return T( std::get< I >( std::forward< Tuple >( t ))... );
}
} // namespace detail

template< class T, class Tuple >
constexpr T make_from_tuple( Tuple && t )
{
  return detail::make_from_tuple_impl< T >( std::forward< Tuple >( t ),
                                            std::make_index_sequence< std::tuple_size< std::remove_reference_t< Tuple > >::value >{} );
}

}
#elif CONSTRUCTOR_PARAM_OPTION==2
  #include "camp/camp.hpp"


namespace camp
{
namespace detail
{
template< class T, class Tuple, idx_t... I >
constexpr T make_from_tuple_impl( Tuple && t, idx_seq< I... > )
{
  return T( get< I >( camp::forward< Tuple >( t ))... );
}
}   // namespace detail

/// Instantiate T from tuple contents, like camp::invoke(tuple,constructor) but
/// functional
template< class T, class Tuple >
constexpr T make_from_tuple( Tuple && tt )
{
  return
    detail::
      make_from_tuple_impl< T >( camp::forward< Tuple >( tt ),
                                 make_idx_seq_t< tuple_size< type::ref::rem< Tuple > >::value >{} );
}
}

#endif

namespace geosx
{

/**
 * @namespace finiteElement Contains the finite element implementation.
 */
namespace finiteElement
{


/**
 * @brief Call a lambda function (callback) with an integral_constant
 *        conversion of the input integral type to allow for static
 *        dispatch.
 * @tparam INTEGRAL_TYPE The type of integer passed in @p input.
 * @tparam LAMBDA The type of @p lambda to execute.
 * @param input The integer to convert to an integral_constant.
 * @param lambda The generic lambda function (takes the integral_constant as
 *               a parameter) that will be executed.
 *
 * Implements a switchyard to convert the value of @p input to an
 * integral_constant<@p INTEGRAL_TYPE, @p input>, and pass that to @p lambda.
 * This allows a runtime @p input to be dispatched as a compile time constant.
 * Note that @p LAMBDA must be a generic lambda that takes in a single `auto`
 * parameter and then converts the value to an INTEGRAL_TYPE. For instance:
 *
 *     int value = 1;
 *     integralTypeDispatch( 1, [&]( auto const constValueType )
 *     {
 *       static constexpr int constValue = decltype( constValueType )::value;
 *
 *       func< constValue >(...);
 *     };
 */
template< typename INTEGRAL_TYPE, typename LAMBDA >
void
integralTypeDispatch( INTEGRAL_TYPE const input,
                      LAMBDA && lambda )
{
  switch( input )
  {
    case 1:
    {
      lambda( std::integral_constant< INTEGRAL_TYPE, 1 >() );
      break;
    }
//    case 4:
//    {
//      lambda( std::integral_constant< INTEGRAL_TYPE, 4 >() );
//      break;
//    }
//    case 5:
//    {
//      lambda( std::integral_constant< INTEGRAL_TYPE, 5 >() );
//      break;
//    }
//    case 6:
//    {
//      lambda( std::integral_constant< INTEGRAL_TYPE, 6 >() );
//      break;
//    }
    case 8:
    {
      lambda( std::integral_constant< INTEGRAL_TYPE, 8 >() );
      break;
    }
    default:
      GEOSX_ERROR( "integralTypeDispatch() is not implemented for value of: "<<input );
  }
}

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
/**
 * @class KernelBase
 * @brief Define the base interface for finite element kernels.
 * @tparam SUBREGION_TYPE The type of subregion that the kernel will act on.
 * @tparam CONSTITUTIVE_TYPE The type of constitutive model present in the
 *                           subregion.
 * @tparam NUM_TEST_SUPPORT_POINTS_PER_ELEM The number of test space support
 *                                          points per element.
 * @tparam NUM_TRIAL_SUPPORT_POINTS_PER_ELEM The number of trial space support
 *                                           points per element.
 * @tparam NUM_DOF_PER_TEST_SP The number of DOF per test support point.
 * @tparam NUM_DOF_PER_TRIAL_SP The number of DOF per trial support point.
 *
 * ### General KernelBase Description
 *
 * KernelBase defines an interface for implementing finite element kernels
 * that will be callable by the family of kernel launch functions. Specific
 * physics kernels may or may not derive from KernelBase, but must follow
 * the same interface in order to be callable from the generic launching
 * functions.
 *
 * The template parameters of KernelBase should be duplicated as part of the
 * interface, EXCEPT for @p NUM_DOF_PER_TEST_SP and @p NUM_DOF_PER_TRIAL_SP.
 * These values should be set internally by the physics solver since each
 * physics discretization will have a constant intrinsic value for these
 * quantities. For example, when solving or the heat equation with scalar
 * temperature as the primary variable at the support point, these will have
 * a value of 1. In contrast, when solving a solid mechanics problem, with
 * vector displacement as the primary variable at the support point, these
 * will have a value of 3. Note that the interface provided by
 * geosx::finiteElement::RegionBasedKernelApplication will construct a
 * kernel assuming only the first 4 template arguments.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          int NUM_TEST_SUPPORT_POINTS_PER_ELEM,
          int NUM_TRIAL_SUPPORT_POINTS_PER_ELEM,
          int NUM_DOF_PER_TEST_SP,
          int NUM_DOF_PER_TRIAL_SP >
class KernelBase
{
public:
  /// Compile time value for the number of test function support points per
  /// element.
  static constexpr int numTestSupportPointsPerElem  = NUM_TEST_SUPPORT_POINTS_PER_ELEM;

  /// Compile time value for the number of trial function support points per
  /// element.
  static constexpr int numTrialSupportPointsPerElem = NUM_TRIAL_SUPPORT_POINTS_PER_ELEM;

  /// Compile time value for the number of degrees of freedom per test function
  /// support point.
  static constexpr int numDofPerTestSupportPoint    = NUM_DOF_PER_TEST_SP;

  /// Compile time value for the number of degrees of freedom per trial
  /// function support point.
  static constexpr int numDofPerTrialSupportPoint   = NUM_DOF_PER_TRIAL_SP;

  /**
   * @brief Constructor
   * @param elementSubRegion Reference to the SUBREGION_TYPE(class template
   *                         parameter) object.
   * @param inputConstitutiveUpdate The constitutive update object.
   * @param finiteElementSpace Placeholder for the finite element space object,
   *                           which currently doesn't do much.
   */
  KernelBase( SUBREGION_TYPE const & elementSubRegion,
              FiniteElementBase const * const finiteElementSpace,
              CONSTITUTIVE_TYPE * const inputConstitutiveType ):
    m_elemsToNodes( elementSubRegion.nodeList().toViewConst() ),
    m_elemGhostRank( elementSubRegion.ghostRank() ),
    m_constitutiveUpdate( inputConstitutiveType->createKernelWrapper() ),
    m_finiteElementSpace( finiteElementSpace )
  {}



  struct StackVariables
  {};

  /**
   * @brief Performs the setup phase for the kernel.
   * @tparam STACK_VARIABLE_TYPE The type of StackVariable that holds the stack
   *                             variables. This is most likely a defined in a
   *                             type that derives from KernelBase.
   * @param k The element index.
   * @param stack The StackVariable object that hold the stack variables.
   *
   * ### KernelBase::setup() Description
   *
   * The operations typically found in setup are thing such as the collection
   * of global data into local stack storage.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );
    GEOSX_UNUSED_VAR( stack );
  }

  /**
   * @brief Performs a state update at a quadrature point.
   * @tparam STACK_VARIABLE_TYPE The type of StackVariable that holds the stack
   *                             variables. This is most likely a defined in a
   *                             type that derives from KernelBase.
   * @param k The element index.
   * @param q The quadrature point index.
   * @param stack The StackVariable object that hold the stack variables.
   *
   * ### KernelBase::quadraturePointStateUpdate() Description
   *
   * The operations found here are the mapping from the support points to the
   * quadrature point, calculation of gradients, etc. From this data the
   * state of the constitutive model is updated if required by the physics
   * package.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointStateUpdate( localIndex const k,
                                   localIndex const q,
                                   StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );
    GEOSX_UNUSED_VAR( q );
    GEOSX_UNUSED_VAR( stack );
  }

  /**
   * @brief Form the element local Jacobian matrix.
   * @tparam STACK_VARIABLE_TYPE The type of StackVariable that holds the stack
   *                             variables. This is most likely a defined in a
   *                             type that derives from KernelBase.
   * @param k The element index.
   * @param q The quadrature point index.
   * @param stack The StackVariable object that hold the stack variables.
   *
   * ### KernelBase::quadraturePointJacobianContribution() Description
   *
   * The results of quadraturePointStateUpdate are used to form the local
   * element Jacobian matrix.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointJacobianContribution( localIndex const k,
                                            localIndex const q,
                                            StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );
    GEOSX_UNUSED_VAR( q );
    GEOSX_UNUSED_VAR( stack );
  }

  /**
   * @brief Calculates the element local Residual vector.
   * @tparam STACK_VARIABLE_TYPE The type of StackVariable that holds the stack
   *                             variables. This is most likely a defined in a
   *                             type that derives from KernelBase.
   * @param k The element index.
   * @param q The quadrature point index.
   * @param stack The StackVariable object that hold the stack variables.
   *
   * ### KernelBase::quadraturePointResidualContribution() Description
   *
   * The results of quadraturePointStateUpdate are used to form the local
   * element residual vector.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointResidualContribution( localIndex const k,
                                            localIndex const q,
                                            StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );
    GEOSX_UNUSED_VAR( q );
    GEOSX_UNUSED_VAR( stack );
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @tparam STACK_VARIABLE_TYPE The type of StackVariable that holds the stack
   *                             variables. This is most likely a defined in a
   *                             type that derives from KernelBase.
   * @param k The element index.
   * @param stack The StackVariable object that hold the stack variables.
   *
   * ### KernelBase::complete() Description
   *
   * The operations typically found in complete are the mapping of the local
   * Jacobian and Residual into the global Jacobian and Residual.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );
    GEOSX_UNUSED_VAR( stack );
    return 0;
  }


  /**
   * @brief Kernel Launcher.
   * @tparam POLICY The RAJA policy to use for the launch.
   * @tparam NUM_QUADRATURE_POINTS The number of quadrature points per element.
   * @tparam KERNEL_TYPE The type of Kernel to execute.
   * @param numElems The number of elements to process in this launch.
   * @param kernelComponent The instantiation of KERNEL_TYPE to execute.
   * @return The maximum residual.
   *
   * This is a generic launching function for all of the finite element kernels
   * that follow the interface set by KernelBase.
   */
  template< typename POLICY,
            int NUM_QUADRATURE_POINTS,
            typename KERNEL_TYPE >
  static
  typename std::enable_if< std::is_same< POLICY, serialPolicy >::value ||
                           std::is_same< POLICY, parallelHostPolicy >::value, real64 >::type
  Launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOSX_MARK_FUNCTION;
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxResidual( 0 );

    forAll< POLICY >( numElems,
                      [=] ( localIndex const k )
    {
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( k, stack );
      for( integer q=0; q<NUM_QUADRATURE_POINTS; ++q )
      {
        kernelComponent.quadraturePointStateUpdate( k, q, stack );

        kernelComponent.quadraturePointJacobianContribution( k, q, stack );

        kernelComponent.quadraturePointResidualContribution( k, q, stack );
      }
      if( kernelComponent.m_elemGhostRank[k] < 0 )
      {
        maxResidual.max( kernelComponent.complete( k, stack ) );
      }
    } );
    return maxResidual.get();
  }

  /**
   * @brief Kernel Launcher.
   * @tparam POLICY The RAJA policy to use for the launch.
   * @tparam NUM_QUADRATURE_POINTS The number of quadrature points per element.
   * @tparam KERNEL_TYPE The type of Kernel to execute.
   * @param numElems The number of elements to process in this launch.
   * @param kernelComponent The instantiation of KERNEL_TYPE to execute.
   * @return The maximum residual.
   *
   * This is a generic launching function for all of the finite element kernels
   * that follow the interface set by KernelBase.
   */
  template< typename POLICY,
            int NUM_QUADRATURE_POINTS,
            typename KERNEL_TYPE >
  static
  typename std::enable_if< !( std::is_same< POLICY, serialPolicy >::value ||
                              std::is_same< POLICY, parallelHostPolicy >::value ), real64 >::type
  Launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOSX_MARK_FUNCTION;
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxResidual( 0 );

    forAll< POLICY >( numElems,
                      [=] GEOSX_DEVICE ( localIndex const k )
    {
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( k, stack );
      for( integer q=0; q<NUM_QUADRATURE_POINTS; ++q )
      {
        kernelComponent.quadraturePointStateUpdate( k, q, stack );

        kernelComponent.quadraturePointJacobianContribution( k, q, stack );

        kernelComponent.quadraturePointResidualContribution( k, q, stack );
      }
      if( kernelComponent.m_elemGhostRank[k] < 0 )
      {
        maxResidual.max( kernelComponent.complete( k, stack ) );
      }
    } );
    return maxResidual.get();
  }

protected:
  /// The element to nodes map.
  typename SUBREGION_TYPE::NodeMapType::base_type::ViewTypeConst const m_elemsToNodes;

  /// The element ghost rank array.
  arrayView1d< integer const > const m_elemGhostRank;

  /// The constitutive update object used to update the constitutive state,
  /// and extract constitutive data.
  typename CONSTITUTIVE_TYPE::KernelWrapper const m_constitutiveUpdate;

  /// The finite element space/discretization object for the element type in
  /// the SUBREGION_TYPE.
  FiniteElementBase const * m_finiteElementSpace;
};


//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

/**
 * @brief Performs a loop over specific regions (by type and name) and calls
 *        a kernel launch on the subregions with compile time knowledge of
 *        sub-loop bounds such as number of nodes and quadrature points per
 *        element.
 * @tparam POLICY The RAJA launch policy to pass to the kernel launch.
 * @tparam CONSTITUTIVE_BASE The common base class for constitutive
 *                           pass-thru/dispatch which gives the kernel launch
 *                           compile time knowledge of the constitutive model.
 *                           This is achieved through a call to the
 *                           #constitutive::ConstitutivePassThru function which
 *                           should have a specialization for CONSTITUTIVE_BASE
 *                           implemented in order to perform the compile time
 *                           dispatch.
 * @tparam REGION_TYPE The type of region to loop over. TODO make this a
 *                     parameter pack?
 * @tparam KERNEL_TEMPLATE The type of template for the physics kernel, which
 *                         conforms to the interface specified by KernelBase.
 * @tparam KERNEL_CONSTRUCTOR_PARAMS The template parameter pack to hold the
 *                                   parameter pack (i.e. custom) parameters
 *                                   that are sent to the constructor for the
 *                                   @p KERNEL_TEMPLATE.
 * @param mesh The MeshLevel object.
 * @param targetRegions The names of the target regions(of type @p REGION_TYPE)
 *                      to apply the @p KERNEL_TEMPLATE.
 * @param constitutiveNames The names of the constitutive models present in the
 *                          Region.
 * @param feDiscretization A pointer to the finite element discretization/space
 *                         object.
 * @param kernelConstructorParams The parameter list for corresponding to the
 *                                parameter @p KERNEL_CONSTRUCTOR_PARAMS that
 *                                are passed to the @p KERNEL_TEMPLATE
 *                                constructor.
 * @return The maximum contribution to the residual, which may be used to scale
 *         the residual.
 *
 * Loops over all regions Applies/Launches a kernel specified by the @p KERNEL_TEMPLATE through
 * #KernelBase::Launch.
 */
template< typename POLICY,
          typename CONSTITUTIVE_BASE,
          typename REGION_TYPE,
          template< typename SUBREGION_TYPE,
                    typename CONSTITUTIVE_TYPE,
                    int NUM_TEST_SUPPORT_POINTS_PER_ELEM,
                    int NUM_TRIAL_SUPPORT_POINTS_PER_ELEM > class KERNEL_TEMPLATE,
          typename ... KERNEL_CONSTRUCTOR_PARAMS >
static
real64 RegionBasedKernelApplication( MeshLevel & mesh,
                                     arrayView1d< string const > const & targetRegions,
                                     arrayView1d< string const > const & constitutiveNames,
                                     FiniteElementDiscretization const * const feDiscretization,
                                     KERNEL_CONSTRUCTOR_PARAMS && ... kernelConstructorParams )
{

  real64 maxResidualContribution = 0;

  NodeManager & nodeManager = *(mesh.getNodeManager());
  EdgeManager & edgeManager = *(mesh.getEdgeManager());
  FaceManager & faceManager = *(mesh.getFaceManager());
  ElementRegionManager & elementRegionManager = *(mesh.getElemManager());


#if CONSTRUCTOR_PARAM_OPTION==0
  using CONSTRUCTOR_PARAMS = typename KERNEL_TEMPLATE< CellElementSubRegion,
                                                       constitutive::NullModel,
                                                       1,
                                                       1 >::ConstructorParams;

  CONSTRUCTOR_PARAMS kernelParams( std::forward< KERNEL_CONSTRUCTOR_PARAMS >( kernelConstructorParams )... );

#elif CONSTRUCTOR_PARAM_OPTION==1
  auto kernelConstructorParamsTuple = std::forward_as_tuple( std::forward< KERNEL_CONSTRUCTOR_PARAMS >( kernelConstructorParams )... );
#elif CONSTRUCTOR_PARAM_OPTION==2
  auto kernelConstructorParamsTuple = camp::forward_as_tuple( std::forward< KERNEL_CONSTRUCTOR_PARAMS >( kernelConstructorParams )... );
#endif


  elementRegionManager.forElementSubRegions< REGION_TYPE >( targetRegions,
                                                            [&] ( localIndex const targetRegionIndex,
                                                                  auto & elementSubRegion )
  {
    localIndex const numElems = elementSubRegion.size();
    typedef TYPEOFREF( elementSubRegion ) SUBREGIONTYPE;


    FiniteElementBase const * finiteElementSpace = nullptr;
    if( feDiscretization != nullptr )
    {
      finiteElementSpace = ( feDiscretization->getFiniteElement( elementSubRegion.GetElementTypeString() ) ).get();
    }

    localIndex const
    numQuadraturePointsPerElem = finiteElementSpace == nullptr ?
                                 1 :
                                 finiteElementSpace->n_quadrature_points();

    constitutive::ConstitutiveBase * constitutiveRelation = nullptr;
    constitutive::NullModel * nullConstitutiveModel = nullptr;
    if( targetRegionIndex <= constitutiveNames.size()-1 )
    {
      constitutiveRelation = elementSubRegion.template getConstitutiveModel( constitutiveNames[targetRegionIndex] );
    }
    else
    {
      nullConstitutiveModel = elementSubRegion.template RegisterGroup< constitutive::NullModel >( "nullModelGroup" );
      constitutiveRelation = nullConstitutiveModel;
    }

    constitutive::ConstitutivePassThru< CONSTITUTIVE_BASE >::Execute( constitutiveRelation,
                                                                      [&]( auto * const castedConstitutiveRelation )
    {
      using CONSTITUTIVE_TYPE = TYPEOFPTR( castedConstitutiveRelation );
      integralTypeDispatch( elementSubRegion.numNodesPerElement(), [&]( auto const NNPE )
      {
        integralTypeDispatch( numQuadraturePointsPerElem, [&]( auto const NQPPE )
        {
          static constexpr int NUM_NODES_PER_ELEM = decltype( NNPE )::value;
          static constexpr int NUM_QUADRATURE_POINTS = decltype( NQPPE )::value;


          using KERNEL_TYPE = KERNEL_TEMPLATE< SUBREGIONTYPE,
                                               CONSTITUTIVE_TYPE,
                                               NUM_NODES_PER_ELEM,
                                               NUM_NODES_PER_ELEM >;

#if CONSTRUCTOR_PARAM_OPTION==0
          KERNEL_TYPE kernelComponent( nodeManager,
                                       edgeManager,
                                       faceManager,
                                       elementSubRegion,
                                       finiteElementSpace,
                                       castedConstitutiveRelation,
                                       kernelParams );

#elif CONSTRUCTOR_PARAM_OPTION==1

          auto temp = std::forward_as_tuple( nodeManager,
                                             edgeManager,
                                             faceManager,
                                             elementSubRegion,
                                             finiteElementSpace,
                                             castedConstitutiveRelation );

          auto fullKernelComponentConstructorArgs = std::tuple_cat( temp,
                                                                    kernelConstructorParamsTuple );

          KERNEL_TYPE kernelComponent  = std::make_from_tuple< KERNEL_TYPE >( fullKernelComponentConstructorArgs );

#elif CONSTRUCTOR_PARAM_OPTION==2
          auto temp = camp::forward_as_tuple( nodeManager,
                                              edgeManager,
                                              faceManager,
                                              elementSubRegion,
                                              finiteElementSpace,
                                              castedConstitutiveRelation );
          auto fullKernelComponentConstructorArgs = camp::tuple_cat_pair( temp,
                                                                          kernelConstructorParamsTuple );
          KERNEL_TYPE kernelComponent  = camp::make_from_tuple< KERNEL_TYPE >( fullKernelComponentConstructorArgs );

#endif

          // We can replace all the shenanigans with params when we can use c++20
          // perfect forwarding of the parameter pack in the lambda capture.
          //        KERNEL_TYPE kernelComponent( nodeManager,
          //                                     edgeManager,
          //                                     faceManager,
          //                                     elementSubRegion,
          //                                     finiteElementSpace,
          //                                     castedConstitutiveRelation,
          //                                     std::forward< KERNEL_CONSTRUCTOR_PARAMS >( kernelConstructorParams )...
          // );


          maxResidualContribution = std::max( maxResidualContribution,
                                              KERNEL_TYPE::template Launch< POLICY,
                                                                            NUM_QUADRATURE_POINTS
                                                                            >( numElems,
                                                                               kernelComponent ) );
        } );
      } );
    } );

    if( nullConstitutiveModel )
    {
      elementSubRegion.deregisterGroup( "nullModelGroup" );
    }

  } );

  return maxResidualContribution;
}
}
}



#endif /* GEOSX_FINITEELEMENT_KERNELBASE_HPP_ */
