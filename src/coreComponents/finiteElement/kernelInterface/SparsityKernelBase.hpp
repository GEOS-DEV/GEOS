/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "ImplicitKernelBase.hpp"

/**
 * @file SparsityKernelBase.hpp
 */

#ifndef GEOSX_FINITEELEMENT_SPARSITYKERNELBASE_HPP_
#define GEOSX_FINITEELEMENT_SPARSITYKERNELBASE_HPP_



namespace geosx
{

namespace finiteElement
{

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
/**
 * @class SparsityKernelBase
 * @brief Define the base interface for implicit finite element kernels.
 * @copydoc geosx::finiteElement::KernelBase
 *
 * ### SparsityKernelBase Description
 * Provide common kernels for generation of CRS Sparsity patterns.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE,
          int NUM_DOF_PER_TEST_SP,
          int NUM_DOF_PER_TRIAL_SP >
class SparsityKernelBase : public ImplicitKernelBase< SUBREGION_TYPE,
                                                      CONSTITUTIVE_TYPE,
                                                      FE_TYPE,
                                                      NUM_DOF_PER_TEST_SP,
                                                      NUM_DOF_PER_TRIAL_SP >
{
public:
  /// Alias for the base class. (i.e. #geosx::finiteElement::ImplicitKernelBase)
  using Base = ImplicitKernelBase< SUBREGION_TYPE,
                                   CONSTITUTIVE_TYPE,
                                   FE_TYPE,
                                   NUM_DOF_PER_TEST_SP,
                                   NUM_DOF_PER_TRIAL_SP >;


  using typename Base::StackVariables;
  using Base::m_dofRankOffset;

  using Base::setup;

  /**
   * @brief Constructor
   * @param nodeManager Reference to the NodeManager object.
   * @param edgeManager Reference to the EdgeManager object.
   * @param faceManager Reference to the FaceManager object.
   * @param targetRegionIndex Index of the region the subregion belongs to.
   * @param inputDofNumber The dof number for the primary field.
   * @param rankOffset dof index offset of current rank
   * @param inputSparsity The sparsity pattern to fill.
   * @copydoc geosx::finiteElement::KernelBase::KernelBase
   */
  SparsityKernelBase( NodeManager const & nodeManager,
                      EdgeManager const & edgeManager,
                      FaceManager const & faceManager,
                      localIndex const targetRegionIndex,
                      SUBREGION_TYPE const & elementSubRegion,
                      FE_TYPE const & finiteElementSpace,
                      CONSTITUTIVE_TYPE & inputConstitutiveType,
                      arrayView1d< globalIndex const > const & inputDofNumber,
                      globalIndex const rankOffset,
                      SparsityPattern< globalIndex > & inputSparsity ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDofNumber,
          rankOffset,
          CRSMatrixView< real64, globalIndex const >(),
          arrayView1d< real64 >() ),
    m_sparsity( inputSparsity )
  {}


  /**
   * @copydoc geosx::finiteElement::KernelBase::complete
   *
   * In this implementation, only the matrix values are inserted, making this
   * implementation appropriate for generating the sparsity pattern.
   */
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );

    for( localIndex r=0; r<stack.numRows; ++r )
    {
      localIndex const row = stack.localRowDofIndex[r] - m_dofRankOffset;
      if( row < 0 || row >= m_sparsity.numRows() ) continue;
      for( localIndex c=0; c<stack.numCols; ++c )
      {
        m_sparsity.insertNonZero( row,
                                  stack.localColDofIndex[c] );
      }
    }
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
            typename KERNEL_TYPE >
  static
  real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    GEOSX_MARK_FUNCTION;

    // launch the kernel
    forAll< POLICY >( numElems,
                      [=] ( localIndex const k )
    {
      // allocate the stack variables
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( k, stack );

      kernelComponent.complete( k, stack );

    } );
    return 0;
  }
private:
  SparsityPattern< globalIndex > & m_sparsity;
};


//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

template < typename POLICY, 
           typename SUBREGION_TYPE, 
           typename CONSTITUTIVE_TYPE,
           typename FE_TYPE, 
           template < typename, typename, typename > class KERNEL_TEMPLATE >
real64 buildSparsityAndInvoke( localIndex const numElems,
                               NodeManager const & nodeManager, 
                               EdgeManager const & edgeManager,
                               FaceManager const & faceManager,
                               localIndex const targetRegionIndex,
                               SUBREGION_TYPE const & elementSubRegion,
                               FE_TYPE const & finiteElementSpace,
                               CONSTITUTIVE_TYPE & inputConstitutiveType,
                               arrayView1d< globalIndex const > const & inputDofNumber,
                               globalIndex const rankOffset,
                               SparsityPattern< globalIndex > & inputSparsityPattern )
{
    // this requires the underlying kernel be fully specified.. which should result in it compiling, thus sparsity needs to be jitted as well
    // if we could move the numDofPerXXX out of the kernel classes and have a different mechanism to specify them without full specifying the kernel
    //  class we could get away without jitting the sparsity as well
    using Kernel = KERNEL_TEMPLATE< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >;
    using SPARSITY_KERNEL_TYPE = SparsityKernelBase< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE, Kernel::numDofPerTestSupportPoint, Kernel::numDofPerTrialSupportPoint >;
    SPARSITY_KERNEL_TYPE kernel ( nodeManager,
                                  edgeManager,
                                  faceManager,
                                  targetRegionIndex,
                                  elementSubRegion,
                                  finiteElementSpace,
                                  inputConstitutiveType,
                                  inputDofNumber,
                                  rankOffset,
                                  inputSparsityPattern );
    return SPARSITY_KERNEL_TYPE::template kernelLaunch< POLICY, SPARSITY_KERNEL_TYPE >( numElems, kernel );
}

/**
 * @brief Helper struct to define a specialization of
 *   #::geosx::finiteElement::SparsityKernelBase that may be used to generate the sparsity pattern.
 * @tparam KERNEL_TEMPLATE Templated class that defines the physics kernel.
 *   Most likely derives from SparsityKernelBase.
 */
template< template< typename,
                    typename,
                    typename > class KERNEL_TEMPLATE >
class SparsityKernelDispatchTemplate
{
public:

  /**
   * @brief Constructor.
   * @param inputDofNumber An array containing the input degree of freedom numbers.
   * @param rankOffset The global rank offset.
   * @param inputSparsityPattern The local sparsity pattern.
   */
  SparsityKernelDispatchTemplate( arrayView1d< globalIndex const > const & inputDofNumber,
                                  globalIndex const rankOffset,
                                  SparsityPattern< globalIndex > & inputSparsityPattern ):
    m_inputDofNumber( inputDofNumber ),
    m_rankOffset( rankOffset ),
    m_inputSparsityPattern( inputSparsityPattern )
  {}

  /**
   * @brief Return a new instance of @c SparsityKernelBase specialized for @c KERNEL_TEMPLATE.
   * @tparam SUBREGION_TYPE The type of of @p elementSubRegion.
   * @tparam CONSTITUTIVE_TYPE The type of @p inputConstitutiveType.
   * @tparam FE_TYPE The type of @p finiteElementSpace.
   * @param nodeManager The node manager.
   * @param edgeManager The edge manager.
   * @param faceManager The face manager.
   * @param targetRegionIndex The target region index.
   * @param elementSubRegion The sub region on which to generate the sparsity.
   * @param finiteElementSpace The finite element space.
   * @param inputConstitutiveType The constitutive relation.
   * @return A new instance of @c SparsityKernelBase specialized for @c KERNEL_TEMPLATE.
   */
  template< typename POLICY, 
            typename SUBREGION_TYPE,
            typename CONSTITUTIVE_TYPE,
            typename FE_TYPE >
  real64 invoke( localIndex const numElems,
                 NodeManager const & nodeManager,
                 EdgeManager const & edgeManager,
                 FaceManager const & faceManager,
                 localIndex const targetRegionIndex,
                 SUBREGION_TYPE const & elementSubRegion,
                 FE_TYPE const & finiteElementSpace,
                 CONSTITUTIVE_TYPE & inputConstitutiveType )
  {
    return buildSparsityAndInvoke< POLICY, 
                                   SUBREGION_TYPE,
                                   CONSTITUTIVE_TYPE, 
                                   FE_TYPE,
                                   KERNEL_TEMPLATE >( numElems,
                                                      nodeManager,
                                                      edgeManager,
                                                      faceManager,
                                                      targetRegionIndex,
                                                      elementSubRegion,
                                                      finiteElementSpace,
                                                      inputConstitutiveType,
                                                      m_inputDofNumber,
                                                      m_rankOffset,
                                                      m_inputSparsityPattern );
  }

private:

  /// The input degree of freedom numbers.
  arrayView1d< globalIndex const > const & m_inputDofNumber;
  /// The global rank offset.
  globalIndex const m_rankOffset;
  /// The local sparsity pattern.
  SparsityPattern< globalIndex > & m_inputSparsityPattern;
};


jitti::CompilationInfo getSparsityCompilationInfo( const string & );

/**
 * @brief Helper struct to define a specialization of
 *   #::geosx::finiteElement::SparsityKernelBase that may be used to generate the sparsity pattern.
 * @tparam KERNEL_TEMPLATE Templated class that defines the physics kernel.
 *   Most likely derives from SparsityKernelBase.
 */
template < const char * NAME, const char * HEADER >
class SparsityKernelDispatchJIT
{
public:

  /**
   * @brief Constructor.
   * @param inputDofNumber An array containing the input degree of freedom numbers.
   * @param rankOffset The global rank offset.
   * @param inputSparsityPattern The local sparsity pattern.
   */
  SparsityKernelDispatchJIT( arrayView1d< globalIndex const > const & inputDofNumber,
                             globalIndex const rankOffset,
                             SparsityPattern< globalIndex > & inputSparsityPattern ):
    m_inputDofNumber( inputDofNumber ),
    m_rankOffset( rankOffset ),
    m_inputSparsityPattern( inputSparsityPattern )
  {}

  /**
   * @brief Return a new instance of @c SparsityKernelBase specialized for @c KERNEL_TEMPLATE.
   * @tparam SUBREGION_TYPE The type of of @p elementSubRegion.
   * @tparam CONSTITUTIVE_TYPE The type of @p inputConstitutiveType.
   * @tparam FE_TYPE The type of @p finiteElementSpace.
   * @param nodeManager The node manager.
   * @param edgeManager The edge manager.
   * @param faceManager The face manager.
   * @param targetRegionIndex The target region index.
   * @param elementSubRegion The sub region on which to generate the sparsity.
   * @param finiteElementSpace The finite element space.
   * @param inputConstitutiveType The constitutive relation.
   * @return A new instance of @c SparsityKernelBase specialized for @c KERNEL_TEMPLATE.
   */
  template< typename POLICY, 
            typename SUBREGION_TYPE,
            typename CONSTITUTIVE_TYPE,
            typename FE_TYPE >
  real64 invoke( localIndex const numElems,
                 NodeManager const & nodeManager,
                 EdgeManager const & edgeManager,
                 FaceManager const & faceManager,
                 localIndex const targetRegionIndex,
                 SUBREGION_TYPE const & elementSubRegion,
                 FE_TYPE const & finiteElementSpace,
                 CONSTITUTIVE_TYPE & inputConstitutiveType )
  {
    string header( HEADER );
    jitti::CompilationInfo info = getSparsityCompilationInfo( header );
    info.templateParams = LvArray::system::demangleType< POLICY >() + ", " +
                          LvArray::system::demangleType< SUBREGION_TYPE >() + ", " + 
                          LvArray::system::demangleType< CONSTITUTIVE_TYPE >() + ", " +
                          LvArray::system::demangleType< FE_TYPE >() + ", " +
                          string( NAME );
    // Unfortunately can't just decltype(&buildSparsityAndInvoke) since we can't fully specify the function template
    using JIT_SPARSITY_DISPATCH = real64(*)( localIndex const,
                                             NodeManager const &, 
                                             EdgeManager const &,
                                             FaceManager const &,
                                             localIndex const,
                                             SUBREGION_TYPE const &,
                                             FE_TYPE const &,
                                             CONSTITUTIVE_TYPE &,
                                             arrayView1d< globalIndex const > const &,
                                             globalIndex const,
                                             SparsityPattern< globalIndex > & );

    string outputDir( STRINGIZE( JITTI_OUTPUT_DIR ) );
    outputDir += "/";
    static jitti::Cache< JIT_SPARSITY_DISPATCH > buildCache( time(NULL), outputDir );
    if ( MpiWrapper::commRank( ) == 0 )
    {
      buildCache.getOrLoadOrCompile( info );
    }
    MpiWrapper::barrier( );
    // check if the library with the function is available
    if ( buildCache.tryGet( info ) == nullptr )
    {
      // if not, refresh by reading the filesystem to find the new library on the first iteration
      buildCache.refresh( );
    }
    auto & jitSparsityDispatch = buildCache.getOrLoad( info );
    return jitSparsityDispatch( numElems,
                                nodeManager,
                                edgeManager,
                                faceManager,
                                targetRegionIndex,
                                elementSubRegion,
                                finiteElementSpace,
                                inputConstitutiveType,
                                m_inputDofNumber,
                                m_rankOffset,
                                m_inputSparsityPattern );

  }

private:

  /// The input degree of freedom numbers.
  arrayView1d< globalIndex const > const & m_inputDofNumber;
  /// The global rank offset.
  globalIndex const m_rankOffset;
  /// The local sparsity pattern.
  SparsityPattern< globalIndex > & m_inputSparsityPattern;
};

#if JITTI == 1
  #define SparsityKernelDispatch SparsityKernelDispatchJIT
#else
  #define SparsityKernelDispatch SparsityKernelDispatchTemplate
#endif

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
/**
 * @brief Fills matrix sparsity.
 * @tparam REGION_TYPE The type of region to loop over. TODO make this a
 *                     parameter pack?
 * @tparam KERNEL_TEMPLATE The type of template for the physics kernel, which
 *                         conforms to the interface specified by KernelBase.
 * @param mesh The MeshLevel object.
 * @param targetRegions The names of the target regions(of type @p REGION_TYPE)
 *                      to apply the @p KERNEL_TEMPLATE.
 * @param discretizationName The name of the finite element discretization.
 * @param inputDofNumber The global degree of freedom numbers.
 * @param rankOffset Offset of dof indices on curren rank.
 * @param inputSparsityPattern The local sparsity pattern to fill.
 * @return 0
 *
 * Fills matrix sparsity using information from physics specific implementation
 * of #geosx::finiteElement::KernelBase interface using the
 * #geosx::finiteElement::ImplicitKernelBase alias to specialize
 * #geosx::finiteElement::SparsityKernelBase to conform with the template
 * pattern specified in the physics kernels.
 */
template< typename REGION_TYPE,
          typename KERNEL_WRAPPER >
static
real64 fillSparsity( MeshLevel & mesh,
                     arrayView1d< string const > const & targetRegions,
                     string const & discretizationName,
                     arrayView1d< globalIndex const > const & inputDofNumber,
                     globalIndex const rankOffset,
                     SparsityPattern< globalIndex > & inputSparsityPattern )
{
  GEOSX_MARK_FUNCTION;

  KERNEL_WRAPPER kernelDispatch( inputDofNumber, rankOffset, inputSparsityPattern );

  regionBasedKernelApplication< serialPolicy,
                                constitutive::NullModel,
                                REGION_TYPE >( mesh,
                                               targetRegions,
                                               discretizationName,
                                               arrayView1d< string const >(),
                                               kernelDispatch );

  return 0;
}

}
}



#endif /* GEOSX_FINITEELEMENT_SPARSITYKERNELBASE_HPP_ */
