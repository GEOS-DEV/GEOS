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
 * @file SolidMechanicsLagrangianFEMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_KERNELS_FEM_HPP_
#define GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_KERNELS_FEM_HPP_

#include "finiteElement/kernelInterface/RegionLoopSparsity.hpp"

namespace geosx
{

class LaplaceFEMKernel : public finiteElement::RegionLoopSparsity
{
public:
  using Base = finiteElement::RegionLoopSparsity;
  static constexpr int numTestDofPerSP = 1;
  static constexpr int numTrialDofPerSP = 1;

  //***************************************************************************
  /**
   * @class Parameters
   */
  struct Parameters : public Base::Parameters
  {
    Parameters( string const & fieldName ):
      Base::Parameters(),
      m_fieldName{'\0'}
    {
      fieldName.copy( m_fieldName, fieldName.size() );
    }

    char m_fieldName[100];
  };


  //***************************************************************************
  template< int NUM_TEST_SUPPORT_POINTS_PER_ELEM,
            int NUM_TRIAL_SUPPORT_POINTS_PER_ELEM >
  struct StackVariables : Base::StackVariables< NUM_TEST_SUPPORT_POINTS_PER_ELEM*numTestDofPerSP,
                                                NUM_TRIAL_SUPPORT_POINTS_PER_ELEM*numTrialDofPerSP >
  {
public:
    using StackVariablesBase = Base::StackVariables< NUM_TEST_SUPPORT_POINTS_PER_ELEM*numTestDofPerSP,
                                                     NUM_TRIAL_SUPPORT_POINTS_PER_ELEM*numTrialDofPerSP >;

    using StackVariablesBase::numRows;
    using StackVariablesBase::numCols;
    static constexpr int numNodes = NUM_TEST_SUPPORT_POINTS_PER_ELEM;


    GEOSX_HOST_DEVICE
    StackVariables():
      StackVariablesBase(),
      primaryField_local{ 0.0 }
    {}

    real64 primaryField_local[numNodes];
  };


  //***************************************************************************
  template< typename SUBREGION_TYPE,
            typename CONSTITUTIVE_TYPE,
            int NUM_NODES_PER_ELEM,
            int >
  using SparsityKernels = Base::Kernels< SUBREGION_TYPE,
                                         CONSTITUTIVE_TYPE,
                                         NUM_NODES_PER_ELEM,
                                         NUM_NODES_PER_ELEM,
                                         1,
                                         1 >;

  //***************************************************************************
  template< typename SUBREGION_TYPE,
            typename CONSTITUTIVE_TYPE,
            int NUM_NODES_PER_ELEM,
            int >
  class Kernels : public Base::Kernels< SUBREGION_TYPE,
                                        CONSTITUTIVE_TYPE,
                                        NUM_NODES_PER_ELEM,
                                        NUM_NODES_PER_ELEM,
                                        1,
                                        1 >
  {
public:
    using KernelsBase = Base::Kernels< SUBREGION_TYPE,
                                      CONSTITUTIVE_TYPE,
                                      NUM_NODES_PER_ELEM,
                                      NUM_NODES_PER_ELEM,
                                      1,
                                      1 >;

    static constexpr int numNodesPerElem = NUM_NODES_PER_ELEM;

    using KernelsBase::m_dofNumber;
    using KernelsBase::m_matrix;
    using KernelsBase::m_rhs;
    using KernelsBase::elemsToNodes;
    using KernelsBase::elemGhostRank;
    using KernelsBase::constitutiveUpdate;
    using KernelsBase::m_finiteElementSpace;
    using KernelsBase::Launch;

    using StackVars = StackVariables< numNodesPerElem,
                                      numNodesPerElem >;

    Kernels( arrayView1d< globalIndex const > const & inputDofNumber,
             ParallelMatrix & inputMatrix,
             ParallelVector & inputRhs,
             NodeManager const & nodeManager,
             SUBREGION_TYPE const & elementSubRegion,
             FiniteElementBase const * const finiteElementSpace,
             CONSTITUTIVE_TYPE * const inputConstitutiveType,
             Parameters const & parameters ):
      KernelsBase( inputDofNumber,
                  inputMatrix,
                  inputRhs,
                  nodeManager,
                  elementSubRegion,
                  finiteElementSpace,
                  inputConstitutiveType,
                  parameters ),
      m_primaryField( nodeManager.template getReference< array1d< real64 > >( parameters.m_fieldName )),
      dNdX( elementSubRegion.template getReference< array3d< R1Tensor > >( dataRepository::keys::dNdX )),
      detJ( elementSubRegion.template getReference< array2d< real64 > >( dataRepository::keys::detJ ) )//,
    {}

    arrayView1d< real64 const > const m_primaryField;
    arrayView3d< R1Tensor const > const dNdX;
    arrayView2d< real64 const > const detJ;


    template< typename STACK_VARIABLE_TYPE >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void preKernel( localIndex const k,
                    STACK_VARIABLE_TYPE & stack ) const
    {
      for( localIndex a=0; a<NUM_NODES_PER_ELEM; ++a )
      {
        localIndex const localNodeIndex = elemsToNodes( k, a );

        stack.primaryField_local[ a ] = m_primaryField[ localNodeIndex ];
        stack.localRowDofIndex[a] = m_dofNumber[localNodeIndex];
        stack.localColDofIndex[a] = m_dofNumber[localNodeIndex];
      }
    }

    template< typename PARAMETERS_TYPE,
              typename STACK_VARIABLE_TYPE,
              typename DYNAMICS_LAMBDA = std::function< void( localIndex, localIndex) > >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void stiffnessKernel( localIndex const k,
                          localIndex const q,
                          PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM( parameters ),
                          STACK_VARIABLE_TYPE & stack ) const
    {
      for( localIndex a=0; a<NUM_NODES_PER_ELEM; ++a )
      {
        for( localIndex b=0; b<NUM_NODES_PER_ELEM; ++b )
        {
          stack.localJacobian[ a ][ b ] += Dot( dNdX( k, q, a ), dNdX( k, q, b ) ) * detJ( k, q );
        }
      }
    }

    template< typename PARAMETERS_TYPE, typename STACK_VARIABLE_TYPE >
    //GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    real64 postKernel( localIndex const GEOSX_UNUSED_PARAM(k),
                       PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM( parameters ),
                       STACK_VARIABLE_TYPE & stack ) const
    {
      for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
      {
        for( localIndex b = 0; b < NUM_NODES_PER_ELEM; ++b )
        {
          stack.localResidual[ a ] += stack.localJacobian[ a ][ b ] * stack.primaryField_local[ b ];
        }
      }

      m_matrix.add( stack.localRowDofIndex,
                    stack.localColDofIndex,
                    &(stack.localJacobian[0][0]),
                    stack.numRows,
                    stack.numCols );

      m_rhs.add( stack.localRowDofIndex,
                 stack.localResidual,
                 stack.numRows );

      return 1.0;
    }

  };

};

} // namespace geosx
#endif // GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_KERNELS_FEM_HPP_
