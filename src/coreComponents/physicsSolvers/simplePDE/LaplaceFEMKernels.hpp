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

#include "../PhysicsLoopInterface.hpp"

namespace geosx
{

class LaplaceFEMKernel
{
public:
  using Base = physicsLoopInterface::FiniteElementRegionLoop;

  struct Parameters : public Base::Parameters
  {
    Parameters( string const & fieldName ):
      Base::Parameters(),
      m_fieldName( fieldName )
    {}

    string m_fieldName;
  };


  template< int NUM_NODES_PER_ELEM, int NUM_DOF_PER_NODE >
  struct StackVariables : Base::StackVariables< NUM_NODES_PER_ELEM, NUM_DOF_PER_NODE >
  {
  public:
    using StackVariablesBase = Base::StackVariables< NUM_NODES_PER_ELEM, NUM_DOF_PER_NODE >;

    using StackVariablesBase::numNodesPerElem;
    using StackVariablesBase::numDofPerNode;
    using StackVariablesBase::ndof;


//      GEOSX_HOST_DEVICE
    StackVariables():
      StackVariablesBase(),
      primaryField_local{ 0.0 }
    {}

    real64 primaryField_local[NUM_NODES_PER_ELEM];
  };


  template< typename SUBREGION_TYPE,
            typename CONSTITUTIVE_TYPE >
  class Kernels : public Base::Kernels< SUBREGION_TYPE, CONSTITUTIVE_TYPE>
  {
  public:
    using KernelBase = Base::Kernels< SUBREGION_TYPE, CONSTITUTIVE_TYPE>;
    using KernelBase::m_dofNumber;
    using KernelBase::m_matrix;
    using KernelBase::m_rhs;
    using KernelBase::elemsToNodes;
    using KernelBase::constitutiveUpdate;
    using KernelBase::m_finiteElementSpace;




    Kernels( arrayView1d< globalIndex const > const & inputDofNumber,
             ParallelMatrix & inputMatrix,
             ParallelVector & inputRhs,
             NodeManager const & nodeManager,
             SUBREGION_TYPE const & elementSubRegion,
             FiniteElementBase const * const finiteElementSpace,
             CONSTITUTIVE_TYPE & inputConstitutiveType,
             Parameters const & parameters ):
      KernelBase( inputDofNumber,
                  inputMatrix,
                  inputRhs,
                  nodeManager,
                  elementSubRegion,
                  finiteElementSpace,
                  inputConstitutiveType,
                  inputConstitutiveType.createKernelWrapper() ),
      m_primaryField( nodeManager.template getReference< array1d< real64 > >( parameters.m_fieldName )),
      dNdX(elementSubRegion.template getReference< array3d< R1Tensor > >( dataRepository::keys::dNdX )),
      detJ(elementSubRegion.template getReference< array2d< real64 > >( dataRepository::keys::detJ ) )//,
    {}

    arrayView1d< real64 const > const m_primaryField;

    arrayView3d< R1Tensor const > const dNdX;
    arrayView2d< real64 const > const detJ;


    template< typename STACK_VARIABLE_TYPE >
  //    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void preKernel( localIndex const k,
                    STACK_VARIABLE_TYPE & stack ) const
    {
      for( localIndex a=0; a<STACK_VARIABLE_TYPE::numNodesPerElem; ++a )
      {
        localIndex const localNodeIndex = elemsToNodes( k, a );

        stack.primaryField_local[ a ] = m_primaryField[ localNodeIndex ];

        for( int i=0; i<3; ++i )
        {
          stack.elementLocalDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
        }
      }
    }

    template< typename PARAMETERS_TYPE,
              typename STACK_VARIABLE_TYPE,
              typename DYNAMICS_LAMBDA = std::function<void( localIndex, localIndex)> >
  //    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void stiffnessKernel( localIndex const k,
                          localIndex const q,
                          PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM(parameters),
                          STACK_VARIABLE_TYPE & stack ) const
    {
      for( localIndex a=0; a<STACK_VARIABLE_TYPE::numNodesPerElem; ++a )
      {
        for( localIndex b=0; b<STACK_VARIABLE_TYPE::numNodesPerElem; ++b )
        {
          stack.localJacobian[ a ][ b ] += Dot( dNdX( k, q, a ), dNdX( k, q, b ) ) * detJ( k, q );
        }
      }
    }

    template< typename PARAMETERS_TYPE, typename STACK_VARIABLE_TYPE >
  //    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    real64 postKernel( PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM( parameters ),
                       STACK_VARIABLE_TYPE & stack ) const
    {
      for( localIndex a = 0; a < STACK_VARIABLE_TYPE::numNodesPerElem; ++a )
      {
        for( localIndex b = 0; b < STACK_VARIABLE_TYPE::numNodesPerElem; ++b )
        {
          stack.localResidual[ a ] += stack.localJacobian[ a ][ b ] * stack.primaryField_local[ b ];
        }
      }

      m_matrix.add( stack.elementLocalDofIndex,
                    stack.elementLocalDofIndex,
                    &(stack.localJacobian[0][0]),
                    stack.ndof,
                    stack.ndof );

      m_rhs.add( stack.elementLocalDofIndex,
                 stack.localResidual,
                 stack.ndof );

      return 1.0;
    }

  };

};

} // namespace geosx
#endif // GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_KERNELS_FEM_HPP_
