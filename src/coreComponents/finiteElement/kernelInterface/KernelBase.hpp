/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
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
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

/**
 * @namespace finiteElement Contains the finite element implementation.
 */
namespace finiteElement
{

/**
 * @brief Compute the values of a finite element field at given quadrature points.
 * 
 * @param[in] basis Operator returning the shape functions values at the desired
 *                  quadrature points.
 *                  `basis ( dof, quad ) = phi_dof ( x_quad )`
 * @param[in] dofs The sets of degrees of freedom or support points.
 * @return Values of the finite element field described by @a dofs
 *         evaluated at quadrature points.
 */
template < typename Basis,
           typename Dofs,
           typename QuadValues,
           std::enable_if_t<
            !is_tensor_basis<Basis> ||
            ( is_tensor_basis<Basis> && 
              get_basis_dim<Basis> == 1 ),
            bool > = true >
GEOSX_HOST_DEVICE
auto interpolateAtQuadraturePoints( Basis const & basis,
                                    Dofs const & dofs)
{
  using T = get_value_type<QuadValues>;

  // Matrix vector product where each thread computes a value
  constexpr size_t num_quads = get_num_quads<Basis>;
  using Result = basis_result<Basis, num_quads>;

  Result q_values;

  forall<Result>([&](size_t quad)
  {
    T res{};
    forall<Dofs>([&](size_t dof)
    {
      res += basis( dof, quad ) * dofs( dof );
    });
    q_values( quad ) = res;
  });

  return q_values;
}

/**
 * @brief Compute the values of a finite element field at given quadrature points.
 * 
 * @param[in] basis Operator returning the shape functions values at the desired
 *                  quadrature points.
 *                  `basis ( dof, quad ) = phi_dof ( x_quad )`
 * @param[in] dofs The sets of degrees of freedom or support points.
 * @return Values of the finite element field described by @a dofs
 *         evaluated at quadrature points.
 */
template < typename Basis,
           typename Dofs,
           typename QuadValues,
           std::enable_if_t<
            is_tensor_basis<Basis> && 
            get_basis_dim<Basis> == 2,
            bool > = true >
GEOSX_HOST_DEVICE
auto interpolateAtQuadraturePoints( Basis const & basis,
                                    Dofs const & dofs)
{
  using T = get_quads_value_type<Basis>
  constexpr size_t num_quads = get_num_quads<Basis>;
  constexpr size_t num_dofs = get_num_dofs<Dofs>;

  // Contraction on the first dimension
  using Tmp = basis_result<Basis, num_quads, num_dofs>;
  Tmp Bu;
  
  foreach_dim<Dofs, 1>([&](size_t dof_y)
  {
    foreach_dim<Tmp, 0>([&](size_t quad_x)
    {
      T res{};
      foreach_dim<Dofs, 0>([&](size_t dof_x)
      {
        res += basis( dof_x, quad_x ) * dofs( dof_x, dof_y );
      });
      Bu( quad_x, dof_y ) = res;
    });
  });

  // Contraction on the second dimension
  using Result = basis_result<Basis, num_quads, num_quads>;
  Result q_values;

  foreach_dim<Result, 0>([&](size_t quad_x)
  {
    foreach_dim<Result, 1>([&](size_t quad_y)
    {
      T res{};
      foreach_dim<Tmp, 1>([&](size_t dof_y)
      {
        res += basis( dof_y, quad_y ) * Bu( quad_x, dof_y );
      });
      q_values( quad_x, quad_y ) = res;
    });
  });

  return q_values;
}

/**
 * @brief Compute the values of a finite element field at given quadrature points.
 * 
 * @param[in] basis Operator returning the shape functions values at the desired
 *                  quadrature points.
 *                  `basis ( dof, quad ) = phi_dof ( x_quad )`
 * @param[in] dofs The sets of degrees of freedom or support points.
 * @return Values of the finite element field described by @a dofs
 *         evaluated at quadrature points.
 */
template < typename Basis,
           typename Dofs,
           typename QuadValues,
           std::enable_if_t<
            is_tensor_basis<Basis> && 
            get_basis_dim<Basis> == 3,
            bool > = true >
GEOSX_HOST_DEVICE
auto interpolateAtQuadraturePoints( Basis const & basis,
                                    Dofs const & dofs)
{
  using T = get_quads_value_type<Basis>
  constexpr size_t num_quads = get_num_quads<Basis>;
  constexpr size_t num_dofs = get_num_dofs<Dofs>;

  // Contraction on the first dimension
  using TmpX = basis_result<Basis, num_quads, num_dofs, num_dofs>;
  TmpX Bu;

  foreach_dim<Dofs, 2>([&](size_t dof_z)
  {
    foreach_dim<Dofs, 1>([&](size_t dof_y)
    {
      foreach_dim<TmpX, 0>([&](size_t quad_x)
      {
        T res{};
        foreach_dim<Dofs, 0>([&](size_t dof_x)
        {
          res += basis( dof_x, quad_x ) * dofs( dof_x, dof_y, dof_z );
        });
        Bu( quad_x, dof_y, dof_z ) = res;
      });
    });
  });

  // Contraction on the second dimension
  using TmpY = basis_result<Basis, num_quads, num_quads, num_dofs>;
  TmpY BBu;
  
  foreach_dim<TmpY, 2>([&](size_t dof_z)
  {
    foreach_dim<TmpY, 0>([&](size_t quad_x)
    {
      foreach_dim<TmpY, 1>([&](size_t quad_y)
      {
        T res{};
        foreach_dim<TmpX, 0>([&](size_t dof_y)
        {
          res += basis( dof_y, quad_y ) * Bu( quad_x, dof_y, dof_z );
        });
        BBu( quad_x, quad_y, dof_z ) = res;
      });
    });
  });

  // Contraction on the third dimension
  using Result = basis_result<Basis, num_quads, num_quads, num_quads>;
  Result q_values;

  foreach_dim<Result, 1>([&](size_t quad_y)
  {
    foreach_dim<Result, 0>([&](size_t quad_x)
    {
      foreach_dim<Result, 2>([&](size_t quad_z)
      {
        T res{};
        foreach_dim<TmpY, 2>([&](size_t dof_z)
        {
          res += basis( dof_z, quad_z ) * BBu( quad_x, quad_y, dof_z );
        });
        q_values( quad_x, quad_y, quad_z ) = res;
      });
    });
  });

  return q_values;
}

/**
 * @brief "Apply" test functions.
 * 
 * @param[in] basis Operator returning the shape functions values at the desired
 *                  quadrature points.
 *                  `basis ( dof, quad ) = phi_dof ( x_quad )`
 * @param[in] q_values Values of the "qFunction" at quadrature points.
 * @return Contribution of the q_values to the degrees of freedom.
 */
template < typename Basis,
           typename QValues >
GEOSX_HOST_DEVICE
auto applyTestFunctions( Basis const & basis,
                         QValues const & q_values )
{
  interpolateAtQuadraturePoints( transpose( basis ), q_values );
}

/**
 * @brief Compute the gradient values of a finite element field at given quadrature points.
 * 
 * @param[in] basis Operator returning the shape functions values at the desired
 *                  quadrature points.
 *                  `basis ( dof, quad ) = phi_dof ( x_quad )`
 * @param[in] dofs The sets of degrees of freedom or support points.
 * @return Gradient values of the finite element field described by @a dofs
 *         evaluated at quadrature points.
 */
template < typename Basis,
           typename Dofs,
           typename QuadValues,
           std::enable_if_t<
            is_tensor_basis<Basis> && 
            get_basis_dim<Basis> == 3,
            bool > = true >
GEOSX_HOST_DEVICE
auto interpolateGradientAtQuadraturePoints( Basis const & basis,
                                            Dofs const & u)
{
  using T = get_quads_value_type<Basis>
  constexpr size_t num_quads = get_num_quads<Basis>;
  constexpr size_t num_dofs = get_num_dofs<Dofs>;

  // Contraction on the first dimension
  using TmpX = basis_result<Basis, num_quads, num_dofs, num_dofs>;
  TmpX Bu, Gu;

  foreach_dim<Dofs, 2>([&](size_t dof_z)
  {
    foreach_dim<Dofs, 1>([&](size_t dof_y)
    {
      foreach_dim<TmpX, 0>([&](size_t quad_x)
      {
        T bu{};
        T gu{};
        foreach_dim<Dofs, 0>([&](size_t dof_x)
        {
          const T val = dofs( dof_x, dof_y, dof_z );
          bu += basis( dof_x, quad_x ) * val;
          gu += gradient( basis )( dof_x, quad_x ) * val;
        });
        Bu( quad_x, dof_y, dof_z ) = bu;
        Gu( quad_x, dof_y, dof_z ) = gu;
      });
    });
  });

  // Contraction on the second dimension
  using TmpY = basis_result<Basis, num_quads, num_quads, num_dofs>;
  TmpY BBu, BGu, GBu;
  
  foreach_dim<TmpY, 2>([&](size_t dof_z)
  {
    foreach_dim<TmpY, 0>([&](size_t quad_x)
    {
      foreach_dim<TmpY, 1>([&](size_t quad_y)
      {
        T bbu{};
        T bgu{};
        T gbu{};
        foreach_dim<TmpX, 0>([&](size_t dof_y)
        {
          const T bu = Bu( quad_x, dof_y, dof_z );
          const T gu = Gu( quad_x, dof_y, dof_z );
          bbu += basis( dof_y, quad_y ) * bu;
          gbu += gradient( basis )( dof_y, quad_y ) * bu;
          bgu += basis( dof_y, quad_y ) * gu;
        });
        BBu( quad_x, quad_y, dof_z ) = bbu;
        GBu( quad_x, quad_y, dof_z ) = gbu;
        BGu( quad_x, quad_y, dof_z ) = bgu;
      });
    });
  });

  // Contraction on the third dimension
  using Result = basis_result<Basis, num_quads, num_quads, num_quads, 3>; // remove magic number?
  Result q_values;

  foreach_dim<Result, 1>([&](size_t quad_y)
  {
    foreach_dim<Result, 0>([&](size_t quad_x)
    {
      foreach_dim<Result, 2>([&](size_t quad_z)
      {
        T gbbu{};
        T bgbu{};
        T bbgu{};
        foreach_dim<TmpY, 2>([&](size_t dof_z)
        {
          const T bbu = BBu( quad_x, quad_y, dof_z );
          const T gbu = GBu( quad_x, quad_y, dof_z );
          const T bgu = BGu( quad_x, quad_y, dof_z );
          gbbu += gradient( basis )( dof_z, quad_z ) * bbu;
          bgbu += basis( dof_z, quad_z ) * gbu;
          bbgu += basis( dof_z, quad_z ) * bgu;
        });
        q_values( quad_x, quad_y, quad_z, 0 ) = bbgu;
        q_values( quad_x, quad_y, quad_z, 1 ) = bgbu;
        q_values( quad_x, quad_y, quad_z, 2 ) = gbbu;
      });
    });
  });

  return q_values;
}

/**
 * @brief "Apply" gradient of the test functions.
 * 
 * @param[in] basis Operator returning the shape functions values at the desired
 *                  quadrature points.
 *                  `basis ( dof, quad ) = phi_dof ( x_quad )`
 * @param[out] q_values Values of the "qFunction" at quadrature points.
 * @return Contribution of the q_values to the degrees of freedom.
 */
template < typename Basis,
           typename Qvalues,
           std::enable_if_t<
            is_tensor_basis<Basis> && 
            get_basis_dim<Basis> == 3,
            bool > = true >
auto applyGradientTestFunctions( Basis const & basis,
                                 Qvalues const & q_values )
{
  using T = get_quads_value_type<Basis>
  constexpr size_t num_quads = get_num_quads<Basis>;
  constexpr size_t num_dofs = get_num_dofs<Dofs>;

  // Contraction on the first dimension
  using TmpX = basis_result<Basis, num_dofs, num_quads, num_quads>;
  TmpX Gqx, Bqy, Bqz;

  foreach_dim<QValues, 2>([&](size_t quad_z)
  {
    foreach_dim<QValues, 1>([&](size_t quad_y)
    {
      foreach_dim<TmpX, 0>([&](size_t dof_x)
      {
        T gqx{};
        T bqy{};
        T bqz{};
        foreach_dim<QValues, 0>([&](size_t quad_x)
        {
          // Using gradient at quadrature points prevents us from storing so many tmp while also reducing FLOPs
          const T qx = q_values( quad_x, quad_y, quad_z, 0 );
          const T qy = q_values( quad_x, quad_y, quad_z, 1 );
          const T qz = q_values( quad_x, quad_y, quad_z, 2 );
          const T b = basis( dof_x, quad_x );
          const T g = gradient( basis )( dof_x, quad_x );
          gqx += g * qx;
          bqy += b * qy;
          bqz += b * qz;
        });
        Gqx( dof_x, quad_y, quad_z ) = gqx;
        Bqy( dof_x, quad_y, quad_z ) = bqy;
        Bqz( dof_x, quad_y, quad_z ) = bqz;
      });
    });
  });

  // Contraction on the second dimension
  using TmpY = basis_result<Basis, num_dofs, num_dofs, num_quads>;
  TmpY BGqx, GBqy, BBqz;
  
  foreach_dim<TmpY, 2>([&](size_t quad_z)
  {
    foreach_dim<TmpY, 0>([&](size_t dof_x)
    {
      foreach_dim<TmpY, 1>([&](size_t dof_y)
      {
        T bgqx{};
        T gbqy{};
        T bbqz{};
        foreach_dim<TmpX, 0>([&](size_t quad_y)
        {
          const T gqx = Gqx( dof_x, quad_y, quad_z );
          const T bqy = Bqy( dof_x, quad_y, quad_z );
          const T bqz = Bqz( dof_x, quad_y, quad_z );
          const T b = basis( dof_y, quad_y );
          const T g = gradient( basis )( dof_y, quad_y );
          bgqx += b * gqx;
          gbqy += g * bqy;
          bbqz += b * bqz;
        });
        BGqx( dof_x, dof_y, quad_z ) = bgqx;
        GBqy( dof_x, dof_y, quad_z ) = gbqy;
        BBqz( dof_x, dof_y, quad_z ) = bbqz;
      });
    });
  });

  // Contraction on the third dimension
  using Result = basis_result<Basis, num_dofs, num_dofs, num_dofs>;
  Result dofs;

  foreach_dim<Result, 1>([&](size_t dof_y)
  {
    foreach_dim<Result, 0>([&](size_t dof_x)
    {
      foreach_dim<Result, 2>([&](size_t dof_z)
      {
        T res{};
        foreach_dim<TmpY, 2>([&](size_t quad_z)
        {
          const T bgqx = BGqx( dof_x, dof_y, quad_z );
          const T gbqy = GBqy( dof_x, dof_y, quad_z );;
          const T bbqz = BBqz( dof_x, dof_y, quad_z );;
          const T b = basis( dof_z, quad_z );
          const T g = gradient( basis )( dof_z, quad_z )
          res += b * bgqx + b * gbqy + g * bbqz;
        });
        dofs( dof_x, dof_y, dof_z ) = res;
      });
    });
  });

  return dofs;  
}

/** @brief Basic non-tensor basis that stores its values. */
template < typename FiniteElement >
class StoredNonTensorBasis
{
private:
  constexpr localIndex num_dofs = get_num_dofs<FiniteElement>;
  constexpr localIndex num_quads = get_num_quads<FiniteElement>;
  real64 const data[ num_quads ][ num_dofs ];

public:
  GEOSX_HOST_DEVICE
  StoredNonTensorBasis()
  {
    real64 basis_functions_at_xq[ num_dofs ];
    for (localIndex q = 0; q < num_quads; q++)
    {
      FiniteElement::calcN( q, basis_functions_at_xq );
      for (localIndex d = 0; d < num_dofs; d++)
      {
        data[ q ][ d ] = basis_functions_at_xq[ d ];
      }
    }
  }

  GEOSX_HOST_DEVICE
  real64 operator()( localIndex dof, localIndex quad ) const
  {
    return data[ quad ][ dof ];
  }
};

template < typename Basis, typename MeshBasis, typename Dofs, typename MeshDofs>
void laplaceKernel(Basis const & basis, MeshBasis const & mesh_basis, Dofs const &u, MeshDofs const x, Dofs & v)
{
  for ( size_t k = BlockIdx.x ; k < num_elems; k+=BlockDim.x )
  {
    auto local_dofs = get_local_elem_dofs(k, u);
    auto mesh_local_dofs = get_local_elem_dofs(k, x);
    auto grad_u_q = interpolateGradientAtQuadraturePoints( basis, u );
    auto J = interpolateGradientAtQuadraturePoints( mesh_basis, x );
    auto q_function = [](auto weight, auto J, auto grad_u_q )
    {
      auto inv_J = inverse( J );
      return weight * det( inv_J ) * inv_J * transpose( inv_J ) * u_q * square ( grad_u_q );
    };
    auto d = apply_q_function(q_function, basis.weight, J, grad_u_q);
    // serial apply_q_function code:
    // for (size_t q = 0; q < num_q_pts; q++)
    // {
    //   auto inv_J = inverse( J(q) );
    //   d(q) = basis.weight(q) * det( inv_J ) * inv_J * transpose( inv_J ) * grad_u_q( q );
    // }
    auto local_v = applyGradientTestFunctions(basis, d);
    gather(v, k, local_v);
  }
}

template < size_t Order >
class LagrangeBasis;

template <>
class LagrangeBasis<1> : public LagrangeBasis1 { };

template < typename FiniteElement >
class StoredLagrangeTensorBasis
{
private:
  constexpr localIndex num_dofs = get_num_1d_dofs<FiniteElement>;
  constexpr localIndex num_quads = get_num_1d_quads<FiniteElement>;
  real64 const data[ num_quads ][ num_dofs ];

public:
  GEOSX_HOST_DEVICE
  constexpr StoredLagrangeTensorBasis()
  {
    for (localIndex q = 0; q < num_quads; q++)
    {
      for (localIndex d = 0; d < num_dofs; d++)
      {
        data[ q ][ d ] = LagrangeBasis<num_dofs-1>::value( d, LagrangeBasis<num_quads-1>::parentSupportCoord( q ) );
      }
    }
  }

  GEOSX_HOST_DEVICE
  constexpr real64 operator()( localIndex dof, localIndex quad ) const
  {
    return data[ quad ][ dof ];
  }
};

template < typename FiniteElement >
class OnTheFlyLagrangeTensorBasis
{
private:
  constexpr localIndex num_dofs = get_num_1d_dofs<FiniteElement>;
  constexpr localIndex num_quads = get_num_1d_quads<FiniteElement>;

public:
  GEOSX_HOST_DEVICE
  constexpr OnTheFlyLagrangeTensorBasis() { }

  GEOSX_HOST_DEVICE
  constexpr real64 operator()( localIndex dof, localIndex quad ) const
  {
    return LagrangeBasis<num_dofs-1>::value( dof, LagrangeBasis<num_quads-1>::parentSupportCoord( quad ) );
  }
};

/** @brief A basic tensor class statically sized. */
template < typename T, size_t... Dims >
StaticTensor
{
private:
  T data[ 1 * ... * Dims ];

public:
  template < typename... Args >
  T const & operator()( Args... indices )
  {
    static_assert( sizeof...(Args) == sizeof...(Dims),
                   "Wrong number of arguments" );
    return data[ 0 ]; // TODO
  }
};

template < typename Tensor,
           size_t Dim,
           typename LambdaFn,
           std::enable_if_t<
             is_serial_tensor<Tensor>,
             bool > = true >
GEOSX_HOST_DEVICE
void foreach_dim( LambdaFn && fn )
{
  for (int i = 0; i < get_tensor_size<Tensor, Dim>; i++)
  {
    fn(i);
  }
}

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
          typename FE_TYPE,
          int NUM_DOF_PER_TEST_SP,
          int NUM_DOF_PER_TRIAL_SP >
class KernelBase
{
public:
  /// Compile time value for the number of test function support points per
  /// element.
  static constexpr int maxNumTestSupportPointsPerElem  = FE_TYPE::maxSupportPoints;

  /// Compile time value for the number of trial function support points per
  /// element.
  static constexpr int maxNumTrialSupportPointsPerElem = FE_TYPE::maxSupportPoints;

  /// Compile time value for the number of degrees of freedom per test function
  /// support point.
  static constexpr int numDofPerTestSupportPoint    = NUM_DOF_PER_TEST_SP;

  /// Compile time value for the number of degrees of freedom per trial
  /// function support point.
  static constexpr int numDofPerTrialSupportPoint   = NUM_DOF_PER_TRIAL_SP;

  /// Compile time value for the number of quadrature points per element.
  static constexpr int numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

  /**
   * @brief Constructor
   * @param elementSubRegion Reference to the SUBREGION_TYPE(class template
   *                         parameter) object.
   * @param inputConstitutiveType The constitutive object.
   * @param finiteElementSpace Placeholder for the finite element space object,
   *                           which currently doesn't do much.
   */
  KernelBase( SUBREGION_TYPE const & elementSubRegion,
              FE_TYPE const & finiteElementSpace,
              CONSTITUTIVE_TYPE & inputConstitutiveType ):
    m_elemsToNodes( elementSubRegion.nodeList().toViewConst() ),
    m_elemGhostRank( elementSubRegion.ghostRank() ),
    m_constitutiveUpdate( inputConstitutiveType.createKernelUpdates() ),
    m_finiteElementSpace( finiteElementSpace )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables allocated on the stack.
   *
   * ### ImplicitKernelBase::StackVariables Description
   *
   * Contains variables that will be allocated on the stack of the main kernel.
   * This will typically consist of local arrays to hold data mapped from the
   * global data arrays, and/or local storage for the residual and jacobian
   * contributions.
   */
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
   * ### KernelBase::quadraturePointKernel() Description
   *
   * The operations found here are the mapping from the support points to the
   * quadrature point, calculation of gradients, etc. From this data the
   * state of the constitutive model is updated if required by the physics
   * package.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
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
   * @return The maximum contribution to the residual.
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
   * @return The maximum residual contribution.
   *
   * This is a generic launching function for all of the finite element kernels
   * that follow the interface set by KernelBase.
   */
  //START_kernelLauncher
  template< typename POLICY,
            typename KERNEL_TYPE >
  static
  real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    GEOSX_MARK_FUNCTION;

    // Define a RAJA reduction variable to get the maximum residual contribution.
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxResidual( 0 );

    forAll< POLICY >( numElems,
                      [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( k, stack );
      for( integer q=0; q<numQuadraturePointsPerElem; ++q )
      {
        kernelComponent.quadraturePointKernel( k, q, stack );
      }
      maxResidual.max( kernelComponent.complete( k, stack ) );
    } );
    return maxResidual.get();
  }
  //END_kernelLauncher

protected:
  /// The element to nodes map.
  traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType::base_type > const m_elemsToNodes;

  /// The element ghost rank array.
  arrayView1d< integer const > const m_elemGhostRank;

  /// The constitutive update object used to update the constitutive state,
  /// and extract constitutive data.
  typename CONSTITUTIVE_TYPE::KernelWrapper const m_constitutiveUpdate;

  /// The finite element space/discretization object for the element type in
  /// the SUBREGION_TYPE.
  FE_TYPE const & m_finiteElementSpace;
};

/**
 * @class KernelFactory
 * @brief Used to forward arguments to a class that implements the KernelBase interface.
 * @tparam KERNEL_TYPE The template class to construct, should implement the KernelBase interface.
 * @tparam ARGS The arguments used to construct a @p KERNEL_TYPE in addition to the standard arguments.
 */
template< template< typename SUBREGION_TYPE,
                    typename CONSTITUTIVE_TYPE,
                    typename FE_TYPE > class KERNEL_TYPE,
          typename ... ARGS >
class KernelFactory
{
public:

  /**
   * @brief Initialize the factory.
   * @param args The arguments used to construct a @p KERNEL_TYPE in addition to the standard arguments.
   */
  KernelFactory( ARGS ... args ):
    m_args( args ... )
  {}

  /**
   * @brief Create a new kernel with the given standard arguments.
   * @tparam SUBREGION_TYPE The type of @p elementSubRegion.
   * @tparam CONSTITUTIVE_TYPE The type of @p inputConstitutiveType.
   * @tparam FE_TYPE The type of @p finiteElementSpace.
   * @param nodeManager The node manager.
   * @param edgeManager The edge manager.
   * @param faceManager The face manager.
   * @param targetRegionIndex The target region index.
   * @param elementSubRegion The subregion to execute on.
   * @param finiteElementSpace The finite element space.
   * @param inputConstitutiveType The constitutive relation.
   * @return A new kernel constructed with the given arguments and @c ARGS.
   */
  template< typename SUBREGION_TYPE, typename CONSTITUTIVE_TYPE, typename FE_TYPE >
  KERNEL_TYPE< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE > createKernel(
    NodeManager & nodeManager,
    EdgeManager const & edgeManager,
    FaceManager const & faceManager,
    localIndex const targetRegionIndex,
    SUBREGION_TYPE const & elementSubRegion,
    FE_TYPE const & finiteElementSpace,
    CONSTITUTIVE_TYPE & inputConstitutiveType )
  {
    camp::tuple< NodeManager &,
                 EdgeManager const &,
                 FaceManager const &,
                 localIndex const,
                 SUBREGION_TYPE const &,
                 FE_TYPE const &,
                 CONSTITUTIVE_TYPE & > standardArgs { nodeManager,
                                                      edgeManager,
                                                      faceManager,
                                                      targetRegionIndex,
                                                      elementSubRegion,
                                                      finiteElementSpace,
                                                      inputConstitutiveType };

    auto allArgs = camp::tuple_cat_pair( standardArgs, m_args );
    return camp::make_from_tuple< KERNEL_TYPE< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE > >( allArgs );
  }

private:
  /// The arguments to append to the standard kernel constructor arguments.
  camp::tuple< ARGS ... > m_args;
};


//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

//START_regionBasedKernelApplication
/**
 * @brief Performs a loop over specific regions (by type and name) and calls a kernel launch on the subregions
 *   with compile time knowledge of sub-loop bounds such as number of nodes and quadrature points per element.
 * @tparam POLICY The RAJA launch policy to pass to the kernel launch.
 * @tparam CONSTITUTIVE_BASE The common base class for constitutive pass-thru/dispatch which gives the kernel
 *   launch compile time knowledge of the constitutive model. This is achieved through a call to the
 *   ConstitutivePassThru function which should have a specialization for CONSTITUTIVE_BASE implemented in
 *   order to perform the compile time dispatch.
 * @tparam SUBREGION_TYPE The type of subregion to loop over. TODO make this a parameter pack?
 * @tparam KERNEL_FACTORY The type of @p kernelFactory, typically an instantiation of @c KernelFactory, and
 *   must adhere to that interface.
 * @param mesh The MeshLevel object.
 * @param targetRegions The names of the target regions(of type @p SUBREGION_TYPE) to apply the @p KERNEL_TEMPLATE.
 * @param finiteElementName The name of the finite element.
 * @param constitutiveStringName The key to the constitutive model name found on the Region.
 * @param kernelFactory The object used to construct the kernel.
 * @return The maximum contribution to the residual, which may be used to scale the residual.
 *
 * @details Loops over all regions Applies/Launches a kernel specified by the @p KERNEL_TEMPLATE through
 * #::geosx::finiteElement::KernelBase::kernelLaunch().
 */
template< typename POLICY,
          typename CONSTITUTIVE_BASE,
          typename SUBREGION_TYPE,
          typename KERNEL_FACTORY >
static
real64 regionBasedKernelApplication( MeshLevel & mesh,
                                     arrayView1d< string const > const & targetRegions,
                                     string const & finiteElementName,
                                     string const & constitutiveStringName,
                                     KERNEL_FACTORY & kernelFactory )
{
  GEOSX_MARK_FUNCTION;
  // save the maximum residual contribution for scaling residuals for convergence criteria.
  real64 maxResidualContribution = 0;

  NodeManager & nodeManager = mesh.getNodeManager();
  EdgeManager & edgeManager = mesh.getEdgeManager();
  FaceManager & faceManager = mesh.getFaceManager();
  ElementRegionManager & elementRegionManager = mesh.getElemManager();

  // Loop over all sub-regions in regions of type SUBREGION_TYPE, that are listed in the targetRegions array.
  elementRegionManager.forElementSubRegions< SUBREGION_TYPE >( targetRegions,
                                                               [&constitutiveStringName,
                                                                &maxResidualContribution,
                                                                &nodeManager,
                                                                &edgeManager,
                                                                &faceManager,
                                                                &kernelFactory,
                                                                &finiteElementName]
                                                                 ( localIndex const targetRegionIndex, auto & elementSubRegion )
  {
    localIndex const numElems = elementSubRegion.size();

    // Get the constitutive model...and allocate a null constitutive model if required.

    constitutive::ConstitutiveBase * constitutiveRelation = nullptr;
    constitutive::NullModel * nullConstitutiveModel = nullptr;
    if( elementSubRegion.template hasWrapper< string >( constitutiveStringName ) )
    {
      string const & constitutiveName = elementSubRegion.template getReference< string >( constitutiveStringName );
      constitutiveRelation = &elementSubRegion.template getConstitutiveModel( constitutiveName );
    }
    else
    {
      nullConstitutiveModel = &elementSubRegion.template registerGroup< constitutive::NullModel >( "nullModelGroup" );
      constitutiveRelation = nullConstitutiveModel;
    }

    // Call the constitutive dispatch which converts the type of constitutive model into a compile time constant.
    constitutive::ConstitutivePassThru< CONSTITUTIVE_BASE >::execute( *constitutiveRelation,
                                                                      [&maxResidualContribution,
                                                                       &nodeManager,
                                                                       &edgeManager,
                                                                       &faceManager,
                                                                       targetRegionIndex,
                                                                       &kernelFactory,
                                                                       &elementSubRegion,
                                                                       &finiteElementName,
                                                                       numElems]
                                                                        ( auto & castedConstitutiveRelation )
    {
      FiniteElementBase &
      subRegionFE = elementSubRegion.template getReference< FiniteElementBase >( finiteElementName );

      finiteElement::dispatch3D( subRegionFE,
                                 [&maxResidualContribution,
                                  &nodeManager,
                                  &edgeManager,
                                  &faceManager,
                                  targetRegionIndex,
                                  &kernelFactory,
                                  &elementSubRegion,
                                  numElems,
                                  &castedConstitutiveRelation] ( auto const finiteElement )
      {
        auto kernel = kernelFactory.createKernel( nodeManager,
                                                  edgeManager,
                                                  faceManager,
                                                  targetRegionIndex,
                                                  elementSubRegion,
                                                  finiteElement,
                                                  castedConstitutiveRelation );

        using KERNEL_TYPE = decltype( kernel );

        // Call the kernelLaunch function, and store the maximum contribution to the residual.
        maxResidualContribution =
          std::max( maxResidualContribution,
                    KERNEL_TYPE::template kernelLaunch< POLICY, KERNEL_TYPE >( numElems, kernel ) );
      } );
    } );

    // Remove the null constitutive model (not required, but cleaner)
    if( nullConstitutiveModel )
    {
      elementSubRegion.deregisterGroup( "nullModelGroup" );
    }

  } );

  return maxResidualContribution;
}
//END_regionBasedKernelApplication

} // namespace finiteElement
} // namespace geosx



#endif /* GEOSX_FINITEELEMENT_KERNELBASE_HPP_ */
