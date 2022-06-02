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
 * @file quadrature_operator.hpp
 */

#ifndef GEOSX_TENSOR_QUAD_OP
#define GEOSX_TENSOR_QUAD_OP

namespace geosx
{

namespace tensor
{

/// Represent an operator that goes from dofs to data at quadrature points.
template <typename QData, typename Basis>
struct QuadratureOperator
{
   const QData qdata;
   const Basis basis;

   GEOSX_HOST_DEVICE
   QuadratureOperator(const QData& qdata, const Basis& basis)
   : qdata(qdata), basis(basis)
   { }
};

template <typename QData,
          typename Basis,
          std::enable_if_t<
             is_qdata<QData> &&
             is_basis<Basis>,
             bool> = true >
GEOSX_HOST_DEVICE
auto operator*(const QData& qdata, const Basis& basis)
{
   return QuadratureOperator<QData,Basis>(qdata, basis);
}

/// Represent an operator that goes from data at quadrature points to dofs.
template <typename QData, typename Basis>
struct TransposeQuadratureOperator
{
   const QData qdata;
   const Basis basis;

   GEOSX_HOST_DEVICE
   TransposeQuadratureOperator(const QData& qdata, const Basis& basis)
   : qdata(qdata), basis(basis)
   { }
};

template <typename QData,
          typename Basis,
          std::enable_if_t<
             is_qdata<QData> &&
             is_basis<Basis>,
             bool> = true >
GEOSX_HOST_DEVICE
auto operator*(const Basis& basis, const QData& qdata)
{
   return TransposeQuadratureOperator<QData,Basis>(qdata, basis);
}

} // namespace tensor

} // namespace geosx

#endif // GEOSX_TENSOR_QUAD_OP
