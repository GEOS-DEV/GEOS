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
 * @file matrix_multiplication.hpp
 */

#ifndef GEOSX_TENSOR_MATRIX_MULT
#define GEOSX_TENSOR_MATRIX_MULT

#include "common/GEOS_RAJA_Interface.hpp"
#include "tensor/tensor_traits.hpp"
#include "tensor/utilities/foreach.hpp"

namespace geosx
{

namespace tensor
{

template <typename Matrix,
          typename Vector,
          std::enable_if_t<
             get_tensor_rank<Matrix> == 2 &&
             get_tensor_rank<Vector> == 1 &&
             is_serial_tensor<Matrix> &&
             is_serial_tensor<Vector>,
             bool> = true >
GEOSX_HOST_DEVICE inline
auto operator*(const Matrix &M, const Vector &u)
{
   using Scalar = get_tensor_type<Vector>;
   constexpr int NCols = get_tensor_size<1,Matrix>;
   using Res = ResultTensor<Vector,NCols>;
   constexpr int Rows = 0;
   constexpr int Cols = 1;
   const int NCols_r = M.template Size<1>();
   Res v(NCols_r);
   Foreach<Rows>(M,[&](int row){
      Scalar res = 0;
      Foreach<Cols>(M, [&](int col){
          res += M(row,col) * u(col);
      });
      v(row) = res;
   });
}

} // namespace tensor

} // namespace geosx

#endif // GEOSX_TENSOR_MATRIX_MULT
