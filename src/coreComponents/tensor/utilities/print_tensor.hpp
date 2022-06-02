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
 * @file print_tensor.hpp
 */

#ifndef GEOSX_TENSOR_PRINT
#define GEOSX_TENSOR_PRINT

#include "tensor/tensor.hpp"
#include <iostream>

namespace geosx
{

namespace tensor
{

/// Function to print tensors
template <typename C, typename L>
std::ostream& operator<<(std::ostream &os, const TensorBase<C,L> &t)
{
   ForallDims<Tensor>::Apply(t,[&](auto... idx)
   {
      os << "value(" << idx... << ")= " << t(idx...) << ", ";
   });
   os << std::endl;
   return os;
}

} // namespace tensor

} // geosx namespace

#endif // GEOSX_TENSOR_PRINT
