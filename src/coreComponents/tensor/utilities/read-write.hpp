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
 * @file read-write.hpp
 */

#ifndef GEOSX_TENSOR_READWRITE
#define GEOSX_TENSOR_READWRITE

#include "tensor/tensor.hpp"

namespace geosx
{

namespace tensor
{

/// Generate a Tensor that be read on device
template <typename C, typename L>
auto Read(const TensorBase<C,L> &t)
{
   return TensorBase<ReadContainer<T>,Layout>(t.ReadData(),t);
}

/// Generate a Tensor that be writen on device (read is unsafe)
template <typename C, typename L>
auto Write(TensorBase<C,L> &t)
{
   return TensorBase<DeviceContainer<T>,Layout>(t.WriteData(),t);
}

/// Generate a Tensor that be read and writen on device
template <typename C, typename L>
auto ReadWrite(TensorBase<C,L> &t)
{
   return TensorBase<DeviceContainer<T>,Layout>(t.ReadWriteData(),t);
}

} // namespace tensor

} // geosx namespace

#endif // GEOSX_TENSOR_READWRITE
