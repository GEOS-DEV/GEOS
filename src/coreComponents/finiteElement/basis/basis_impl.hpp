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
 * @file basis_impl.hpp
 */

#ifndef GEOSX_BASIS_IMPL
#define GEOSX_BASIS_IMPL

#include "tensor/tensor_types.hpp"
#include "tensor/utilities/pow.hpp"

namespace geosx
{

namespace tensor
{

// ALL THIS SHOULD BE REWRITTEN...
// TODO Maybe remove this class?
// TODO maybe D before Q?
template <int Dim, bool IsTensor, typename TensorType>
class BasisTensor : public TensorType
{
public:
   GEOSX_HOST_DEVICE
   BasisTensor(int quads, int dofs): TensorType(quads,dofs) { }

   GEOSX_HOST_DEVICE
   BasisTensor(double *shared_mem, int quads, int dofs)
      : TensorType(shared_mem,quads,dofs) { }
};

/// Represent the rank 2 tensor containing B1d or G1d with dynamic sizes
template <int Dim>
using DynamicBasisTensor = BasisTensor<Dim,true,DynamicDTensor<2>>;
template <int Dim>
using DynamicSharedBasisTensor = BasisTensor<Dim,true,DeviceDTensor<2>>;

/// Represent the rank 2 tensor containing B1d or G1d with static sizes
template <int Dim, int Q, int D>
using StaticBasisTensor = BasisTensor<Dim,true,StaticDTensor<Q,D>>;
template <int Dim, int Q, int D>
using StaticSharedBasisTensor = BasisTensor<Dim,true,StaticPointerDTensor<Q,D>>;

/// Represent the rank 2 tensor containing B or G with dynamic sizes
template <int Dim>
using DynamicBasisNonTensor = BasisTensor<
   Dim, false, DynamicDTensor<2,pow(DynamicMaxSize,2*Dim)>>;

/// Represent the rank 2 tensor containing B or G with static sizes
template <int Dim, int Q, int D>
using StaticBasisNonTensor = BasisTensor<Dim,false,StaticDTensor<Q,D>>;

} // tensor namespace

} // geosx namespace

#endif // GEOSX_BASIS_IMPL
