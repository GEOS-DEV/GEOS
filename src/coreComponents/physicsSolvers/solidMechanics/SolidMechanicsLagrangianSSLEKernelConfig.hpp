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

#ifndef GEOSX_SOLIDMECHANICSLAGRANGIANSSLEKERNELCONFIG_HPP
#define GEOSX_SOLIDMECHANICSLAGRANGIANSSLEKERNELCONFIG_HPP

#define USE_ELEM_PATCHES 1

#if USE_ELEM_PATCHES

  #define ELEM_PATCH_MAX_ELEM 64
  #define ELEM_PATCH_MAX_NODE 128
  #define ELEM_PATCH_VIZ 0
  #define ELEM_PATCH_REORDER_NODES 1

  #define SSLE_USE_PATCH_KERNEL 1
  #if SSLE_USE_PATCH_KERNEL
    #define SSLE_PATCH_KERNEL_SHARED_FVEC 0
  #endif

#endif

#endif //GEOSX_SOLIDMECHANICSLAGRANGIANSSLEKERNELCONFIG_HPP
