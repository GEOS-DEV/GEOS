/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file CellBlock.hpp
 */

#ifndef GEOSX_SOLIDMECHANICSLAGRANGIANSSLEKERNELCONFIG_HPP
#define GEOSX_SOLIDMECHANICSLAGRANGIANSSLEKERNELCONFIG_HPP

#define INLINE_STRESS_UPDATE 1
#define STORE_NODE_DATA_LOCALLY 0

#define SSLE_USE_PATCH_KERNEL 1

#if SSLE_USE_PATCH_KERNEL
#define SSLE_PATCH_KERNEL_MAX_ELEMS 64
#define SSLE_PATCH_KERNEL_MAX_NODES 128
#define SSLE_PATCH_KERNEL_SHARED_FVEC 0
#define SSLE_PATCH_KERNEL_REORDER_NODES 1
#define SSLE_PATCH_KERNEL_VIZ_OUTPUT 0
#endif

#endif //GEOSX_SOLIDMECHANICSLAGRANGIANSSLEKERNELCONFIG_HPP
