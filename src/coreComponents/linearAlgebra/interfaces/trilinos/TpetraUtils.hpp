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
 * @file TpetraUtils.hpp
 */
#ifndef GEOSX_TPETRAUTILS_HPP
#define GEOSX_TPETRAUTILS_HPP

#include "common/GeosxConfig.hpp"

#ifdef GEOSX_USE_MPI
#include <Teuchos_DefaultMpiComm.hpp>

/// Alias for specific Tpetra::Comm implementation used.
using Tpetra_Comm = Teuchos::MpiComm< int >;
#else
#include <Teuchos_DefaultSerialComm.hpp>

/// Alias for specific Tpetra::Comm implementation used.
using Tpetra_Comm = Teuchos::SerialComm< int >;
#endif

#endif //GEOSX_TPETRAUTILS_HPP
